#!/usr/bin/env python
# Writing sequence into file per bead. Quicker by making index file
#
#-----------------------------------------------------------------------------
# Author : Zhuobing Peng
# Email  : pengzhuobing@genomics.cn
# Update : March 2021 (Since Feb 2021)
#-----------------------------------------------------------------------------

import sys
import os
import time
import gzip
import glob
import shutil
import logging
import argparse
from multiprocessing import Pool

if os.name != 'nt':
    # import the file lock package if Linux system
    import fcntl


def get_parser():
    parser = argparse.ArgumentParser(description="task generation")
    parser.add_argument("--r1", help="read1 fq file", required=True)
    parser.add_argument("--r2", help="read2 fq file", required=True)
    parser.add_argument("-o", dest="out", help="output dir", required=True)
    parser.add_argument("-d", dest="IDlst", help="individual barcode id file (make sure the id sorted)", required=True)
    parser.add_argument("-t", dest="thread", help="threads to use for compressing", type=int, default=8)
    parser.add_argument("-f", dest="fmt", help="format [fq|fa|fq+fa|fq.gz]",choices=["fq", "fa", "fq.gz"], default="fq")
    parser.add_argument("-m", dest="mode", help="Assembly mode [megahit | spades | idba]", default="megahit")
    parser.add_argument("-c", dest="cpu", help="RCA cpu", type=int, default=1)
    parser.add_argument("--mem", help="RCA mem", type=int, default=10000000000)
    parser.add_argument("--mAsm", help="Assembly mem", type=int, default=500)

    return parser



parser = get_parser()
args = parser.parse_args()
in1 = args.r1
in2 = args.r2
IDlst = args.IDlst
thread = args.thread
sdir = args.out
fmt = args.fmt
mode = args.mode
cpu = args.cpu
mem = args.mem
mAsm = args.mAsm


tdir = "/dev/shm"
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(filename)s[line:%(lineno)d] - %(levelname)s: %(message)s')
logger = logging.getLogger("Assemble_RCA_Task")

result1 = os.path.join(sdir, "summary.BI.%s.clip.all.fasta" % mode)
result2 = os.path.join(sdir, "summary.BI.%s.clip.metadata.tsv" % mode)
result3 = os.path.join(sdir, "summary.BI.%s.contig.tsv" % mode)
result4 = os.path.join(sdir, "summary.BI.%s.contig.fasta" % mode)
outn = "BI.{amode}.contig.fasta".format(amode = mode)

cmp_template_map = {
    "TEST": '''
mkdir %(outdir)s/%(amode)s/
# assemble
echo metabbq RCAasm.sh %(amode)s %(sam)s BI %(Bar)s BI %(mem)s %(cpu)s %(sam)s/primers.cfg > %(outdir)s//%(amode)s/final.contigs.fa 
cat %(outdir)s//%(amode)s/final.contigs.fa > %(outdir)s/%(amode)s/%(out)s
    ''',
    "contig.fasta": '''
# assemble
metabbq RCAasm.sh %(amode)s /dev/shm BI %(Bar)s BI %(mem)s %(cpu)s %(sam)s/primers.cfg
# RCA
metabbq clusterHelper asmPick -l %(mAsm)s -b %(Bar)s -i %(outdir)s//%(amode)s/final.contigs.fa > %(outdir)s/%(amode)s/%(out)s
    ''',
    "contig.tsv": '''
# assemble
metabbq RCAasm.sh %(amode)s /dev/shm BI %(Bar)s BI %(mem)s %(cpu)s %(sam)s/primers.cfg
# RCA
awk -v bc=%(Bar)s '/^>/{{print bc"\\t"$0}}' %(outdir)s/%(amode)s/final.contigs.fa | sed 's/>//;s/ flag=/\\t/;s/ multi=/\\t/;s/ len=/\\t/' > %(outdir)s/%(amode)s/%(out)s
    ''',
    "clip.metadata.tsv": '''
# assemble
metabbq RCAasm.sh %(amode)s /dev/shm BI %(Bar)s BI %(mem)s %(cpu)s %(sam)s/primers.cfg
# RCA
awk -v b=$i -F ' +|\\t|_|=' '/multi=/{{print b"\\t"$2"\\t"$3"\\t"$5"\\t"$7"\\t"$9"\\t"$10"\\t"$11"\\t"$12"\\t"$13}}'  %(outdir)s/%(amode)s/RCAclip.log > %(outdir)s/%(amode)s/%(out)s
    ''',
    "clip.all.fasta": '''
# assemble
metabbq RCAasm.sh %(amode)s /dev/shm BI %(Bar)s BI %(mem)s %(cpu)s %(sam)s/primers.cfg
# RCA
awk -v b=$i 'FNR%2==1{{sub(">",">"b"_",$0);id=$0}}FNR%2==0{{print id"\\n"$0}}' %(outdir)s/%(amode)s/RCAclip.fa > %(outdir)s/%(amode)s/%(out)s
    ''',
     "task": '''
# mkdir %(outdir)s/%(amode)s/
# assemble
metabbq RCAasm.sh %(amode)s /dev/shm BI %(Bar)s BI %(mem)s %(cpu)s %(sam)s/primers.cfg
# RCA contig.fasta
metabbq clusterHelper asmPick -l %(mAsm)s -b %(Bar)s -i %(outdir)s/%(amode)s/final.contigs.fa > %(outdir)s/%(amode)s/BI.contig.fasta
# RCA clip.all.fasta
awk -v b=%(Bar)s 'FNR%%2==1{{sub(">",">"b"_",$0);id=$0}}FNR%%2==0{{print id"\\n"$0}}' %(outdir)s/%(amode)s/RCAclip.fa > %(outdir)s/%(amode)s/BI.clip.all.fasta
# RCA clip.metadata.tsv
awk -v b=%(Bar)s -F ' +|\\t|_|=' '/multi=/{{print b"\\t"$2"\\t"$3"\\t"$5"\\t"$7"\\t"$9"\\t"$10"\\t"$11"\\t"$12"\\t"$13}}'  %(outdir)s/%(amode)s/RCAclip.log > %(outdir)s/%(amode)s/BI.clip.metadata.tsv
# RCA contig.tsv
awk -v bc=%(Bar)s '/^>/{{print bc"\\t"$0}}' %(outdir)s/%(amode)s/final.contigs.fa | sed 's/>//;s/ flag=/\\t/;s/ multi=/\\t/;s/ len=/\\t/' > %(outdir)s/%(amode)s/BI.contig.tsv

    ''',
}


variants = {
    "out": outn, # run_type
    "amode": mode,
    "sam": sdir,
    "Bar": "1",
    "mem": mem,
    "cpu": cpu,
    "mAsm": mAsm,
    "outdir": "./"
}

current_template = cmp_template_map.get("task")
# logger.info("current cmd template: %s" % current_template)
if not current_template:
    logger.error("Error run type.")
    sys.exit(0)


def run_task(shell, path):
    if os.name == 'nt':
        os.system(shell) # windows
        if os.system(shell) != 0:
            raise Exception("Thread exit nonzero")
    else:
        logger.info("Shell: %s, started." % shell)
        os.system("sh %s" % shell)  # linux
        if os.system("sh %s" % shell) != 0:
             raise Exception("Thread exit nonzero")

    return path


def mkdir(path):
    folder = os.path.exists(path)
    if not folder:
        os.makedirs(path)
        logger.info("new folder: %s" % path)
    else:
        logger.info( "There is this folder: %s" % path)


def Idvd(idvd):
    bcs = {}
    bcr = {}
    with open(idvd, 'r') as DD:
        for line in DD:
            line.strip()
            line = line.split("\t")
            bcs[line[0]] = line[1]
            bcr[line[0]] = line[2]
    return bcs, bcr


class TaskGenerate:

    def get_flock(self):
        # Acquire a file lock, the file is operated by other processes and will wait for the lock to be released
        logger.info("Locking...")
        if os.name == 'nt':
            #return directly if windows system
            return
        fcntl.flock(self.output_fh, fcntl.LOCK_EX)
        logger.info("Locked successfully.")

    def release_flock(self):
        if os.name == 'nt':
            #return directly if windows system
            return
        fcntl.flock(self.output_fh, fcntl.LOCK_UN)
        logger.info("Lock released successfully.")

    def __init__(self, thread=8):
        self.thread = thread
        self.result_file = result1
        self.p = None
        self.output_fh = open(self.result_file, 'w')
        self.output_fh2 = open(result2, 'w')
        self.output_fh3 = open(result3, 'w')
        self.output_fh4 = open(result4, 'w')
        self.current_task_number = 0

    def init(self):
        pass


    def run(self, target_func):
        def error_callback(data):
            logger.info("======== run error: %s" % data)
            # 某个样本运行失败，直接退出运行
            self.p.terminate()

        def success(data):
            logger.info("Task complete, current_task_number: %s, return data: %s" % (self.current_task_number, data))
            self.current_task_number -= 1
            if glob.glob(data):
                self.get_flock()
                try:
                    #----BI.contig.tsv
                    outf = os.path.join(data, mode, "BI.contig.tsv")
                    oo = open(outf)
                    res = oo.read()
                    self.output_fh.write(res)
                    oo.close()
                    #---BI.clip.all.fasta
                    outc = os.path.join(data, mode, "BI.clip.all.fasta")
                    oc = open(outc)
                    resc = oc.read()
                    self.output_fh2.write(resc)
                    oc.close()
                    #----BI.clip.metadata.tsv
                    outm = os.path.join(data, mode, "BI.clip.metadata.tsv")
                    om = open(outm)
                    resm = om.read()
                    self.output_fh3.write(resm)
                    om.close()
                    # ----BI.contig.fasta
                    outa = os.path.join(data, mode, "BI.contig.fasta")
                    oa = open(outa)
                    resa = oa.read()
                    self.output_fh4.write(resa)
                    oa.close()
                except Exception as e:
                    logger.error("Collection result error, %s. %s" % (data, str(e)))
                self.release_flock() # Unlock
                shutil.rmtree(data) # Delete the corresponding directory
            else:
                logger.info("no such file")

        self.p = Pool(self.thread)
        max_pending_task = self.thread * 2
        logger.info("specifcBarcode mode. Start.")

        if in1.endswith(".gz"):
            fh1 = gzip.open(in1, "r")
            fh2 = gzip.open(in2, "r")
        else:
            fh1 = open(in1, "r")
            fh2 = open(in2, "r")

        bcs, bcr = Idvd(IDlst)
        for key in bcs:
            # Ensure that the maximum number of delivery tasks does not exceed max_pending_task
            while self.current_task_number > max_pending_task:
                logger.info("waiting--------")
                time.sleep(5)
            # oDir = tdir + "/Assemble_BI/" + str('BI' + "%08d" % int(key))
            oDir = os.path.join(tdir,"Assemble_BI", str('BI' + "%08d" % int(key)))
            mkdir(oDir)
            out1 = open(oDir + "/sort.1.fq", "w")
            out2 = open(oDir + "/sort.2.fq", "w")
            bc = bcs[key]
            rpb = bcr[key]
            ee = bc.replace("_", "")
            rd = 0
            while True:
                line = fh1.readline()
                if not line:
                    break
                iID = line
                barN = iID.split("/")[1]
                barN = barN.replace("_","")
                if ee < barN:
                    logger.info("error: ", ee, "is less than last bb, Make sure the IDs are ordered.")
                    break
                if barN < ee:
                    fh1.readline()
                    fh1.readline()
                    fh1.readline()
                    continue
                rd += 1
                if fmt == "fa":
                    out1.writelines(line)
                    out1.writelines(fh1.readline())
                    fh1.readline()
                    fh1.readline()
                if fmt == "fq":
                    out1.writelines(line)
                    out1.writelines(fh1.readline())
                    out1.writelines(fh1.readline())
                    out1.writelines(fh1.readline())
                    out2.writelines(fh2.readline())
                    out2.writelines(fh2.readline())
                    out2.writelines(fh2.readline())
                    out2.writelines(fh2.readline())
                if rd == int(rpb):
                    batchName = '/batch.assemble.BI' + str("%08d" % int(key)) + '.sh'
                    shellN = oDir + batchName # shell文件
                    fhs = open(shellN, 'w')
                    variants["Bar"] = key
                    variants["outdir"] = oDir
                    bsh = current_template % variants
                    fhs.writelines(bsh)
                    fhs.close()
                    self.p.apply_async(run_task, args=(shellN, oDir), callback=success,
                                       error_callback=error_callback) # 多线程执行多个shell文件
                    self.current_task_number += 1
                    break
            out1.close()
            out2.close()
        fh1.close()
        fh2.close()
        self.p.close()
        self.p.join()

    def generate_next_task(self):
        pass

if __name__ == '__main__':
    task_generator = TaskGenerate(thread)
    task_generator.run(run_task)
