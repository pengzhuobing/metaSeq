#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 11 16:59:32 2017

@author: songz
"""
from __future__ import print_function
from __future__ import division

# an iterator oobject for reading a single sequence file
# It will NOT check the format of the file, either can it deal with multiple line FASTA file.
class sequence(object):
    def __init__(self, filePath, fastx='a', gz=False, trunk_size=1):
        self.fastx = fastx
        self.gzip = gz
        self.size = trunk_size
        if self.gzip:
            import gzip
            self.file = gzip.open(filePath, 'rt')
        else:
            self.file = open(filePath, 'r')
        if fastx == 'a':
            self.n = 2
        elif fastx == 'q':
            self.n = 4
        else:
            print('Please specify the right format, "a" for FASTA and "q" for FASTQ.')
    
    def __iter__(self):
        return self
    
    def __next__(self):
        record = []
        for i in range(self.n):
            line = self.file.readline().strip('\n')
            if line:
                record.append(line)
            else:
                raise StopIteration
        record[0] = record[0][1:]
        return record

# Same iterator but read in multiple record into memory at once
class sequence_trunk(object):
    def __init__(self, filePath, fastx='a', gz=False, trunk_size=2):
        self.fastx = fastx
        self.gzip = gz
        self.trunk_size = trunk_size
        if self.gzip:
            import gzip
            self.file = gzip.open(filePath, 'rt')
        else:
            self.file = open(filePath, 'r')
        if fastx == 'a':
            self.n = 2
        elif fastx == 'q':
            self.n = 4
        else:
            print('Please specify the right format, "a" for FASTA and "q" for FASTQ.')
            self.n = 1
    
    def __iter__(self):
        return self
    
    def __next__(self):
        record_trunk = []
        for record in range(self.trunk_size):
            record = []
            for i in range(self.n):
                line = self.file.readline().strip('\n')
                if line:
                    record.append(line)
                else:
                    if len(record_trunk) > 0:
                        return record_trunk
                    else:
                        raise StopIteration
            record[0] = record[0][1:]
            record_trunk.append(record)
        return record_trunk


# Iterator for two files
# It only work for files with ABSOLUTELY corresponding record.
class sequence_twin(object):
    def __init__(self, file_r1, file_r2, fastx='a', gz=False):
        self.fastx = fastx
        self.gzip = gz
        if self.gzip:
            import gzip
            self.r1 = gzip.open(file_r1, 'rt')
            self.r2 = gzip.open(file_r2, 'rt')
        else:
            self.r1 = open(file_r1, 'r')
            self.r2 = open(file_r2, 'r')
        if fastx == 'a': self.n = 2
        elif fastx == 'q': self.n = 4
        else:
            print('Please specify the right format, "a" for FASTA and "q" for FASTQ.')
            self.n = 1
    
    def __iter__(self):
        return self
    
    def __next__(self):
        record = [[],[]]
        for i in range(self.n):
            line_r1 = self.r1.readline().strip('\n')
            line_r2 = self.r2.readline().strip('\n')
            if line_r1:
                record[0].append(line_r1)
                record[1].append(line_r2)
            else:
                raise StopIteration
        record[0][0] = record[0][0][1:]
        record[1][0] = record[1][0][1:]
        return record[0], record[1]


class sequence_twin_trunk(object):
    def __init__(self, file_r1, file_r2, fastx='a', gz=False, trunk_size=2):
        self.fastx = fastx
        self.gzip = gz
        if self.gzip:
            import gzip
            self.r1 = gzip.open(file_r1, 'rt')
            self.r2 = gzip.open(file_r2, 'rt')
        else:
            self.r1 = open(file_r1, 'r')
            self.r2 = open(file_r2, 'r')
        if fastx == 'a': self.n = 2
        elif fastx == 'q': self.n = 4
        else:
            print('Please specify the right format, "a" for FASTA and "q" for FASTQ.')
            self.n = 1
        self.trunk_size = trunk_size

    def __iter__(self):
        return self
    
    def __next__(self):
        r1_trunk = []
        r2_trunk = []
        for record in range(self.trunk_size):
            r1 = []
            r2 = []
            for i in range(self.n):
                line_r1 = self.r1.readline().strip('\n')
                line_r2 = self.r2.readline().strip('\n')
                if line_r1:
                    r1.append(line_r1)
                    r2.append(line_r2)
                else:
                    if len(r1_trunk) > 0:
                        return r1_trunk, r2_trunk
                    else:
                        raise StopIteration
            r1[0] = r1[0][1:]
            r2[0] = r2[0][1:]
            r1_trunk.append(r1)
            r2_trunk.append(r2)
        return r1_trunk, r2_trunk


# This function is still under test, it reads in file bytes, should be a bit faster than read in by line
# Need to fins a way to remove header
class sequence_fastq_bytes(object):
    def __init__(self, filePath, size, ):
        self.file = open(filePath, 'r')
        self.size = int(size)
        self.tail = ''
    def __iter__(self):
        return self
    
    def __next__(self):
        content = self.file.read(self.size)
        #print(self.tail)
        if content:
            content = self.tail + content
            self.tail = ''
            content = content.split('\n')
            self.tail = content[-1]
            tail_n = (len(content) - 1) % 4 # Set the line step to 4, set to 2 for FASTA
            if tail_n > 0:
                self.tail = '\n'.join(content[:-1][-1*tail_n::]) + '\n' + self.tail
            else:
                self.tail = self.tail
            #print(self.tail)
            content = content[:(len(content) - tail_n - 1)]
            content = [content[x:x+4] for x in range(0, len(content), 4)]
            return content
        else:
            if self.tail:
                content = self.tail.split('\n')
                content = [content[x:x+4] for x in range(0, len(content), 4)]
                return content
            else:
                raise StopIteration


# Write the content to a fastq file
def write_seqs(seq_content, filePath, fastx='a', gz=False, mode='w'):
    count = 0
    if fastx == 'a':
        n = 2
        header = '>'
    elif fastx == 'q':
        n = 4
        header = '@'
    else:
        n = 1
        header = ''
    
    if not gz:
        f = open(filePath, mode)
    else:
        import gzip
        f = gzip.open(filePath, mode)
    for record in seq_content:
        label = header + record[0]                
        for line in [label] + record[1:]:
            f.write('%s\n' % line)
            count += 1
    f.close()
    return count


# Parse a sorted stLFR FASTA data by bead (barcode)
# Return a tuple containing all the short sequences in the bead
# Currently only support FASTA since all QC should be at the upstream
class stlfr_bead(object):
    def __init__(self, filePath):
        self.stop = False
        with open(filePath, 'r') as f:
            line1 = f.readline()
            self.barcode = line1.split('-')[0]
        self.file = open(filePath, 'r')
        self.current_bead = []
    
    def __iter__(self):
        return self
    
    def __next__(self):
        close = False
        while close == False:
            current_bead = self.current_bead
            record = []
            for i in range(2): # Read in the next sequence
                line = self.file.readline().strip('\n')
                if line:
                    record.append(line)
                else:
                    if self.stop:
                        raise StopIteration
                    else:
                        self.stop = True
                        current_bead = (self.barcode, tuple((tuple(i) for i in current_bead)))
                        return current_bead
            barcode = record[0].split('-')[0]
            if barcode == self.barcode:
                current_bead.append(record)
            else:
                current_bead = (self.barcode, tuple((tuple(i) for i in current_bead)))
                self.barcode = barcode
                self.current_bead = [record]
                close = True
                return current_bead


#%%
# Return reverse compliment of a sequence
# This part is got from Stakoverflow
#(https://stackoverflow.com/questions/19570800/reverse-complement-dna) by corinna
# Works for Python 3
def revcomp(seq):
    return seq.translate(str.maketrans('ACGTacgtRYMKrymkVBHDvbhdN', 'TGCAtgcaYRKMyrkmBVDHbvdhN'))[::-1]