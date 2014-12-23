__author__ = 'german'

class Sequences():
    seqs = {}
    def __init__(self, filename):
        with open(filename) as f:
            name_of_seq = None
            seq = ""
            for line in f:
                if line.startswith(">"):
                    if name_of_seq:
                        self.seqs[name_of_seq] = seq
                    name_of_seq = line.strip()
                    seq = ""
                else:
                    seq += line.strip().upper()
            self.seqs[name_of_seq] = seq
