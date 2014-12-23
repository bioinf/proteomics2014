__author__ = 'german'

class align_table():
    seq1 = None
    seq2 = None
    table = None
    predacessors = None
    def __init__(self, seq1, seq2):
        self.table = [[0 for j in range(len(seq2) + 1)] for i in range(len(seq1) + 1)]
        self.predacessors = [[[0,0,""] for j in range(len(seq2) + 1)] for i in range(len(seq1) + 1)]
        self.seq1 = seq1
        self.seq2 = seq2

    def __init__(self, seq1, seq2, par):
        if par == 0:
            self.table = [[0 for j in range((seq2) + 1)] for i in range((seq1) + 1)]
            self.predacessors = [[[0,0,""] for j in range((seq2) + 1)] for i in range((seq1) + 1)]
            self.predacessors[0][0] = [-1, -1]
        if par == 1:
            self.table = [[0 for j in range(len(seq2) + 1)] for i in range(len(seq1) + 1)]
            self.predacessors = [[[0,0,""] for j in range(len(seq2) + 1)] for i in range(len(seq1) + 1)]
            self.seq1 = seq1
            self.seq2 = seq2


    def set_pred_top(self, i, j):
        self.predacessors[i][j] = (i - 1, j, "top")

    def set_pred_left(self, i, j):
        self.predacessors[i][j] = (i, j - 1, "left")

    def set_pred_diag(self, i, j):
        self.predacessors[i][j] = (i - 1, j - 1, "diag")

