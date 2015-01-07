

class Blosum(object):
    def __init__(self, filename):
        self.matrix = {}

        with open(filename, 'r') as input:
            line = input.readline()
            while line[0] =='#':
                line = input.readline()

            tokens = line.split()
            line = input.readline()

            while line:
                blocks = line.split()
                self.matrix[blocks[0]] = {}
                for i in range(len(tokens)):
                    self.matrix[blocks[0]][tokens[i]] = int(blocks[i + 1])
                line = input.readline()

    def keys(self):
        return self.matrix.keys()

    def get_elem(self, a, b):
        return self.matrix[a][b]