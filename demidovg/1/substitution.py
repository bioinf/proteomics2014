from collections import defaultdict

class Substitution_matrix():
    matrix = defaultdict(dict)

    def __init__(self, filename):
        with open(filename) as f:
            top = True
            top_string = ""
            for line in f:
                if not line.startswith("#") and top == True:
                    top_string = line.split()
                    top = False
                elif not line.startswith("#") and top == False:
                    weights = line.split()
                    for i in xrange(1,len(weights)):
                        self.matrix[top_string[i - 1]][weights[0]] = int(weights[i])

    def wt(self, letter_one, letter_two):
        try:
            return self.matrix[letter_one][letter_two]
        except KeyError:
            return self.matrix["*"]["*"]


