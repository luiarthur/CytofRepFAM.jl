import numpy as np

class Popultion:
    def __init__(self):
        self.all = dict()
        self.size = 0

    def keys(self):
        return self.all.keys()

    def label(self, subpop):
        b = self.__binarize(subpop)

        if b not in self.keys():
            self.size += 1
            self.all[b] = self.size 

        return self.all[b]

    def __binarize(self, subpop):
        subpop_char_list = [str(int(s)) for s in subpop]
        return ''.join(subpop_char_list)

# Example usage:
# pop = Popultion()
# pop.label([1,0,1])  # 1
# pop.label([1,0,1])  # 1
# pop.label([1,1,1])  # 2
# pop.label([1,1,1,1])  # 3
