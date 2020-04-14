#!/usr/bin/env python

import os, re, collections

localDir = os.path.dirname(__file__)
absDir = os.path.join(os.getcwd(), localDir)

class Suggestion(object):
    def __init__(self, taxonomy):
        file_name = "spell_file_%s.txt" % taxonomy
        #NWORDS = train(words(open('spell_file.txt').read()))
        self.NWORDS = self.train(openos.path.join(absDir, file_name)).read().split('\n'))
        self.alphabet = 'abcdefghijklmnopqrstuvwxyz'

    def words(self, text): return re.findall('[a-z]+', text.lower()) 

    def train(self, features):
        model = collections.defaultdict(lambda: 1)
        for f in features:
            model[f] += 1
        return model

    def edits1(self, word):
        n = len(word)
        return set([word[0:i]+word[i+1:] for i in range(n)] +                     # deletion
                   [word[0:i]+word[i+1]+word[i]+word[i+2:] for i in range(n-1)] + # transposition
                   [word[0:i]+c+word[i+1:] for i in range(n) for c in self.alphabet] + # alteration
                   [word[0:i]+c+word[i:] for i in range(n+1) for c in self.alphabet])  # insertion

    def known_edits2(self, word):
        return set(e2 for e1 in self.edits1(word) for e2 in self.edits1(e1) if e2 in self.NWORDS)

    def known(self, words): return set(w for w in words if w in self.NWORDS)

    def correct(self, word):
        candidates = self.known([word]) or self.known(self.edits1(word)) or self.known_edits2(word) or [word]
        return candidates
        #return max(candidates, key=lambda w: self.NWORDS[w]) # Only the best result

if __name__ == '__main__':
    import sys
    sugg = Suggestion('ncbi')
    print(sugg.correct(sys.argv[1]))
