#! /usr/bin/env python

import sys
import hompol
import math
from Bio.Seq import translate
from BioExt import untranslate, randgene
from BioExt.stats import bbinom
from scipy.stats import rv_discrete

# A script essentially. Takes the input file containing a sequence we're looking to model
# and an output gene length
def writeBothToFile(seqFileName, newSeqLength, randFileName, transFileName):
    randOutFile = open(randFileName, "w")
    randOutFile.write("> Randomly generated gene of length %d \n" % int(newSeqLength))
    transOutFile = open(transFileName, "w")
    transOutFile.write("> AA translation and untranslation of randomly generated gene of length %d\n" % int(newSeqLength))

    sequence = generateRandomGene(seqFileName, newSeqLength)
    randOutFile.write(sequence)
    seqPerm = untranslate(translate(sequence))
    transOutFile.write(seqPerm)

    randOutFile.close()
    transOutFile.close()

# Generate the random gene using a combination of my homopol counting module 
# and Lance's everything else. 
def generateRandomGene(seqFileName, newSeqLength):
    f = open(seqFileName, 'r')
    characterList = []
    for line in f:
        if line[0] != '>':
            characterList += list(line.rstrip('\n'))
    homopolCounts = hompol.getCounts(characterList)[1:]
    
    n = len(homopolCounts) - 1
    b = bbinom(*bbinom.fit(homopolCounts))
    randseq = randgene(hompol.roundTo(newSeqLength, 3), b.ppf)
    return randseq

def fastaParser(readFileName):
    reads = []
    readFile = open(readFileName, 'r')
    tempRead = []
    for line in readFile:
        if line[0] == ">":
            if tempRead != []:
                reads.append(tempRead)
                tempRead = []
        else: 
            tempRead += list(line.rstrip('\n'))
    if tempRead != []:
        reads.append(tempRead)
    return reads

def stripLowerCase(read):
    stripped = []
    for char in read:
        if char.isupper():
            stripped.append(char)
    return stripped


if __name__ == "__main__":
    writeBothToFile(*sys.argv[1:])



