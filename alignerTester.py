#! /usr/bin/env python3.2

from homPolDist import writeBothToFile, fastaParser, stripLowerCase
from subprocess import call, check_output
from Bio.Seq import translate
from math import log10
from numpy import median,std
from Bio.SubsMat.MatrixInfo import blosum62 as blosum
import binascii
import sys

blosum62 = {}
for (a,b),val in blosum.items():
    oldT = (a,b)
    newT = (b,a)
    blosum62[tuple(newT)] = val
    blosum62[tuple(oldT)] = val
randGeneFN = "randGene.fasta"
transGeneFN = "transGene.fasta"
# Trans gene no homopolymer errors
tnhpReadsFN = "tnhpReads.fasta"
randSffFN = "randFlow.sff"
transSffFN = "transFlow.sff"
randReadsFN = "randReads.fasta"
transReadsFN = "transReads.fasta"
randCodonAlnInFile = "randCodonAlnInFile.fasta"
transCodonAlnInFile = "transCodonAlnInFile.fasta"
tnhpCodonAlnInFile = "tnhpCodonAlnInFile.fasta"
randCodonAlnOutFile = "randCodonAlnOutFile.fasta"
transCodonAlnOutFile = "transCodonAlnOutFile.fasta"
tnhpCodonAlnOutFile = "tnhpCodonAlnOutFile.fasta"
randNucBF = "randNuc.bf"
transNucBF = "transNuc.bf"
tnhpNucBF = "tnhpNuc.bf"
randCodonBF = "randCodon.bf"
transCodonBF = "transCodon.bf"
tnhpCodonBF = "tnhpCodon.bf"

aaDict = {'A': 0,  'R': 1,  'N': 2,  'D': 3,  'C': 4,  'Q': 4,  'E': 5,  
    'G': 6,  'H': 7,  'I': 8,  'L': 9,  'K': 10,  'M': 11,  'F': 12,  
    'P': 13,  'S': 14,  'T': 15, 'W': 16,  'Y': 17,  'V': 18,  'B': 19,  
    'Z': 20,  'X': 21, '*': 22}

def writeThreeDatasets(distSourceFile, geneSize, readNumber):
# build the random and the untranslated genes
    writeBothToFile(distSourceFile, geneSize, 
            randGeneFN, transGeneFN)
# generate homopolymer errored reads from the random gene
    call("clonesim -c %d %s | kitsim | flowsim -o %s" 
            % (int(readNumber), randGeneFN, randSffFN), 
            shell=True)
    call("clonesim -c %d %s | kitsim | flowsim -o %s" 
            % (int(readNumber), transGeneFN, transSffFN),
            shell=True)
    call("clonesim -c %d %s | kitsim > %s" 
            % (int(readNumber), transGeneFN, tnhpReadsFN),
            shell=True)
# translate these SFF files to fasta files.
    call("flower %s -f=%s "
            % (randSffFN, randReadsFN), shell=True)
    call("flower %s -f=%s "
            % (transSffFN, transReadsFN), shell=True)

# write out the alignment scores for the three methods
def scoreMethods(distSourceFile, geneSize, readNumber):
    writeThreeDatasets(distSourceFile, geneSize, readNumber)
# get read files into batch files
    writeHyphyBatchFiles()

    """
    call("~/Software/hyphy/hyphy/HYPHYMP %s" % randCodonBF, shell=True)
    print("")
    print("")
    print("")
    print("")
    """
    randResults = check_output("~/Software/hyphy/hyphy/HYPHYMP %s" % randCodonBF, shell=True)
    transResults = check_output("~/Software/hyphy/hyphy/HYPHYMP %s" % transCodonBF, shell=True)
    tnhpResults = check_output("~/Software/hyphy/hyphy/HYPHYMP %s" % tnhpCodonBF, shell=True)
    for line in transResults.split():
        print(line)
    print("")
    """
    print(transResults)
    print("")
    print(tnhpResults)
    print("")
    """

    randCodonAlignments = parseCodonResults(randResults)
    transCodonAlignments = parseCodonResults(transResults)
    tnhpCodonAlignments = parseCodonResults(tnhpResults)

    randNucAlignments = parseNucResults(randResults)
    transNucAlignments = parseNucResults(transResults)
    tnhpNucAlignments = parseNucResults(tnhpResults)

    """
    print(randNucAlignments)
    print("")
    print(transNucAlignments)
    print("")
    print(tnhpNucAlignments)
    print("")
    """

    randCodonScoreList = scoreAlignmentList(randCodonAlignments)
    transCodonScoreList = scoreAlignmentList(transCodonAlignments)
    tnhpCodonScoreList = scoreAlignmentList(tnhpCodonAlignments)

    randNucScoreList = scoreAlignmentList(randNucAlignments)
    transNucScoreList = scoreAlignmentList(transNucAlignments)
    tnhpNucScoreList = scoreAlignmentList(tnhpNucAlignments)

    """
    for rqTuple in transCodonAlignments:
        gaplessTuple = gaplessStrings(rqTuple)
        print(gaplessTuple[0])
        print("")
        print(gaplessTuple[1])
        print("")
        print("")
    print("")
    print("nuc: ")
    print("")
    for rqTuple in transNucAlignments:
        gaplessTuple = gaplessStrings(rqTuple)
        print(gaplessTuple[0])
        print("")
        print(gaplessTuple[1])
        print("")
        print("")
    print("codon scores: ")
    print(transCodonScoreList)
    print("")
    print("nuc scores: ")
    print(transNucScoreList)
    print("")
    """

    fullCodonList = []
    fullCodonList.extend(randCodonScoreList)
    fullCodonList.extend(transCodonScoreList)
    fullCodonList.extend(tnhpCodonScoreList)

    fullNucList = []
    fullNucList.extend(randNucScoreList)
    fullNucList.extend(transNucScoreList)
    fullNucList.extend(tnhpNucScoreList)

    medianCodonScore = median(fullCodonList)
    medianNucScore = median(fullNucList)

    for i in range(len(transCodonAlignments)):
        print("Nuc: ", transNucAlignments[i])
        print("Sco: ", transNucScoreList[i])
        print("Cod: ", transCodonAlignments[i])
        print("Sco: ", transCodonScoreList[i])



    print("Cod: ", fullCodonList)
    print(medianCodonScore)
    print(std(fullCodonList))

    print("Nuc:", fullNucList)
    print(medianNucScore)
    print(std(fullNucList))

    print("Codon: ")
    print("rand: ")
    print(median(randCodonScoreList), " ", std(randCodonScoreList))
    print("trans: ")
    print(median(transCodonScoreList), " ", std(transCodonScoreList))
    print("tnhp: ")
    print(median(tnhpCodonScoreList), " ", std(tnhpCodonScoreList))
    print("")

    print("Nuc: ")
    print("rand: ")
    print(median(randNucScoreList), " ", std(randNucScoreList))
    print("trans: ")
    print(median(transNucScoreList), " ", std(transNucScoreList))
    print("tnhp: ")
    print(median(tnhpNucScoreList), " ", std(tnhpNucScoreList))

def performCodonAlignment(distSourceFile, geneSize, readNumber):
    writeThreeDatasets(distSourceFile, geneSize, readNumber)
# write datasets for codon alignment
    call("cat %s %s > %s" % (randGeneFN, randReadsFN, randCodonAlnInFile), shell=True)
    call("cat %s %s > %s" % (transGeneFN, transReadsFN, transCodonAlnInFile), shell=True)
    call("cat %s %s > %s" % (transGeneFN, tnhpReadsFN, tnhpCodonAlnInFile), shell=True)
# use codonAligner to align
    call("codonAligner %s %s" % (randCodonAlnInFile, randCodonAlnOutFile), shell=True)
    call("codonAligner %s %s" % (transCodonAlnInFile, transCodonAlnOutFile), shell=True)
    call("codonAligner %s %s" % (tnhpCodonAlnInFile, tnhpCodonAlnOutFile), shell=True)


def writeHyphyBatchFiles():
    randReads = fastaParser(randReadsFN)
    transReads = fastaParser(transReadsFN)
    tnhpReads = fastaParser(tnhpReadsFN)

    randRef = fastaParser(randGeneFN)
    transRef = randRef
    tnhpRef = randRef
#transRef = fastaParser(transGeneFN)
#tnhpRef = fastaParser(transGeneFN)

    randCodonB = generateCodonAndNucB(randRef, randReads)
    transCodonB = generateCodonAndNucB(transRef, transReads)
    tnhpCodonB = generateCodonAndNucB(tnhpRef, tnhpReads)

    """
    randNucB = generateCodonAndNucB(randRef, randReads)
    transNucB = generateCodonAndNucB(transRef, transReads)
    tnhpNucB = generateCodonAndNucB(tnhpRef, tnhpReads)
    """

    writeBF(randCodonB, randCodonBF)
    writeBF(transCodonB, transCodonBF)
    writeBF(tnhpCodonB, tnhpCodonBF)

    """
    writeBF(randNucB, randNucBF)
    writeBF(transNucB, transNucBF)
    writeBF(tnhpNucB, tnhpNucBF)
    """

def _quicksize(value):
    return int(log10(max(1, value)))

def generateCodonAndNucB(ref, reads):  
    optsI = 0
    codonB = ""
    codonB += "opts = {};\n"
    codonB += "opts[ "
    codonB += '"%s%d"' % ('0' * (_quicksize(len(reads) + 3) - _quicksize(optsI)), optsI)
    codonB += " ] = \"No\";\n"
    optsI += 1
    codonB += "opts[ "
    codonB += '"%s%d"' % ('0' * (_quicksize(len(reads) + 3) - _quicksize(optsI)), optsI)
    codonB += " ] = \""
    optsI += 1
    codonB += ''.join(ref[0])
    codonB += "\";\n"
    codonB += "opts[ "
    codonB += '"%s%d"' % ('0' * (_quicksize(len(reads) + 3) - _quicksize(optsI)), optsI)
    codonB += "] = \"0.0\";\n"
    optsI += 1
    codonB += "opts[ "
    codonB += '"%s%d"' % ('0' * (_quicksize(len(reads) + 3) - _quicksize(optsI)), optsI)
    codonB += " ] = \""
    optsI += 1
    codonB += str(len(reads))
    codonB += "\";\n"
#index = 4
    for read in reads:
        codonB += "opts[ "
#        codonB += str(index)
        codonB += '"%s%d"' % ('0' * (_quicksize(len(reads) + 3) - _quicksize(optsI)), optsI)
        optsI += 1
#        index += 1
        codonB += " ] = \""
        codonB += ''.join(stripLowerCase(read))
        codonB += "\";\n"

    codonB += "alnopts = {};\n"
    codonB += "alnopts[ \"SEQ_ALIGN_CHARACTER_MAP\" ] = \"ACGT\";\n"
    codonB += "alnopts[ \"SEQ_ALIGN_SCORE_MATRIX\" ] = {{ 5, -4, -4, -4}, {-4, 5, -4, -4}, {-4, -4, 5, -4}, {-4, -4, -4, 5}};\n"
    codonB += "alnopts[ \"SEQ_ALIGN_GAP_OPEN\" ] = 40;\n"
    codonB += "alnopts[ \"SEQ_ALIGN_GAP_OPEN2\" ] = 40;\n"
    codonB += "alnopts[ \"SEQ_ALIGN_GAP_EXTEND\" ] = 1;\n"
    codonB += "alnopts[ \"SEQ_ALIGN_GAP_EXTEND2\" ] = 1;\n"
    codonB += "alnopts[ \"SEQ_ALIGN_AFFINE\" ] = 1;\n"
    codonB += "alnopts[ \"SEQ_ALIGN_NO_TP\" ] = 1;\n"

    refPullingString = ""
    refPullingString += "inseqs = {{ opts[ "
    refPullingString += '"%s%d"' % ('0' * (_quicksize(len(reads) + 3) - _quicksize(1)), 1)
    refPullingString += " ], "
    for i in range(len(reads)):
        codonB += refPullingString
        codonB += "opts[ "
        codonB += '"%s%d"' % ('0' * (_quicksize(len(reads) + 3) - _quicksize(4+i)), i+4)
        codonB += " ]"
        codonB += " }};\n" 
        codonB += "AlignSequences( alnres, inseqs, alnopts );\n"
        codonB += "fprintf( stdout, \"\" + alnres + \"\\n\" );\n"

    codonB += "ExecuteAFile( \"codonaligner.bf\", opts );\n"
    codonB += "fprintf( stdout, _cdnaln_outstr + \"\\n\" );"
    return codonB

def writeBF(batchString, bfName):
    bf = open(bfName, 'w')
    bf.write(batchString)
    bf.close()

def parseCodonResults(resultsString):
    alignmentList = []
    resultsList = resultsString.split()
    codonStart = 0
    codonEnd = 0
    for i,line in enumerate(resultsList):
        for c in str(line, encoding='ascii'):
            if c == "}":
                codonStart = i + 1
    for i,line in enumerate(resultsList):
        for c in str(line, encoding='ascii'):
            if c == "[":
                codonEnd = i
    j = codonStart
    while j < codonEnd:
        ref = str(resultsList[j], encoding='ascii')
        query = str(resultsList[j+1], encoding='ascii')
        alignmentList.append([ref, query])
        j += 2
    return(alignmentList)

def parseNucResults(resultsString):
    alignmentList = []
    resultsList = resultsString.split()
    nucStart = -1
    nucEnd = -1
    for i,line in enumerate(resultsList):
        for c in str(line, encoding='ascii'):
            if c == "{":
                nucStart = i
        if nucStart != -1:
            break
    for i,line in enumerate(resultsList):
        for c in str(line, encoding='ascii'):
            if c == "}":
                nucEnd = i + 1
    j = nucStart
    while j < (nucEnd -2):
        ref = ''.join(c.upper() for c in str(resultsList[j], encoding='ascii') if c.isalpha() or c == "-")
        query = ''.join(c.upper() for c in str(resultsList[j + 1], encoding='ascii') if c.isalpha() or c == "-")
        alignmentList.append([ref, query])
        j += 3
    return(alignmentList)

def scoreAlignmentList(alignmentList):
    scoreList = []
    queriesWithStops = []
    minScore = 0
    for al in alignmentList:
        stopCodons = 0
        queryScore = 0
        alN = gaplessStrings(al)
        ref = alN[0]
        query = alN[1]
        refAA = translate(ref)
#print(refAA)
        queryAA = translate(query)
#print(queryAA)
        for c in refAA:
            if c == '*':
                stopCodons += 1
        for c in queryAA:
            if c == '*':
                stopCodons += 1
        if stopCodons == 0:
            queryScore = sum(score_pairwise(refAA, queryAA, blosum62, -5, -1))
#print(queryScore)
            scoreList.append(queryScore)
            if minScore > queryScore:
                minScore = queryScore
        else:
            queriesWithStops.append(stopCodons)
#print("stopCodon")
#print("")
    for stopCodonCount in queriesWithStops:
        scoreList.append(stopCodonCount * minScore)
    return scoreList

def gaplessStrings(stringT):
    j = 0
    tbrR = ""
    tbrQ = ""
    for cI in range(0, min(len(stringT[0]),len(stringT[1]))-2, 3):
        if (stringT[0][cI] == '-' and stringT[0][cI+1] == '-' and stringT[0][cI+2] == '-'):
            pass
        elif (stringT[1][cI] == '-' and stringT[1][cI+1] == '-' and stringT[1][cI+2] == '-'):
            pass
        else:
            tbrR += stringT[0][cI]
            tbrR += stringT[0][cI + 1]
            tbrR += stringT[0][cI + 2]
            tbrQ += stringT[1][cI]
            tbrQ += stringT[1][cI + 1]
            tbrQ += stringT[1][cI + 2]
    tbrR = tbrR.replace('-', 'N')
    tbrQ = tbrQ.replace('-', 'N')
    return (tbrR, tbrQ)

def score_pairwise(seq1, seq2, matrix, gap_s, gap_e, gap = True):
    for A,B in zip(seq1, seq2):
        diag = ('-'==A) or ('-'==B)
        yield (gap_e if gap else gap_s) if diag else matrix[(A,B)]
        gap = diag

if __name__ == "__main__":
    scoreMethods(*sys.argv[1:])
