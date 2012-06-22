#! /usr/bin/env python

import math

def checkArray(array, index):
    while index >= len(array):
        array.append(0)

def getCounts(listOfCharacters):
    listOfCounts = []
    currentCharacter = 'A'
    currentCount = 0
    for c in listOfCharacters:
        if c == currentCharacter:
            currentCount += 1
        else:
            checkArray(listOfCounts, currentCount)
            listOfCounts[currentCount] += 1
            currentCount = 1
            currentCharacter = c
    listOfCounts[0] = 0
    return listOfCounts

def average(array):
    arraySum = 0
    numberOfSubstrings = 0
    for i,value in enumerate(array):
        arraySum += i*value
        numberOfSubstrings += value
    return arraySum/numberOfSubstrings

def averageOfSquares(array):
    arraySum = 0
    numberOfSubstrings = 0
    for i,value in enumerate(array):
        arraySum += (i*value)**2
        numberOfSubstrings += value
    return arraySum/numberOfSubstrings
        
def variance(array):
    arraySum = 0
    numberOfSubstrings = 0
    mu = average(array)
    for i,value in enumerate(array):
        for j in range(value):
            arraySum += (i-mu)**2
        numberOfSubstrings += value
    return arraySum/numberOfSubstrings

def logAverage(array):
    arraySum = 0
    numberOfSubstrings = 0
    for i,value in enumerate(array):
        for j in range(value):
            arraySum += math.log(i)
        numberOfSubstrings += value
    return arraySum/numberOfSubstrings

def logVariance(array):
    arraySum = 0
    numberOfSubstrings = 0
    lnAv = logAverage(array)
    for i,value in enumerate(array):
        for j in range(value):
            arraySum += (math.log(i) - lnAv)**2
        numberOfSubstrings += value
    return arraySum/numberOfSubstrings
        
def estimateExponentialPDF(array):
    PDF = [0]*20
    arrayAverage = average(array)
    expLambda = 1/arrayAverage
    for i in range(len(PDF)):
        PDF[i] = expLambda * math.exp(-1 * expLambda * i)
    return PDF 
        
def estimateExponentialCDF(array):
    CDFvalues = [0]*20
    arrayAverage = average(array)
    expLambda = 1/arrayAverage
    for i in range(len(CDFvalues)):
        CDFvalues[i] = 1 - math.exp(-1 * expLambda * i)
    thresholdValues = [0]*20
    for i in range(1, len(thresholdValues)):
        thresholdValues[i] = CDFvalues[i] - CDFvalues[i-1]
    return thresholdValues

def estimateLogNormalThresholds(countArray):
    PDFvalues = [0]*20
    mu = logAverage(countArray)
    sigma = math.sqrt(logVariance(countArray))
    for i in range(1,len(PDFvalues)):
        expon = -1 * ((math.log(i) - mu)**2)
        expon /= 2 * (sigma**2)
        coef = 1/(i * sigma * math.sqrt(2*math.pi))
        PDFvalues[i] = coef * math.exp(expon)
    return PDFvalues

def estimateNormPDF(countArray):
    PDF = [0] * 20
    mu = average(countArray)
    sigma = math.sqrt(variance(countArray))
    for i in range(1, len(PDF)):
        coef = 1 / (math.sqrt(2 * math.pi * sigma**2))
        expn = ((i - mu)**2 * -1)/(2 * sigma**2)
        PDF[i] = coef * math.exp(expn)
    return PDF

def empiricalPDF(countArray):
    PDF = [0] * 20
    numSubst = 0
    for i, value in enumerate(countArray):
        PDF[i] = value
        numSubst += value
    for i in range(len(PDF)):
        PDF[i] /= numSubst
    return PDF

def nChoosek(n, k):
    return math.factorial(n)/(math.factorial(k)*math.factorial(n-k))

def betaBinomialPDF(countArray):
    PDF = [0] * 20
    numSubst = sum(countArray)
    n = numSubst
    for i, value in enumerate(countArray):
        PDF[i] = value
    m1 = average(countArray)
    m2 = averageOfSquares(countArray)

    alphaHat = numSubst * m1 - m2
    alphaHat /= numSubst * ((m2 / m1) - m1 - 1) + m1
    betaHat = (n - m1) * (n - m2/m1)
    betaHat /= n * (m2/m1 - m1 - 1) + m1

    for i, value in enumerate(PDF):
        PDF[i] = nChoosek(len(PDF), i)

    return PDF

def roundTo(x, base=3):
    return int(base * round(float(x)/base))
def roundDownTo(x, base=3):
    return int((x/base)*base)






