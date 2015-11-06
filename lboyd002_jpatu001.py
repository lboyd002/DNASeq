#-------------------------------------------------------------------------------
# Author:      amir
# Created:     10/25/2015
#
# Instructions:
#
# 1) Make sure to rename the file (studentNetId.py) to your netId. (Do not include your first name, last name ... or any extra character)
# 2) To run the program type the following statement in the command line:  
#       -) python studentNetId.py DNASeq1FilePath DNASeq2FilePath OutputFilePath                                                                   
#    where  DNASeq1FilePath is the path to the file that contains First DNA sequence (e.g. DNASeq1.txt)
#           DNASeq2FilePath is the path to the file that contains Second DNA sequence (e.g. DNASeq2.txt)
#           OutputFilePath is the path that the output is goint to be saved (e.g. out.txt)
# 3) Make sure to replace FirstName, LastName, SectionNumber, NetId in studentInfo with your information
# 4) You may add as many functions as you want to this program
# 5) The core function in your program is DNASeqAlignment function, where it takes three arguments (DNASeq1,DNASeq2,outputPath) and 
#    computes the similarityScore, sequenceAlignment1 and sequenceAlignment2. At the end, the function writes the result to the output file (Do not make any changes to the output section).
# 6) sequenceAlignment1 and sequenceAlignment2 are strings and they are composed of following characters: 'A', 'T', 'G', 'C' and '-', Do not include extra space or any other character in the strings.
# 7) Make sure your program works with at least one of the following python versions: (2.7.9, 2.7.8, 2.6.6)
# 8) Once you have tested your program with one of the versions listed above, assign that version number to pythonVersion in studentInfo function
# 9) Make sure to write enough comments in order to make your code easy to understand. 
# 10) Describe your algorithm in ALGORITHM section below (you may add as many lines as you want).
# 11) To understand the flow of the program consider the following example:
#      0) Let say we have DNASeq1.txt file which contains AACCTGACATCTT and DNASeq2.txt file contains CCAGCGTCAACTT
#      1) If we execute the following command in the command line: -) python studentNetId.py DNASeq1.txt DNASeq2.txt out.txt
#      2) input arguments are parsed       
#      3) studentInfo() function will be executed and the output will be saved in out.txt file
#      4) DNASeqAlignment() function will be called
#      5) At the entry of the DNASeqAlignment function, DNASeq1='AACCTGACATCTT' and DNASeq2='CCAGCGTCAACTT'
#      6) You should compute the sequence alignment of DNASeq1 and DNASeq2. Let say the result is as follows:
#       A A C C T G A C - - - - A T C T T
#       | | | | | | | | | | | | | | | | |
#       - - C C A G - C G T C A A - C T T      
#      7) At the end of the DNASeqAlignment function sequenceAlignment1='AACCTGAC----ATCTT', sequenceAlignment2='--CCAG-CGTCAA-CTT', similarityScore=6.25
#      8) In the output section the result is going to be saved in out.txt file
#-------------------------------------------------------------------------------

# ALGORITHM: 
#
#
#
#
#


import os
import sys
import argparse
import decimal

def studentInfo():
    pythonVersion = '2.7.6'
    
    student1FirstName = "Lauren"
    student1LastName = "Boyd"
    student1SectionNumber = "001"
    student1NetId = "lboyd002"
    
    student2FirstName = "Juan Paulo"
    student2LastName = "Patulandong"
    student2SectionNumber = "001"
    student2NetId = "jpatu001"
    
    info = 'Python version: ' + pythonVersion + '\n' + '\n'
    info = info + 'FirstName: ' + student1FirstName + '\n'
    info = info + 'LastName: ' + student1LastName + '\n'
    info = info + 'Section: ' + student1SectionNumber + '\n'
    info = info + 'NetId: ' + student1NetId + '\n' + '\n'
    
    info = info + 'FirstName: ' + student2FirstName + '\n'
    info = info + 'LastName: ' + student2LastName + '\n'
    info = info + 'Section: ' + student2SectionNumber + '\n'
    info = info + 'NetId: ' + student2NetId + '\n' + '\n'
    
    return info

def penalty(x,y,seq1,seq2): # Penalty Definitions
    if (seq1[x] == seq2[y]): 
    	return 1.0 #MATCH
    elif ( (seq1[x] == 'A' and seq2[y] == 'T') or (seq1[x] == 'T' and seq2[y] == 'A') ): 
	    return -0.15 #SUBS T-A 
    elif ((seq1[x] == 'C' and seq2[y] == 'G') or (seq1[x] == 'G' and seq2[y] == 'C')): 
	    return -0.15 #SUBS G-C 
    else: 
	    return -0.1 #SUBS for the rest


def biggest(deletion, insert, match):
    if (deletion > insert and deletion > match):
        Max = deletion
    elif(insert > deletion and insert > match):
        Max = insert
    else:
        Max = match
    return Max

def reversed(str):
    return str[::-1]

def DNASeqAlignment(DNASeq1,DNASeq2,outputPath):
    similarityScore = 0
    sequenceAlignment1 = ''
    sequenceAlignment2 = ''
    
    #########################################################################################
    # Compute new values for similarityScore and sequenceAlignment1 and sequenceAlignment2  #                                                                  #
    #########################################################################################
    
    #end of rows
    end1 = len(DNASeq1) + 1 
    #end of columns
    end2 = len(DNASeq2) + 1 
    
    #Matrix initialization
    c = {}
    #make first spot of matrix 0
    c[0,0] = 0
    #loop though the first row and set values
    for i in range(1,end1):
        c[i,0]= c[i-1,0] - 0.2
    #loop thought the first column and set values
    for j in range(1,end2):
        c[0,j] = c[0, j-1] - 0.2
    # loop through the rows
    for i in range(1,end1):
        #loop through the columns
        for j in range(1,end2):
            deletion = c[i-1,j] - 0.2
            insert = c[i,j-1] - 0.2
            match = c[i-1, j-1] + penalty(i-1, j-1, DNASeq1, DNASeq2)
            
            #use function to check which is the max score and place into matrix
            c[i,j] = biggest(deletion,insert,match)
    
    #printing matrix
    #for i in range(0, end1):
    #    for j in range(0, end2):
    #        print c[i,j],
    #    print '\n'
    
    #similarity score becomes the last entry in the matrix
    similarityScore = c[ end1-1, end2-1]
    
    #incrementors to loop through rows and columns
    i, j =0, 0
    #using the length of the strings to start from the end of the matrix
    endi = end1-1
    endj = end2-1
    #loop through until incrementors hit the end of the column or the row
    while( (i < end1 and endi >0) or (j < end2 and endj > 0)):
        #check if the spot above the index we are looking at is greater
        if(c[endi - 1, endj] > c[endi -1, endj -1] and c[endi - 1, endj] > c[endi, endj -1]):
            #if it is, then skip for the first DNA
            sequenceAlignment1 = sequenceAlignment1 + '_'
            sequenceAlignment2 = sequenceAlignment2 + DNASeq2[endj-1]
            #increment i
            if(i < end1):
                i = i + 1
        #check if spot to the left is greater
        elif(c[endi, endj-1] > c[endi -1, endj -1] and c[endi, endj-1] > c[endi-1, endj -1]):
            #if it is, then skip for the second DNA
            sequenceAlignment1 = sequenceAlignment1 + DNASeq1[endi-1]
            sequenceAlignment2 = sequenceAlignment2 + '_'
            #increment j
            if(j < end2):
                j = j +1
        else:
            sequenceAlignment1 = sequenceAlignment1 + DNASeq1[endi-1]
            sequenceAlignment2 = sequenceAlignment2 + DNASeq2[endj-1]
            if(i < end1 - 1):
                i = i + 1
            if(j < end2 -1):
                j = j +1
        #decrease the backwards index
        endi = endi - 1
        endj = endj - 1
        
    sequenceAlignment1 = reversed(sequenceAlignment1)
    sequenceAlignment2 = reversed(sequenceAlignment2)
    
    #################################  Output Section  ######################################
    result = "Similarity score: " + str(similarityScore) + '\n'
    result = result + "Sequence alignment1: " + sequenceAlignment1 + '\n'
    result = result + "Sequence alignment2: " + sequenceAlignment2 + '\n'
    writeToFile(outputPath,result)
    
def writeToFile(filePath, content):
    with open(filePath,'a') as file:
        file.writelines(content)

def readFile(filePath):
    logLines = ''
    with open(filePath,'r') as file:
        for logText in file:
            logLines = logLines + logText

    uniqueChars = ''.join(set(logLines))
    for ch in uniqueChars:
        if ch not in ('A','a','C','c','G','g','T','t'):
            logLines = logLines.replace(ch,'')
    logLines = logLines.upper()
    return logLines

def removeFile(filePath):
    if os.path.isfile(filePath):
        os.remove(filePath)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='DNA sequence alignment')
    parser.add_argument('DNASeq1FilePath', type=str, help='Path to the file that contains First DNA sequence')
    parser.add_argument('DNASeq2FilePath', type=str, help='Path to the file that contains Second DNA sequence')
    parser.add_argument('OutputFilePath', type=str, help='Path to the output file')
    args = parser.parse_args()
    DNASeq1 = readFile(args.DNASeq1FilePath)
    DNASeq2 = readFile(args.DNASeq2FilePath)
    outputPath = args.OutputFilePath
    removeFile(outputPath)
    writeToFile(outputPath,studentInfo())
    DNASeqAlignment(DNASeq1,DNASeq2,outputPath)
