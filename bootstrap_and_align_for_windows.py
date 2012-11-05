from team_2_optimal_alignment import *
import argparse
from readfasta import readfasta
import time
import random
import datetime

allSequenceData = readfasta("mtDNA.fasta")

while True:
    allSequenceData = readfasta(args.file)

    sequences = []
    genomeIdentifiers = []
    for sequenceData in allSequenceData:
        genomeIdentifiers.append(sequenceData[0])
        sequences.append(sequenceData[2])

    # Normalize all sequence lengths
    maxLength = 0
    for sequence in sequences:
        maxLength = max( maxLength, len(sequence) )
    for i in range(len(sequences)):
        # Add some number of dashes to the end
        sequences[i] += '-' * (maxLength - len(sequences[i]))

    # Fill up the new sequences with random selections from the old one
    bootstrappedSeqs = [ ['\n'] * maxLength for seq in sequences ]
    rangeOfOptions = range(maxLength)
    for i in rangeOfOptions:
        # Choose a column
        chosenColumn = random.choice(rangeOfOptions)
        # In each sequence, replace the next spot with the value of that column
        for seq in range(len(bootstrappedSeqs)):
            bootstrappedSeqs[seq][i] = sequences[seq][chosenColumn] #sequences[] is OOR

    # Put the new, random sequences back into the sequences variable
    for i in range(len(sequences)):
        sequences[i] = "".join(bootstrappedSeqs[i])


    # Write the sequences
    bootstrapID = "bootstrap" + str(int(time.time()) % 1351800000)
    fastaFileName = bootstrapID + ".fasta"
    fastaFile = open(fastaFileName, "w")
    for i in range(len(sequences)):
        fastaFile.write( ">" + allSequenceData[i][0] + "_" + bootstrapID + " "
                         + allSequenceData[i][1] + "\n" + sequences[i] + "\n\n" )
    fastaFile.write("\n")
    fastaFile.close()

    now = datetime.datetime.now()
    print(now.strftime("%Y-%m-%d %H:%M") + ": Wrote the bootstrap (randomized) sequences to disk. (File name "
          + fastaFileName + " in your current working directory.)" )


    fileName = "alignments_" + bootstrapID + ".txt"

    # Get the scores for all alignments of the bootstrapped sequences
    allAlignments = [ [OptimalAlignment() for x in range(len(sequences))]
                      for x in range(len(sequences))]
    allScores = [] # List of tuples: (score, position-in-allAlignments-tuple)
    theCount = 0
    for x in range(len(sequences)):
        for y in range(0, x):
            theAlignment = OptimalAlignment( sequences[x], sequences[y] )
            file = open(fileName, "a")
            file.write( "\nSequences compared: " + str(genomeIdentifiers[x]) + " and " + str(genomeIdentifiers[y]))
            file.write( "\nScore: " + str(theAlignment.getScore()) )
            file.write( "\nConsensus seq: " )
            file.write( "".join(theAlignment.getConsensusSeq()) )
            file.write( "\n" )
            file.close()
            allScores.append( theAlignment.getScore() )
            theCount += 1
            now = datetime.datetime.now()
            print(now.strftime("%Y-%m-%d %H:%M") + ": Finished a pass! Just " + str(91 - theCount) + " to go!")

    allScores.sort() # highest-scoring is first
    file = open(fileName, "a")
    file.write( str(allScores) )
    file.write( "\n" )
    file.close()

    print("\n\nFinished all comparisons!\n\n\n")
