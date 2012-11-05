from team_2_optimal_alignment import *
import argparse
import cProfile
from readfasta import readfasta
import time


allSequenceData = readfasta("mtDNA.fasta")

sequences = []
genomeIdentifiers = []
for sequenceData in allSequenceData:
    genomeIdentifiers.append(sequenceData[0])
    sequences.append(sequenceData[2])

# Get the scores for all alignments
allAlignments = [ [OptimalAlignment() for x in range(len(sequences))]
                  for x in range(len(sequences))]
allScores = [] # List of tuples: (score, position-in-allAlignments-tuple)
for x in range(len(sequences)):
    for y in range(0, x):
        theAlignment = OptimalAlignment( sequences[x], sequences[y] )
        file = open("alignments.fasta", "a")
        file.write( "\nSequences compared: " + str(genomeIdentifiers[x]) + " and " + str(genomeIdentifiers[y]))
        file.write( "\nScore: " + str(theAlignment.getScore()) )
        file.write( "\nConsensus seq: " )
        file.write( "".join(theAlignment.getConsensusSeq()) )
        file.write( "\n" )
        file.close()
        allScores.append( theAlignment.getScore() )

allScores.sort() # highest-scoring is first
file = open("alignments.fasta", "au")
file.write( str(allScores) )
file.write( "\n" )
file.close()
