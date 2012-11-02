'''
A program to build phylogenetic trees.

Written by Tyler Young, 1 November 2012
'''

from team_2_optimal_alignment import *
import argparse
import cProfile
from readfasta import readfasta
import time


def getMultipleSequenceAlignment( inputListOfGenomes ):
    '''
    Uses clustering alignment to build a multiple sequence alignment of the
    input sequences
    @param inputListOfGenomes A list whose elements are FASTA gene sequences from
                              (presumably) different species
    @return ??
    '''
    # Create our own, mutable copy of the list of genomes
    sequences = list(inputListOfGenomes)
    while len(sequences) > 1:
        # Get the scores for all alignments
        allAlignments = [ [OptimalAlignment() for x in range(len(sequences))]
                          for x in range(len(sequences))]
        allScores = [] # List of tuples: (score, position-in-allAlignments-tuple)
        for x in range(len(sequences)):
            for y in range(0, x):
                allAlignments[x][y] = OptimalAlignment( sequences[x], sequences[y] )
                allScores.append( ( allAlignments[x][y].getScore(), (x,y) ) )

        sortedScores = sorted( allScores, key=lambda scoreTuple: scoreTuple[0] )
        sortedScores.reverse() # highest-scoring is first

        # Choose pairwise alignments of sequences and fix those
        # Start with the best alignment
        bestAlignmentCoordinates = ( sortedScores[0][1] )
        bestAlignment = allAlignments[bestAlignmentCoordinates[0]][bestAlignmentCoordinates[1]]
        newSequenceList = [ bestAlignment.getConsensusSeq(), ]

        # To determine the order of sequence addition:
        # 1. Build a best-guess phylogenetic tree
        # 2. Use that tree to determine the order to add successive sequences to the alignment





parser =  argparse.ArgumentParser(description="Prepare a dot plot for two "
                                              + "FASTA sequences")
parser.add_argument("file", type=str,
                    help="The path to the FASTA file "
                         + "containing all sequences to be compared.")
args = parser.parse_args()

allSequenceData = readfasta(args.file)


genomes = [ "GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGTATTTTCG",
            "GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGTATTTTCGTCTGG",
            "GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGTATTTTCGTCTGG",
            "GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGTATTTTCGTCTGGGG",
            "GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGTATTTTCGTCTGGG",
            "GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGTATTTTCGTCTGGGGG",
            "GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGTATTTTCGTCTGG"]
print(time.time())
alignment = OptimalAlignment( allSequenceData[0][2], allSequenceData[1][2] )
print("Score: ", alignment.getScore())
print(time.time())


#print(alignment.getConsensusSeq())

#print()

#getMultipleSequenceAlignment(genomes)

