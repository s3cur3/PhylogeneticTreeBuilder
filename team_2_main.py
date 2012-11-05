'''
A program to build phylogenetic trees.

Written by Tyler Young, 1 November 2012
'''

from team_2_optimal_alignment import *
import argparse
import cProfile
from readfasta import readfasta
import time
from team_2_nj_node import *


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


def getNeighborJoiningSequences( inputFromReadFasta ):
    '''
    Returns a list of nodes which are representations of the sequences in the
    FASTA file.
    '''
    nodes = []
    for i in range(len(inputFromReadFasta)):
        nodes.append( NeighborJoiningNode(inputFromReadFasta[2], inputFromReadFasta[0], inputFromReadFasta[1] ) )

    return nodes


def doNeighborJoining(sequences, distMatrix):
    '''
    Performs neighbor joining on the input sequence, using data in the distance
    matrix.
    '''
    numSeqs = len(sequences)
    distFromAllOthers = [0] * numSeqs

    for i in range(numSeqs):
        for j in range(numSeqs):
            # This is u_i in the algorithm formulation
            distFromAllOthers[i] += distMatrix[i][j]

        distFromAllOthers[i] /= (numSeqs - 2) # divide that sum by n-2

    distPrime = [ [0] * numSeqs for i in numSeqs ]

    minDPrime = 100000000
    pairLoc = (-1,-1)
    pairWithLowestDPrime = (NeighborJoiningNode("?"), NeighborJoiningNode("?"))
    for i in range(numSeqs):
        for j in range(numSeqs):
            distPrime[i][j] = (numSeqs - 2) * distMatrix[i][j] \
                              - distFromAllOthers[i] - distFromAllOthers[j]
            if distPrime[i][j] < minDPrime:
                pairLoc = (i,j)
                pairWithLowestDPrime = (sequences[i], sequences[j])


    newNodeName = "(" + pairWithLowestDPrime[0].getIdentifier() + ", " \
                  + pairWithLowestDPrime[1].getIdentifier() + ")"
    newNode = NeighborJoiningNode( newNodeName )
    distToChild0 = (distMatrix[pairLoc[0]][pairLoc[1]]
                    + distFromAllOthers[pairLoc[0]]
                    - distFromAllOthers[pairLoc[1]])/2
    newNode.addChild( pairWithLowestDPrime[0], distToChild0 )
    distToChild1 = (distMatrix[pairLoc[0]][pairLoc[1]]
                    + distFromAllOthers[pairLoc[1]]
                    - distFromAllOthers[pairLoc[0]])/2
    newNode.addChild( pairWithLowestDPrime[1], distToChild1 )
    newNode.setParent( pairWithLowestDPrime[0].getParent() )

    # Modify the distance matrix: remove the 2 old nodes, add in the new node
    distToNewNode = [0]*len(sequences)
    for k in range(len(sequences)):
        distToNewNode[k] = (distMatrix[pairLoc[0]][k] + distMatrix[pairLoc[1]][k]
                            + distMatrix[pairLoc[0]][pairLoc[1]])/2

    # Remove the 2 old nodes from sequences

    # Add the new node to the sequences

    # Remove the 2 old nodes from the distance matrix

    # Add the new node to the distance matrix using the calculated distances









parser =  argparse.ArgumentParser(description="Prepare a dot plot for two "
                                              + "FASTA sequences")
parser.add_argument("file", type=str,
                    help="The path to the FASTA file "
                         + "containing all sequences to be compared.")
args = parser.parse_args()

allSequenceData = readfasta(args.file)

print(len(allSequenceData))


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

