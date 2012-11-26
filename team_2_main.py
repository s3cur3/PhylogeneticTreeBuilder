'''
A program to build phylogenetic trees.

Written by Tyler Young, 1 November 2012
'''

import argparse
import time
import os
import math
from readfasta import readfasta
from team_2_neighbor_joining import doNeighborJoining
from team_2_pairwise_alignments import writePairwiseAlignmentFile
from team_2_bootstrapping import bootstrapAlignmentsToFile


def clearScreen():
    print(os.system('clear'), chr(13), "  ", chr(13),)


# Process the program arguments
parser =  argparse.ArgumentParser(description="Takes a set of FASTA sequences and "
                                              + "their pairwise alignments, and "
                                              + "bootstraps them to generate trees "
                                              + "with confidence ratings.")
parser.add_argument("fastaFile", type=str,
                    help="The path to the FASTA file "
                         + "containing all sequences to be compared.")
parser.add_argument("-p", "--pairwiseComparisonFile", type=str,
                    required=False, default="",
                    help="The path to the alignment file "
                         + "containing all pairwise alignments and their scores. "
                         + "If this is not present, we will generate all pairwise "
                         + "comparisons (may take a very long time!!).")
parser.add_argument("-b", "--bootstrap", type=int, default=1000,
                    help="The number of bootstrapping iterations to perform. "
                         + "Reasonable values range from one to ten thousand. "
                           "These are computationally cheap (seconds per round)." )
args = parser.parse_args()

allTrees = {} # Simply counts the occurrence of each tree
canonicalFasta = readfasta(args.fastaFile)
canonicalComparisons = args.pairwiseComparisonFile

# Run all pairwise comparisons if necessary
# This will give us the "canonical comparison" file
if (canonicalComparisons == "") or (not os.path.exists(canonicalComparisons)):
    print("Pairwise comparison filing missing or absent. Running all comparisons...")
    if canonicalComparisons == "":
        fastaName = (args.fastaFile).rstrip(".fasta")
        canonicalComparisons = fastaName + ".txt"
    writePairwiseAlignmentFile( canonicalFasta, canonicalComparisons )

# Build the "canonical" phylogenetic tree and add it to the collection of all trees
canonicalFinalNode = doNeighborJoining(canonicalComparisons, canonicalFasta)
allTrees[canonicalFinalNode.getTreeFile()] = 1

# Build a phylogenetic tree many times (do bootstrapping!)
for i in range(args.bootstrap):
    # Read in the pairwise comparison file, and create a new, randomized one
    bootstrapID = "bootstrap" + str(int(time.time()) % 1351800000)
    bootstrapAlignmentsToFile( canonicalComparisons,
                               "bootstrap/alignment_" + bootstrapID + ".txt" )

    finalNode = doNeighborJoining(canonicalComparisons,
                                   canonicalFasta)

    if finalNode.getTreeFile() in allTrees.keys():
        allTrees[finalNode.getTreeFile()] += 1
    else:
        allTrees[finalNode.getTreeFile()] = 1
    clearScreen()
    print("Performing bootstrapping . . .")
    print( "[", "*" * int(math.ceil(i/args.bootstrap*100)),
                " " * int(math.ceil((args.bootstrap-i)/args.bootstrap*100)),
           "]"  )

for key in allTrees.keys():
    print("Tree:")
    print("\t", key)
    print("Frequency:")
    print("\t", allTrees[key])
