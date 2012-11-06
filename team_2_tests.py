'''
Runs a number of tests on the classes and functions in our phylogeny generator
'''

from team_2_optimal_alignment_sensitive import *
from team_2_neighbor_joining import *
from readfasta import readfasta
from team_2_bootstrapping import getBootstrappedSequences

genomes = [ "AATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGTATTTTCG",
            "GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGTATTTTCGTCTGG",
            "GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGTATTTTCGTCTGG",
            "GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGTATTTTCGTCTGGGG",
            "GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGTATTTTCGTCTGGG",
            "GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGTATTTTCGTCTGGGGG",
            "GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGTATTTTCGTCTGG"]
print("Running optimal alignment . . .")
alignment = OptimalAlignment( genomes[0], genomes[1] )
print("Optimal alignment score: ", alignment.getScore(), "\n")
consensus = "".join(alignment.getConsensusSeq())
print("Consensus seq: ", consensus)
print("Consensus seq length: ", len(consensus))

simulatedOutputOfReadFASTA = [ ["0", "someID Some long description", "GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGTATTTTCG"],
                               ["1", "otherID Some long description", "GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGTATTTTCGTCTGG"],
                               ["2", "ID3 Some long description", "GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGTATTTTCGTCTGG"],
                               ["3", "4thID Some long description", "GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGTATTTTCGTCTGGG"],
                               ["4", "fiftheID Some long description", "GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGTATTTTCGTCTGG"]]
nodes = getNeighborJoiningSequences( simulatedOutputOfReadFASTA )
#distMatrix = [ [0, 1, 2, 3],
#               [1, 0, 4, 5],
#               [2, 4, 0, 6],
#               [3, 5, 6, 0], ]
distMatrix = [ [0, 1, 2, 3, 4],
               [1, 0, 5, 6, 7],
               [2, 5, 0, 8, 9],
               [3, 6, 8, 0, 10],
               [4, 7, 9, 10, 0] ]
winningNode = getNeighborJoiningPhylogeny(nodes, distMatrix)
print("Results of neighbor joining:")
print("\tGot:     ", winningNode.getIdentifier() )
print("\tExpected: (3, (1, (4, (0, 2))))")

dir = "/Users/s3cur3/Dropbox/school/CS 325 - Bioinformatics/PhylogeneticTreeBuilder/"
newDistMatrix = constructMatrixFromFile( dir + "alignments.txt",
                                         readfasta(dir + "mtDNA.fasta") )
print("\nResults of matrix construction for canonical: ")
for row in newDistMatrix:
    print(row)
print()

finalNode = doNeighborJoining(dir + "alignments.txt",
                              readfasta(dir + "mtDNA.fasta"))
print("\n\n", finalNode.getTreeFile())

simulatedOutputOfReadFASTA = [ ["0", "someID Some description", "ABCDEFGHIJKLMNOP"*10000],
                               ["1", "otherID Some description", "ABCDEFGHIJKLMNOP"*10000],
                               ["2", "ID3 Some description", "ABCDEFGHIJKLMNOP"*10000],
                               ["3", "4thID Some description", "ABCDEFGHIJKLMNOP"*10000],
                               ["4", "fiftheID Some description", "ABCDEFGHIJKLMNOP"*10000]]
print("\nResults of bootstrapping:")
for seq in getBootstrappedSequences(simulatedOutputOfReadFASTA):
    print(seq[0:100])

realOutputOfReadFASTA = readfasta(dir + "mtDNA.fasta")
print("\nBootstrapping on a real sequence (truncated)")
for seq in getBootstrappedSequences(realOutputOfReadFASTA):
    print(seq[0:100])