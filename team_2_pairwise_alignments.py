'''
A method for performing all pairwise alignments of a set of FASTA sequences.

Tyler Young
'''

from team_2_optimal_alignment_sensitive import *
import datetime, math

def writePairwiseAlignmentFile( allSequenceData, destinationFileName="alignments.fasta" ):
    '''
    Writes a file with all pairwise alignment data. Note that this can take about 15
    minutes per pairwise comparison, for a total of nearly a full day for the full
    set of comparisons on 14 species.
    @param allSequenceData The output of readFasta() on your input FASTA file
    @param destinationFileName The file to which we should write the alignment data
                               If this file already exists, we will simply append
                               to it.
    '''
    now = datetime.datetime.now()
    print("\n", now.strftime("%Y-%m-%d %H:%M"), "Beginning all pairwise alignments.")

    sequences = []
    genomeIdentifiers = []
    for sequenceData in allSequenceData:
        genomeIdentifiers.append(sequenceData[0])
        sequences.append(sequenceData[2])

    theCount = 0
    for x in range(len(sequences)):
        for y in range(0, x):
            theAlignment = OptimalAlignment( sequences[x], sequences[y] )
            file = open(destinationFileName, "a")
            file.write( "\nSequences compared: " + str(genomeIdentifiers[x])
                        + " and " + str(genomeIdentifiers[y]))
            file.write( "\nScore: " + str(theAlignment.getScore()) )
            file.write( "\nConsensus seq: " )
            file.write( "".join(theAlignment.getConsensusSeq()) )
            file.write( "\n" )
            file.close()

            theCount += 1

            # Print the status
            totalComparisons = (len(sequences)*len(sequences))/2 - len(sequences)/2
            now = datetime.datetime.now()
            print(now.strftime("%Y-%m-%d %H:%M:%S"), "Finished an alignment. Just",
                  str(totalComparisons - theCount), "to go!")
            print( "[", "*" * int(math.ceil(theCount/totalComparisons*100)),
                   " " * int(math.ceil((totalComparisons-theCount)/totalComparisons*100)),
                   "]"  )
