'''
A method to do bootstrapping (getting a random sequence from some starting sequences)

Tyler Young
6 November
'''
import random
from collections import namedtuple
from team_2_optimal_alignment_simple import scoreConsensusSequence

PADDING_CHAR = "?"

def bootstrapAlignmentsToFile( alignmentFileName, destinationFileName ):
    '''
    @param alignmentFileName The name of the alignment file that we should randomize
    @param destinationFileName The name of the file to which we should write our
                               randomized (bootstrapped) scoring/alignment data
    '''
    # Create a struct to hold the sequence data, with fields name and sequence
    SequenceData = namedtuple("SequenceData", "name consensusSequence")

    originalAlignmentFile = open(alignmentFileName, 'rU')
    seqNames = ""
    score = 0
    seq = ""

    allSequenceData = [ ]
    for line in originalAlignmentFile:
        line = line.rstrip()
        if line != "" and line[0] != "[":
            label, data = line.split(": ")
            if label == "Sequences compared":
                seqNames = data
            elif label == "Score":
                score = data
            else: # must be the consensus sequence
                seq = data

            if seqNames != "" and seq != "" and score != 0:
                allSequenceData.append( SequenceData(name=seqNames,
                                                     consensusSequence=seq) )
                # Reset the variables
                seqNames = ""
                seq = ""
                score = 0

    # Now that we have all the consensus sequences, perform bootstrapping
    simulatedOutputFromReadFASTA = []
    for sequence in allSequenceData:
        simulatedOutputFromReadFASTA.append( [sequence.name, sequence.name,
                                              sequence.consensusSequence] )
    randomizedSeqs = getBootstrappedSequences( simulatedOutputFromReadFASTA )

    # Write the randomized sequences (and their scores) to disk
    bootstrapAlignmentFile = open(destinationFileName, 'w')
    for i in range(len(simulatedOutputFromReadFASTA)):
        bootstrappedSeq = "".join(randomizedSeqs[i])
        consensusScore = scoreConsensusSequence(bootstrappedSeq, PADDING_CHAR)

        # Write the new alignment data to the file
        bootstrapAlignmentFile.write("Sequences compared: ")
        bootstrapAlignmentFile.write(simulatedOutputFromReadFASTA[i][0])
        bootstrapAlignmentFile.write("\nScore: ")
        bootstrapAlignmentFile.write( str(consensusScore) )
        bootstrapAlignmentFile.write("\nConsensus seq: ")
        bootstrapAlignmentFile.write( bootstrappedSeq )
        bootstrapAlignmentFile.write("\n\n")

    bootstrapAlignmentFile.close()




def getBootstrappedSequences( allSequenceData ):
    '''
    @param allSequenceData The output of the readFasta() method on your input file
    @return A randomized (bootstrapped) set of sequences taken from the original data
    '''
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
        # Add some number of "ignore-me" characters to the end
        sequences[i] += PADDING_CHAR * (maxLength - len(sequences[i]))

    # Fill up the new sequences with random selections from the old one
    bootstrappedSeqs = [['\n'] * maxLength for seq in sequences]
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

    return sequences