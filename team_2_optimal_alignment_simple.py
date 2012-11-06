'''
Functions for performing optimal alignment of 2 sequences
Tyler Young
1 November 2012
'''

GAP_PENALTY = -7
MISMATCH = -2
MATCH = 1

from team_2_fasta_tools import getConsensusLetter
import array

DEFAULT_SCORE = -1000

def scoreConsensusSequence( consensusSequence, paddingChar='' ):
    '''
    Examines the consensus sequence and deduces what alignment must have produced
    it. Uses that information to generate a new score.
    @param consensusSequence The (assumedly randomized) consensus sequence
    @param paddingChar A padding character that might appear in the sequence. If we
                       see this character, we will ignore this location.
    @return The score of the consensus sequence
    '''
    score = 0
    for letter in consensusSequence:
        if letter in "ATCGU":
            score += MATCH
        elif letter in "X-":
            score += GAP_PENALTY
        elif letter == paddingChar:
            score += 0
        else: # was a mismatch
            score += MISMATCH
    return score

class OptimalAlignment:
    def __init__(self, seq0='', seq1=''):
        self.bestScore = DEFAULT_SCORE
        self.seq0 = seq0
        self.seq1 = seq1
        self.consensusSeq = []

        if seq0 != '' and seq1 != '':
            self.optimallyAlign(seq0, seq1)

    def getScore(self):
        return self.bestScore

    def getConsensusSeq(self):
        return self.consensusSeq

    def computeConsensusSeq(self, matrix):
        '''
        Calculates a consensus sequence from the most recent alignment
        Requires that you previously either initialized this object with two
        sequences or you called optimallyAlign().
        @return A list of characters representing the consensus sequence in
                FASTA format
        '''
        width = len(self.seq0)
        height = len(self.seq1)

        revConsensusSeq = [] # build the consensus sequence up in reverse
        x = width
        y = height
        for i in range( max(width, height) ):
            if matrix[x][y] == matrix[x-1][y-1] + MATCH \
                and self.seq0[x-1] == self.seq1[y-1]:
                # Match here
                revConsensusSeq.append(self.seq0[x-1])
            elif matrix[x][y] == matrix[x-1][y-1] + MISMATCH:
                # Mismatch
                revConsensusSeq.append(
                    getConsensusLetter(self.seq0[x-1], self.seq1[y-1]) )
            elif matrix[x][y] == matrix[x-1][y]:
                # Gap in y seq
                y += 1
                revConsensusSeq.append('X')
            else:
                # Gap in x seq
                x += 1
                revConsensusSeq.append('X')
            x -= 1
            y -= 1
        revConsensusSeq.reverse()
        self.consensusSeq = revConsensusSeq


    def optimallyAlign(self, seq0, seq1):
        '''
        Optimally aligns two sequences
        @params seq0, seq1 The two FASTA sequences to align
        @return The score of the optimal alignment
        '''
        if seq0 == self.seq0 and seq1 == self.seq1 \
            and self.bestScore != DEFAULT_SCORE:
            # Already calculated the score for these sequences!
            return self.getScore()
        else:
            self.seq0 = seq0
            self.seq1 = seq1

        width = len(seq0) + 1
        height = len(seq1) + 1
        # Using plain arrays sacrifices some speed in exchange for using
        # about 11% of the memory of a standard list
        matrix = [ array.array('i', range(height)) for x in range(width)]

        # Set up top row
        for x in range(width):
            matrix[x][0] = x*GAP_PENALTY

        # Set up the leftmost column
        for x in range(1, height):
            matrix[0][x] = matrix[0][x-1] + GAP_PENALTY

        for x in range(1, width):
            for y in range(1, height):
                if seq0[x-1] == seq1[y-1]:
                    scoreOfDiagonal = matrix[x-1][y-1] + MATCH
                else:
                    scoreOfDiagonal = matrix[x-1][y-1] + MISMATCH

                matrix[x][y] = max( matrix[x-1][y] + GAP_PENALTY,
                                    matrix[x][y-1] + GAP_PENALTY,
                                    scoreOfDiagonal )


        self.bestScore = matrix[width-1][height-1]
        self.computeConsensusSeq( matrix )
        # let's make it really clear to the garbage collector
        # that we no longer need that gigabyte of memory!
        matrix = []

    def dumpOptimalAlignment(self, seq0, seq1, matrix):
        '''
        Outputs the optimal alignment matrix to the screen for debugging
        '''
        print("      ", end='')
        for x in range(len(seq1)):
            print(seq1[x], end='  ')
        print( "\n ", end='' )
        for x in range(len(matrix)):
            # Print the label
            print(seq0[x-1], end = ' ') if x > 0 else print( ' ', end = '')

            for y in range(len(matrix[x])):
                prefix = '' if matrix[x][y] < 0 else ' '
                print(prefix, matrix[x][y], sep='', end=' ')
            print()
        print()