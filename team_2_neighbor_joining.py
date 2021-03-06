'''
A set of tools for getting a phylogenetic tree from alignment data.
Use the doNeighborJoining() method as the preferred way of getting the phylogeny
directly from your combination of alignment and FASTA data.

Tyler Young
'''

from team_2_fasta_tools import countChanges
from team_2_optimal_alignment_sensitive import *
DEBUGGING = False
ADD_DISTANCES = True

def doNeighborJoining(alignmentData, fastaData):
    '''
    Reads the input alignment data file (generated by either the align_all.py or
    bootstrap_and_align.py programs) and the FASTA data, and generates the neighbor
    joining phylogeny.

    This method encapsulates everything else in the file.
    @param alignmentData A .txt file output by either the align_all.py or
                         bootstrap_and_align.py programs. Contains a series of sets
                         of 3 lines, which are the pairwise alignments of 2
                         sequences together with their scores. These look like:
                                Sequences compared: someID1 and someID2
                                Score: 12345
                                Consensus seq: ABCDEFG
    @param fastaData The output of readFasta() on the .fasta file corresponding to
                     this alignment data.
    @return The final "node", whose name, which looks something like
            ((species1, species2), species 3), represents the taxonomy returned by
            neighbor joining.
    '''
    distanceMatrix = constructMatrixFromFile( alignmentData, fastaData )
    sequenceList = getNeighborJoiningSequences( fastaData )
    return getNeighborJoiningPhylogeny( sequenceList, distanceMatrix )




def constructMatrixFromFile(alignmentData, fastaData):
    '''
    Reads the input file (generated by the team_2_pairwise_alignments.py program)
    into a distance matrix, which can be acted on by getNeighborJoiningPhylogeny().
    @param alignmentData A .txt file output by either the align_all.py or
                         bootstrap_and_align.py programs. Contains a series of sets
                         of 3 lines, which are the pairwise alignments of 2
                         sequences together with their scores. These look like:
                                Sequences compared: someID1 and someID2
                                Score: 12345
                                Consensus seq: ABCDEFG
    @param fastaData The output of readFasta() on the .fasta file corresponding to
                     this alignment data.
    @return A 2-D array suitable for input to doNeighborJoining(). Indices correspond
            to sequences in the list of sequences returned by
            constructNeighborJoiningSequencesFromFile(), and entries are the pairwise
            distances between each sequence.
    '''
    def getPosInMatrix(seqNames, fastaData, bootstrapID):
        '''
        @param seqNames The "Sequences compared" data from the alignment data file.
        @param fastaData The output of readFasta on your .fasta file
        @param bootstrapID Empty string if this is not a bootstrapping run. Else,
                           this should be a string like "bootstrap1234" which
                           identifies this run of the bootstrapper as unique from all
                           others
        @return An (x, y) tuple indicating the position of this comparison in the
                distance matrix
        '''
        name0, name1 = seqNames.split(" and ")
        x = -1
        y = -1
        for i in range(len(fastaData)):
            if bootstrapID != "":
                if fastaData[i][0] == name0 or fastaData[i][0].rstrip("_" + bootstrapID) == name0:
                    x = i
                if fastaData[i][0] == name1 or fastaData[i][0].rstrip("_" + bootstrapID) == name1:
                    y = i
            else:
                if fastaData[i][0] == name0:
                    x = i
                if fastaData[i][0] == name1:
                    y = i
        return (x, y)


    file = open(alignmentData, 'rU')
    bootstrappingID = ""
    if "bootstrap" in alignmentData:
        bootstrappingID = "bootstrap" \
                          + alignmentData.split("alignments_bootstrap")[1].rstrip(".txt")
    seqNames = ""
    score = 0
    seq = ""

    distMatrix = [ [0]*len(fastaData) for i in fastaData ]
    for line in file:
        line = line.rstrip()
        if line != "" and line[0] != "[":
            label, data = line.split(": ")
            if label == "Sequences compared":
                seqNames = data
            elif label == "Score":
                score = data
            else: # must be the consensus sequence
                seq = data

            x, y = getPosInMatrix(seqNames, fastaData, bootstrappingID)
            distMatrix[x][y] = int(score)

    # Now make it symmetric.
    for i in range(len(distMatrix)):
        for j in range(len(distMatrix)):
            if distMatrix[i][j] == 0 and distMatrix[j][i] != 0:
                distMatrix[i][j] = distMatrix[j][i]

    return distMatrix



def getNeighborJoiningSequences(fastaData):
    '''
    Takes FASTA data and transforms it into a list of nodes for the neighbor joining
    algorithm, as required by getNeighborJoiningPhylogeny().
    @param fastaData The output of readFasta() on your FASTA file.
    @return a list of nodes which are representations of the sequences in the
            FASTA file.
    '''
    nodes = []
    for fastaData in fastaData:
        nodes.append( NeighborJoiningNode(fastaData[0], fastaData[2], fastaData[1]) )
    return nodes


def getNeighborJoiningPhylogeny(sequences, distMatrix):
    '''
    Performs neighbor joining on the input sequence, using data in the distance
    matrix.
    @param sequences A list of NeighborJoiningNode returned by the
                     getNeighborJoiningSequences() method.
    @param distMatrix A 2-D array whose indices correspond to sequences in the list
                      of sequences and whose entries are the pairwise distances
                      between each sequence.
    @returnType NeighborJoiningNode
    @return The final "node", whose name, which looks something like
            ((species1, species2), species 3), represents the taxonomy returned by
            neighbor joining.
    '''
    if DEBUGGING:
        print("\nTop of getNeighborJoiningPhylogeny()")
        for row in distMatrix:
            print(row)
        print()

    numSeqs = len(sequences)
    distFromAllOthers = [0] * numSeqs

    # Find the distance from each sequence to all others
    for i in range(numSeqs):
        for j in range(numSeqs):
            # This is u_i in the algorithm formulation
            distFromAllOthers[i] += distMatrix[i][j]

    if DEBUGGING:
        print("\nDistance of each sequence to all others")
        print(distFromAllOthers, "\n")

    distPrime = [ [0] * numSeqs for i in range(numSeqs) ]

    # Construct the Distance_Prime matrix
    minDPrime = 100000000
    maxDPrime = -1000000000
    maxLoc = (-1,-1)
    pairWithMaxDPrime = (NeighborJoiningNode("?"), NeighborJoiningNode("?"))
    pairLoc = (-1,-1)
    pairWithLowestDPrime = (NeighborJoiningNode("?"), NeighborJoiningNode("?"))
    for i in range(numSeqs):
        for j in range(numSeqs):
            distPrime[i][j] = ((numSeqs - 2)*distMatrix[i][j]) - \
                              (distFromAllOthers[i] + distFromAllOthers[j])
            if distPrime[i][j] < minDPrime and i != j:
                minDPrime = distPrime[i][j]
                pairLoc = (i,j)
                pairWithLowestDPrime = (sequences[i], sequences[j])
            if distPrime[i][j] > maxDPrime and i != j:
                maxDPrime = distPrime[i][j]
                maxLoc = (i,j)
                pairWithMaxDPrime = (sequences[i], sequences[j])

    if DEBUGGING:
        print("\nDist prime:")
        for row in distPrime:
            print(row)
        print( "Lowest dprime at: ", pairLoc)
        print("\t(",pairWithLowestDPrime[0].getReadableIdentifier(),", ",
              pairWithLowestDPrime[1].getReadableIdentifier(),")\n",sep="")
        print( "Max dprime at: ", maxLoc)
        print("\t(",pairWithMaxDPrime[0].getReadableIdentifier(),", ",
              pairWithMaxDPrime[1].getReadableIdentifier(),")\n",sep="")

    # Get the distance to the new node's children in years
    if ADD_DISTANCES:
        print("Getting distances between ", pairWithMaxDPrime[0].getIdentifier(), "and",pairWithMaxDPrime[1].getIdentifier() )
        CHANGES_PER_YEAR = 1.7E-8
        consensusSeq = OptimalAlignment( pairWithMaxDPrime[0].getSequence(),
                                         pairWithMaxDPrime[1].getSequence() ).getConsensusSeq()
        delta = countChanges( consensusSeq )
        print("Num changes: ", delta)
        distInYears = (delta*(1/CHANGES_PER_YEAR))/len(consensusSeq)
        distanceString = ":" + str(distInYears)
        print("Optimal alignment begins", consensusSeq[0:100])
        print("Dist in years is ",distInYears)

    else:
        distanceString = ""
        consensusSeq = ""


    # Add a new node to the tree
    newNodeName = "(" + pairWithMaxDPrime[0].getReadableIdentifier() + distanceString + ", "\
                  + pairWithMaxDPrime[1].getReadableIdentifier() + distanceString + ")"
    newNode = NeighborJoiningNode( newNodeName )
    distToChild0 = (distMatrix[maxLoc[0]][maxLoc[1]]
                    + distFromAllOthers[maxLoc[0]]
                    - distFromAllOthers[maxLoc[1]])/2
    newNode.addChild( pairWithMaxDPrime[0], distToChild0 )
    distToChild1 = (distMatrix[maxLoc[0]][maxLoc[1]]
                    + distFromAllOthers[maxLoc[1]]
                    - distFromAllOthers[maxLoc[0]])/2
    newNode.addChild( pairWithMaxDPrime[1], distToChild1 )
    newNode.setParent( pairWithMaxDPrime[0].getParent() )
    newNode.setConsensusSequence(consensusSeq)

    if DEBUGGING:
        print("New node's distance to its children: ", distToChild0, ", ", distToChild1)

    # Modify the distance matrix: remove the 2 old nodes, add in the new node
    # Find the distance to the new node
    distToNewNode = [0]*len(sequences)
    for k in range(len(sequences)):
        distToNewNode[k] = (distMatrix[maxLoc[0]][k] + distMatrix[maxLoc[1]][k]
                            - distMatrix[maxLoc[0]][maxLoc[1]])/2

    if DEBUGGING:
        print("\n", distToNewNode, "\n")

    # Remove the 2 old nodes from sequences
    sequences.remove(pairWithMaxDPrime[0])
    sequences.remove(pairWithMaxDPrime[1])

    # Add the new node to the sequences
    sequences.append(newNode)

    if len(sequences) == 1:
        return sequences[0]
    else:
        # Add the new node to the distance matrix using the calculated distances
        distMatrix.append([0] * len(distMatrix[0]))
        for row in distMatrix: # for each old sequence in the matrix
            row.append(distToNewNode[i])
        for i in range(len(distMatrix[-1])): # for each spot in the new row
            distMatrix[-1][i] = distMatrix[i][-1]


        # Remove the 2 old nodes from the distance matrix
        sortedLocs = sorted(maxLoc)
        del distMatrix[sortedLocs[1]] # Order matters! (modifying list's indices)
        del distMatrix[sortedLocs[0]]
        for i in range(len(distMatrix)):
            del distMatrix[i][sortedLocs[1]]
            del distMatrix[i][sortedLocs[0]]

        if DEBUGGING:
            for row in distMatrix:
                print(row)
            print()
        return getNeighborJoiningPhylogeny(sequences, distMatrix)

class NeighborJoiningNode:
    '''
    A class to represent nodes in a tree for neighbor-joining.
    '''
    def __init__(self, identifier, fastaSequence="", fullDescription="" ):
        self.seq = fastaSequence
        self.id = identifier
        self.desc = fullDescription
        self.children = []
        self.parent = None

    def setParent(self, node):
        '''
        Sets this node's parent to be a particular node
        '''
        self.parent = node

    def addChild(self, node, distToChild):
        '''
        Adds a node to this node's list of children
        '''
        self.children.append( (node, distToChild) )

    def getSequence(self):
        return self.seq

    def getIdentifier(self):
        return self.id

    def getReadableIdentifier(self):
        if "29690980|gb|AY195787.1" in self.id:
            return "H. sapiens Navaho"
        elif "13272920|gb|AF346989.1" in self.id:
            return "H. sapiens Japan"
        elif "13272556|gb|AF346963.1" in self.id:
            return "H. sapiens Aborigine Australia"
        elif "13272836|gb|AF346983.1" in self.id:
            return "H. sapiens German"
        elif "13273130|gb|AF347004.1" in self.id:
            return "H. sapiens PNG Highlands"
        elif "119394815|gb|EF060354.1" in self.id:
            return "H. sapiens Bedouin"
        elif "32348589|gb|AY289094.1" in self.id:
            return "H. sapiens Samoa"
        elif "13272738|gb|AF346976.1" in self.id:
            return "H. sapiens Effik Africa"
        elif "13272878|gb|AF346986.1" in self.id:
            return "H. sapiens Ibo Africa"
        elif "48596203|gb|AY195783.2" in self.id:
            return "H. sapiens San Africa"
        elif "13273270|gb|AF347014.1" in self.id:
            return "H. sapiens Yoruba Africa"
        elif "13273074|gb|AF347000.1" in self.id:
            return "H. sapiens Mkambe Africa"
        elif "253947317|emb|FM865409.1" in self.id:
            return "Neanderthal El sidron"
        elif "643679|dbj|D38113.1|CHPMTB" in self.id:
            return "Chimpanzee"
        elif "5835149|ref|NC_001645.1" in self.id:
            return "Gorilla"
        else:
            return self.getIdentifier()

    def getDescription(self):
        return self.desc

    def getParent(self):
        '''
        Returns the parent of this node
        '''
        return self.parent

    def getChildren(self):
        '''
        Returns a list of tuples; each tuple is a (node, distance) pair indicating
        what node is a child of this one and the distance from this node to the child
        '''
        return self.children

    def getTreeFile(self):
        return self.getIdentifier() + ";"

    def setConsensusSequence(self, consensusSequence):
        self.seq = consensusSequence