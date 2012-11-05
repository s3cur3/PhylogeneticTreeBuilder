'''
A class to represent nodes in a tree for neighbor-joining.
'''

class NeighborJoiningNode:
    def __init__(self, identifier ):
        self.seq = ""
        self.id = identifier
        self.desc = identifier
        self.children = []

    def __init__(self, fastaSequence, identifier, fullDescription ):
        self.seq = fastaSequence
        self.id = identifier
        self.desc = fullDescription
        self.children = []

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

