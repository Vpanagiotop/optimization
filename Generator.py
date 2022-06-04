import numpy as np
from scipy.sparse import coo_matrix, linalg
from Node import Node
from Element import Element
class Generator:
    def __init__(
        self, length, height, elementsNumber={'x': int, 'y': int},
        Emin:float = 1.e-9, Emax:float = 1, pen=3.0
    ):
        self.length = length
        self.height = height
        self.elementsNumber = elementsNumber
        self.grid = [elementsNumber['x'] + 1, elementsNumber['y'] + 1]
        self.elementWidth = length/elementsNumber['x']
        self.elementHeight = height/elementsNumber['y']
        self.Emin = Emin
        self.Emax = Emax
        self.pen = pen
        
        points = []
        xCoordinates = np.linspace(0, length, self.grid[0], dtype=float)
        yCoordinates = np.linspace(0, height, self.grid[1], dtype=float)
        for x in xCoordinates:
            for y in yCoordinates:
                points.append([x,y])
        self.node = assignNodes(points)
        self.dofs = assigndegreesOfFreedom(self.node)
        self.element = assignElements(elementsNumber, self.node)
        self.elementDofs = assignElementDofs(self)
        self.Kel, self.Kmatrix = getElasticStiffness(self)
        self.H, self.Hs = getFilterWeightMatrix(self)
    
    def format(self):
        self.fixedDofs, self.freeDofs = setBoundaryConditions(self)
        self.F = setNodalLoads(self)

    def modify(self, xPhys):
        getElementStiffness(self, xPhys)
        self.Kel, self.Kmatrix = getElasticStiffness(self)

    def solveFE(self):
        Kff = self.Kmatrix[self.freeDofs,:][:,self.freeDofs]
        u=np.zeros((len(self.dofs),1))
        u[self.freeDofs,0] = linalg.spsolve(Kff,self.F[self.freeDofs,0])
        return u

def assignNodes(points):
    node = []
    for i in range(len(points)):
        tag = i
        x, y = points[i][0], points[i][1]
        node.append(Node(tag, x, y)) 
    return node

def assigndegreesOfFreedom(node):
    dof = np.linspace(0,len(node)*2 - 1, len(node)*2, dtype=int)
    for i in range(len(node)):
        node[i].dof = [dof[2 * i], dof[2 * i + 1]]
    return dof

def assignElements(elementsNumber, node):
    element = []
    tag = 0
    for xElement in range(elementsNumber['x']):
        i = xElement * (elementsNumber['y'] + 1)
        for j in range(elementsNumber['y']):
            nodeTag = [i + j, i + elementsNumber['y'] + 1 + j,  i + elementsNumber['y'] + 2 + j, i + j + 1]
            element.append(Element(tag, nodeTag, node))
            tag = tag + 1
    return element

def assignElementDofs(self):
    elementDofs = np.zeros((len(self.element),8), dtype =int)
    for ele in self.element:
        elementDofs[ele.tag] = np.array([ele.dof])
    return elementDofs

def getElasticStiffness(self):
    Kel = np.zeros((self.dofs.size, self.dofs.size), float)
    for ele in self.element:
        for iEle in range(len(ele.dof)):
            i = int(ele.dof[iEle])
            for jEle in range(len(ele.dof)):
                j = int(ele.dof[jEle])
                Kel[i,j] = Kel[i,j] + float(ele.stiffness[iEle, jEle])
    Kmatrix = coo_matrix(Kel).tocsc()
    return Kel, Kmatrix

def setBoundaryConditions(self):
    fixedDofs = np.empty((0,0), dtype=int)
    for node in self.node:
        if node.xFixed == True:
           fixedDofs =  np.append(fixedDofs,node.dof[0])
        if node.yFixed == True:
           fixedDofs =  np.append(fixedDofs,node.dof[1])
    freeDofs = np.setdiff1d(self.dofs,fixedDofs)
    return fixedDofs, freeDofs

def setNodalLoads(self):
    F = np.zeros((len(self.dofs),1))
    for node in self.node:
        if node.Px != 0:
            F[node.dof[0]] = node.Px
        if node.Py!= 0:
            F[node.dof[1]] = node.Py
    return F

def getFilterWeightMatrix(self):
    H = np.zeros((len(self.element), len(self.element)))
    for iEle in self.element:
        i = iEle.tag
        xi = iEle.centerCoords[0]
        yi = iEle.centerCoords[1]
        minBound = max(0, i - 5 * self.elementsNumber['y'])
        maxBound = min(len(self.element), i + 5 * self.elementsNumber['y'])
        for jEle in self.element[minBound:maxBound]:
            j = jEle.tag
            xj = jEle.centerCoords[0]
            yj = jEle.centerCoords[1]
            H[i,j] = max(5 - np.sqrt((xj - xi)**2 + (yj - yi)**2), 0)
    H=coo_matrix(H).tocsc()
    Hs=H.sum(1)
    return H, Hs
            
def getElementStiffness(self, xPhys):
    for ele in self.element:
        ele.stiffness = ele.nominalStiffness * (self.Emin+(xPhys[ele.tag])**self.pen*(self.Emax-self.Emin))
