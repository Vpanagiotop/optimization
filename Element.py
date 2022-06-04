import numpy as np
class Element:
    def __init__(self, tag, nodeTag, node):
        self.tag = tag
        self.node_tag = nodeTag
        self.node = [node[i] for i in nodeTag]
        self.a, self.b, self.r, self.rho, self.m, self.lamda = getElementInfo(self.node)
        self.nominalStiffness = getNominalStiffness(self.r, self.rho, self.m, self.lamda)
        self.stiffness = self.nominalStiffness.copy()
        dof = np.array([], dtype=int)
        self.centerCoords = [self.node[0].x+self.a, self.node[0].y+self.b]
        for i in range(len(self.node)):
            dof = np.append(dof, self.node[i].dof)
        self.dof = dof

def getElementInfo(node):
    a = abs(node[0].x - node[1].x)/2
    b = abs(node[0].y - node[2].y)/2
    r = a/b
    rho = (1 - 0.3)/2
    m = round(3 * (1 + 0.3) / 2,5)
    lamda = round(3 * (1 - 3 *  0.3) / 2,5)
    return a, b, r, rho, m, lamda

def getNominalStiffness(r, rho, m, lamda):
    k = [
        4/r + 4 * rho * r,
        m,
        -4/r + 2 * rho * r,
        -lamda,
        -2/r - 2 * rho * r,
        -m,
        2/r - 4 * rho * r,
        lamda
    ]

    return 1/(12 * (1-0.3**2)) * np.array([
        [k[0], k[1], k[2], k[3], k[4], k[5], k[6], k[7]],
	    [k[1], k[0], k[7], k[6], k[5], k[4], k[3], k[2]],
	    [k[2], k[7], k[0], k[5], k[6], k[3], k[4], k[1]],
        [k[3], k[6], k[5], k[0], k[7], k[2], k[1], k[4]],
        [k[4], k[5], k[6], k[7], k[0], k[1], k[2], k[3]],
        [k[5], k[4], k[3], k[2], k[1], k[0], k[7], k[6]],
        [k[6], k[3], k[4], k[1], k[2], k[7], k[0], k[5]],
        [k[7], k[2], k[1], k[4], k[3], k[6], k[5], k[0]]
    ])
    