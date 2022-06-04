from Generator import Generator
import Plotter
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors

x = 50
model = Generator(length=x, height=x, elementsNumber={'x':x, 'y':x})
volfrac = 0.5
x=volfrac * np.ones(model.elementsNumber['x']*model.elementsNumber['y'],dtype=float)
xold=x.copy()
xPhys=x.copy()
for node in model.node:
    if node.x == 0:
        node.xFixed = True
    elif node.x == model.length and node.y == 0:
        node.yFixed = True
    if node.x == 0 and node.y == model.height:
        node.Py = -1
model.format()

model.modify(xPhys)
u = model.solveFE()
ce = np.array([])
for ele in model.element:
    ce = np.append(ce, (np.dot(u[ele.dof].reshape(1,8), ele.nominalStiffness) * u[ele.dof].reshape(1,8)).sum())
obj=((model.Emin + xPhys**model.pen * (model.Emax - model.Emin)) * ce).sum()
dc = (-model.pen * xPhys**(model.pen-1)*(model.Emax-model.Emin))*ce
Plotter.plotNodes(model.node)
xPhys = np.ones(len(model.element),dtype=float) * 0.5
plt.ion() # Ensure that redrawing is possible
fig,ax = plt.subplots()
im = ax.imshow(-xPhys.reshape((model.elementsNumber['x'],model.elementsNumber['y'])).T, cmap='gray',\
	interpolation='none',norm=colors.Normalize(vmin=-1,vmax=0))
fig.show()

xx =1