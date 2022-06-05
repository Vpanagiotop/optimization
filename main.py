from Generator import Generator
import Plotter
from optimalityCriteria import optimalityCriteria
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors


x = 100
model = Generator(length=100, height=50, elementsNumber={'x':100, 'y':50})
volfrac = 0.5
g = 0
x=volfrac * np.ones(model.elementsNumber['x'] * model.elementsNumber['y'],dtype=float)
xOld=x.copy()
xPhys=x.copy()
for node in model.node:
    if node.x == 0:
        node.xFixed = True
    elif node.x == model.length and node.y == model.height:
        node.yFixed = True
    if node.x == 0 and node.y == 0:
        node.Py = -1
model.format()
loop = 0
change = 1
plt.ion() # Ensure that redrawing is possible
fig,ax = plt.subplots()
im = ax.imshow(-xPhys.reshape((model.elementsNumber['x'],model.elementsNumber['y'])).T, cmap='gray',\
	interpolation='none',norm=colors.Normalize(vmin=-1,vmax=0))
fig.show()
while change>0.01 and loop<2000:
    loop=loop+1
    model.modify(xPhys)
    u = model.solveFE()
    ce = np.array([])
    for ele in model.element:
        ce = np.append(ce, (np.dot(u[ele.dof].reshape(1,8), ele.nominalStiffness) * u[ele.dof].reshape(1,8)).sum())
    obj=((model.Emin + xPhys**model.pen * (model.Emax - model.Emin)) * ce).sum()
    dc = (-model.pen * xPhys**(model.pen-1)*(model.Emax-model.Emin))*ce
    dc[:] = np.asarray((model.H*(x*dc))[np.newaxis].T/model.Hs)[:,0] / np.maximum(0.001,x)
    xOld[:] = x
    (x[:],g) = optimalityCriteria(x, dc, g)
    xPhys = x
    change = np.linalg.norm(x.reshape(len(x),1) - xOld.reshape(len(x),1),np.inf)
    im.set_array(-xPhys.reshape((model.elementsNumber['x'],model.elementsNumber['y'])).T)
    fig.canvas.draw()
    fig.canvas.flush_events()
    plt.pause(0.1)
    # Write iteration history to screen (req. Python 2.6 or newer)
    print("it.: {0} , obj.: {1:.3f} Vol.: {2:.3f}, ch.: {3:.3f}".format(\
                loop,obj,(g+volfrac*len(x))/(len(x)),change))
Plotter.plotNodes(model.node)



xx =1