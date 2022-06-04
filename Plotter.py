import matplotlib.pyplot as plt
import matplotlib.patches as mpatch

def plotNodes(nodes):
    fig2, ax2 = plt.subplots()
    for node in nodes:
        if node.xFixed == True or node.yFixed == True:
            ax2.scatter(node.x, node.y, color='b', s=10)
        # if node.y == 0:
        #     ax2.annotate(node.tag, (node.x, node.y))