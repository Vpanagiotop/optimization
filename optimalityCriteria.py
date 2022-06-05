import numpy as np
def optimalityCriteria(x, dc, g):
    l1 = 0
    l2 = 1e9
    move = 0.2
    # reshape to perform vector operations
    xnew = np.zeros(len(x))
    while (l2 - l1) / (l1 + l2) > 1e-3:
        lmid = 0.5 * (l2 + l1)
        xnew[:]= np.maximum(0.0, np.maximum(x - move, np.minimum(1.0, np.minimum(x + move, x * np.sqrt(-dc/lmid)))))
        gt=g+np.sum((xnew-x))
        if gt>0 :
            l1=lmid
        else:
            l2=lmid
    return (xnew,gt)