class Node:
    def __init__(
        self, tag,
        x:int, y:int,
        xFixed:bool = False,
        yFixed:bool = False,
        Px: int = 0, Py: int = 0
    ):
        self.tag = tag
        self.x = x
        self.y = y
        self.xFixed = xFixed
        self.yFixed = yFixed
        self.Px = Px
        self.Py = Py