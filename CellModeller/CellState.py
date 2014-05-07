
class CellState:
    # Don't show these attributes in gui
    excludeAttr = ['divideFlag']

    def __init__(self, cid):
        self.id = cid
        self.growthRate = 1.0
        self.color = [0.5,0.5,0.5]
        self.divideFlag = False

