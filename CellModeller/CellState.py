from copy import deepcopy

class CellState(dict):
    # Don't show these attributes in gui
    excludeAttr = ['id', 'divideFlag', 'ends','cellAdh']
    
    """dot.notation access to dictionary attributes"""
    __getattr__ = dict.get
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__

    def __init__(self, cid):
        self.id = cid
        self.growthRate = 1.0
        self.color = [0.5,0.5,0.5]
        self.divideFlag = False
        self.cellAge = 0
        self.neighbours = []
        self.effGrowth = 0.0
        self.cellType = 0

    def __deepcopy__(self, memo=None):
        return CellState(deepcopy(dict(self), memo=memo))
