class GridCell:
    def __init__(self, id, pos, volume):
        self.id = id
        self.position = pos
        self.vol = volume
    def getId(self):
        return self.id
    def pos(self):
        return self.position
    def volume(self):
        return self.vol
    
    
class RegularGrid:
    def __init__(self, dim, size, orig):
        self.dim = dim
        self.len = dim[0]*dim[1]*dim[2]
        self.size = size
        self.vol = size[0]*size[1]*size[2]
        self.orig = orig
        self.cells = []
        for i in range(dim[2]):
            for j in range(dim[1]):
                for k in range(dim[0]):
                    pos = self.posFromIdx((k,j,i))
                    idx = k + (j + i*dim[2])*dim[1]
                    self.cells.append( GridCell(idx, pos, self.vol) )

    def getCells(self):
        return self.cells
    
    def listIdxFromPos(self, p, sig):
        idx = self.idxFromPos(p)
        return idx[0] + idx[1]*self.dim[1] + idx[0]*self.dim[1]*self.dim[2] + sig*self.dim[1]*self.dim[2]*self.dim[3] # + idx[2]*self.dim[0]*self.dim[1]
    
    def idxFromPos(self, p):
        return ( (p[0]-self.orig[0])/self.size[0], (p[1]-self.orig[1])/self.size[1], (p[2]-self.orig[2])/self.size[2] )

    def posFromIdx(self, idx):
        return ( idx[0]*self.size[0]+self.orig[0], idx[1]*self.size[1]+self.orig[1], idx[2]*self.size[2]+self.orig[2])
        
    def divide(self, cell, a0, a1):
        pass
    
    def update(self, dt):
        pass