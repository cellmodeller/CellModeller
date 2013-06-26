%module CGSBacteriumExt
%{
#include "CGSBacterium/Solver.h"
%}

%inline %{typedef struct
{
	double x1,y1,z1;
	double x2,y2,z2;
} CellEnds;
%}

%feature("shadow") Cell::ends() %{
def ends(self,*args):
    e = $action(self,*args)
	v = ([e.x1, e.y1, e.z1], [e.x2, e.y2, e.z2])
	return v
%}


%extend Cell {
	PyObject * pos()
	{
		btVector3 p = self->position();
		btScalar * pval(p);
		PyObject * rval = PyList_New(3);
		PyObject *robj;
		for (int i=0; i<3; ++i)
		{
			robj = SWIG_From_double(static_cast< double >(pval[i]));
			PyList_SetItem(rval, i, robj);
		}
		return rval;
	}
};


%extend Cell {
	CellEnds ends()
	{
		btVector3 x1(-1,0,0); btVector3 x2(1,0,0);
		CellEnds ends;
		btTransform& t = self->transform;
		double l = self->length();//+self->radius();
		btVector3 p1 = t*(0.5*l*x1); 
		btVector3 p2 = t*(0.5*l*x2);
		btScalar * p1val(p1); btScalar * p2val(p2);
		ends.x1 = p1val[0]; ends.y1 = p1val[1]; ends.z1 = p1val[2];
		ends.x2 = p2val[0]; ends.y2 = p2val[1]; ends.z2 = p2val[2];
		return ends;
	}
};



%extend Cell {
	PyObject * contacts()
	{
		PyObject * rval = PyList_New(self->ctMap.size());
		PyObject *robj;
		//for (int i=0; i<self->nnbs; ++i)
		int i=0;
		for (std::map<int,Contact>::iterator it=self->ctMap.begin(); it!=self->ctMap.end(); ++it)
		{
			Contact& ct = it->second;
			double color[3];
			if (ct.active)
			{
				color[0]=1.0; color[1]=0.0; color[2]=0.0;
			}else{
				color[0]=0.0; color[1]=1.0; color[2]=0.0;
			}
			btVector3 pA = ct.pointA + self->transform.getOrigin();
			btVector3 pB = pA + ct.normal*10.0;
			btScalar * pAval(pA); btScalar * pBval(pB);
			PyObject * cpair = PyList_New(9);
			for (int j=0; j<3; ++j)
			{
				robj = SWIG_From_double(static_cast< double >(pAval[j]));
				PyList_SetItem(cpair, j, robj);
				robj = SWIG_From_double(static_cast< double >(pBval[j])); 
				PyList_SetItem(cpair, j+3, robj);
				robj = SWIG_From_double(color[j]);
				PyList_SetItem(cpair, j+6, robj);
			}
			PyList_SetItem(rval, i, cpair);
			++i;
		}
		return rval;
	}
};

%extend Cell {
	void setGrowthRate(double r)
	{
		self->lenVel = r;
	}
};
class Cell 
{

public:


	Cell(btVector3& p);
	btVector3& getPos();
	void updateInertiaTensor() ;

	// Interface for cellmodeller
	double length();
	double radius();
	double volume() ;
	btVector3 position();
	void setId(int _id) ;
	int getId() ;
};

/*%extend Solver {
	PyObject * divide(Cell * p, double f1=1.0, double f2=1.0)
	{
		Cell * d1=0, * d2=0;
		self->divideCell(p, f1, f2, d1, d2);
		PyObject * rval = PyList_New(2);
		PyObject *robj1, *robj2;
		robj1 =  SWIG_NewPointerObj(SWIG_as_voidptr(d1), SWIGTYPE_p_Cell, 0 |  0 );
		robj2 =  SWIG_NewPointerObj(SWIG_as_voidptr(d2), SWIGTYPE_p_Cell, 0 |  0 );
		PyList_SetItem(rval, 0, robj1);
		PyList_SetItem(rval, 1, robj2);
		return rval;
	}
};*/

%feature("shadow") Solver::ends(int id) %{
def ends(self,*args):
	e = $action(self,*args)
	v = ([e.x1, e.y1, e.z1], [e.x2, e.y2, e.z2])
	return v
%}

%extend Solver{
	CellEnds ends(int id)
	{
		btVector3 x1(-1,0,0); btVector3 x2(1,0,0);
		CellEnds ends;
                Cell *cell = self->cells[self->cellIdToIdx[id]];
		btTransform& t = cell->transform;
		double l = self->length(id);//+self->radius(id);
		btVector3 p1 = t*(0.5*l*x1); 
		btVector3 p2 = t*(0.5*l*x2);
		btScalar * p1val(p1); btScalar * p2val(p2);
		ends.x1 = p1val[0]; ends.y1 = p1val[1]; ends.z1 = p1val[2];
		ends.x2 = p2val[0]; ends.y2 = p2val[1]; ends.z2 = p2val[2];
		return ends;
	}
};

%extend Solver {
	PyObject * pos(int id)
	{
		btVector3 p = self->cells[self->cellIdToIdx[id]]->position();
		btScalar * pval(p);
		PyObject * rval = PyList_New(3);
		PyObject *robj;
		for (int i=0; i<3; ++i)
		{
			robj = SWIG_From_double(static_cast< double >(pval[i]));
			PyList_SetItem(rval, i, robj);
		}
		return rval;
	}
};

%extend Solver {
	void setGrowthRate(int id, double r)
	{
                Cell *c = self->cells[self->cellIdToIdx[id]];
		c->lenVel = r;
	}
};

class Solver
{

public: 
	Solver();
    void init();
	~Solver();

	void collide(double dt);
	void collideCapsules(double dt);
	double CGSolve(int maxIters, double dt);
	double integrate(double dt);

	// Interface for cellmodeller
	int step(double dt) ;
        void addCell(int id);
	int numCells() ;
	Cell * p_getCell(int i) ;

	%pythoncode %{
	def getCells(self):
		cells = []
		for i in range(self.numCells()):
			cells.append(self.p_getCell(i))
		return cells
	%}
	%pythoncode %{
	def getNeighbours(self):
		return None
	def setRegulator(self, reg):
		pass
	%}
        void divideCell(int pid, int d1id, int d2id, double f1, double f2);

	// cell operations
	double length(int cid); 
	double radius(int cid);
	double volume(int cid);
};
