from proteus import *
from proteus.default_p import *
from proteus.ctransportCoefficients import smoothedHeaviside
from math import *
from vortex2D import *

LevelModelType = NCLS.LevelModel
logEvent = Profiling.logEvent
name=soname+"_ls"

nd=2

class init_cond:
    def __init__(self,L):
        self.radius = 0.15
        self.xc=0.5
        self.yc=0.75
    def uOfXT(self,x,t):
        #return self.radius - math.sqrt((x[0]-self.xc)**2 + (x[1]-self.yc)**2)
        beta = epsCoupez*he
        scaling = 1.
        dist = self.radius - math.sqrt((x[0]-self.xc)**2 + (x[1]-self.yc)**2)
        return dist
        # return scaling*beta*math.tanh(dist/beta)
        #return smoothedHeaviside(epsFactHeaviside*he,self.radius - math.sqrt((x[0]-self.xc)**2 + (x[1]-self.yc)**2))
        #return math.tanh((x[0]-self.xc))

class exact:
    def __init__(self,L):
        self.radius = 0.15
        self.xc=0.5
        self.yc=0.75
    def uOfXT(self,x,t):
        #return self.radius - math.sqrt((x[0]-self.xc)**2 + (x[1]-self.yc)**2)
        beta = epsCoupez*he
        dist = self.radius - math.sqrt((x[0]-self.xc)**2 + (x[1]-self.yc)**2)
        return dist
        return beta*math.tanh(dist/beta)
        #return smoothedHeaviside(epsFactHeaviside*he,self.radius - math.sqrt((x[0]-self.xc)**2 + (x[1]-self.yc)**2))
        #return math.tanh((x[0]-self.xc))

analyticalSolution = {0:exact(L)}

RD_model=None
coefficients = MyCoefficients(epsFact=epsFactHeaviside,checkMass=checkMass,RD_model=RD_model,useMetrics=useMetrics,
                              EDGE_VISCOSITY=EDGE_VISCOSITY, 
                              ENTROPY_VISCOSITY=ENTROPY_VISCOSITY,
                              LUMPED_MASS_MATRIX=LUMPED_MASS_MATRIX, 
                              lambda_coupez=lambda_coupez,
                              pure_redistancing=pure_redistancing,
                              redistancing_tolerance=redist_tolerance*he,
                              epsCoupez=epsCoupez*he,
                              epsFactRedistancing=epsFactRedistancing*he)
                              
coefficients.variableNames=['u']
initialConditions  = {0:init_cond(L)}

#now define the Dirichlet boundary conditions
def getDBC(x,flag):
    #return lambda x,t: 1.0
    pass
 
def zeroInflow(x):
    return lambda x,t: 0.0

dirichletConditions = {0:getDBC}
fluxBoundaryConditions = {0:'outFlow'}

def zeroadv(x,flag):
    return lambda x,t: 0.0

#advectiveFluxBoundaryConditions =  {}
advectiveFluxBoundaryConditions =  {0:zeroadv}
diffusiveFluxBoundaryConditions = {0:{}}
