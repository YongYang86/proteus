from proteus import Domain
from proteus.mprans import NCLS
from proteus import Norms
from proteus import Profiling
import numpy as np

#timeIntegration_ncls = "SSP33"
timeIntegration_ncls = "FE"
lRefinement=3
#end time of simulation
T=2.

fullNewton=False
# ENTROPY VISCOSITY and ART COMPRESSION PARAMETERS
EDGE_VISCOSITY=1
ENTROPY_VISCOSITY=0
POWER_SMOOTHNESS_INDICATOR=2
LUMPED_MASS_MATRIX=0
FCT=0
cK=0.0
# FOR EDGE BASED ENTROPY VISCOSITY
cE=1.0
cMax=0.1
# FOR SUPG 
shockCapturingFactor_ncls=0.2
#Other time parameters
if timeIntegration_ncls == "SSP33":
    timeOrder = 3
else:
    timeOrder = 1

runCFL = 0.1#0.3,0.185,0.125 for dgp1,dgp2,dgpk(3)
lag_shockCapturing_ncls=True
#if True uses PETSc solvers
parallel = False
linearSmoother = None
#compute mass balance statistics or not
checkMass=False
#number of space dimensions
nd=2
#time integration, not relevant if using BDF with cfl timestepping
rtol_u = {0:1.0e-4}
atol_u = {0:1.0e-4}
rtol_res = {0:1.0e-4}
atol_res = {0:1.0e-4}
#
#spatial approximation orders
cDegree_ncls=0
pDegree_ncls=1
useHex=False
useMetrics=0.0
#
#spatial quadrature orders
quad_order = 2*pDegree_ncls+1
#parallel partitioning info
from proteus import MeshTools
partitioningType = MeshTools.MeshParallelPartitioningTypes.node
#spatial mesh
#tag simulation name to level of refinement
nn=nnx=(2**lRefinement)*10+1
nny=(nnx-1)/10+1
nnz=1
he=1.0/(nnx-1.0)
L=[1.0,1.0]

unstructured=False #True for tetgen, false for tet or hex from rectangular grid
box=Domain.RectangularDomain(L=(1.0,0.1),
                             x=(0.0,0.0),
                             name="box");
box.writePoly("box")
domain = box
#number of output time steps
nDTout = 10
#smoothing factors
#eps
epsFactHeaviside=epsFactDirac=epsFact_ncls=1.5 #1.5
#
if useMetrics:
    shockCapturingFactor_ncls=0.5
    lag_shockCapturing_ncls=True

#use absolute tolerances on al models
atolRedistance = max(1.0e-12,0.1*he)
atolConservation = max(1.0e-12,0.001*he**2)
atolVolumeOfFluid= max(1.0e-12,0.001*he**2)
atolLevelSet     = max(1.0e-12,0.001*he**2)
#controls
linearSolverConvergenceTest = 'r-true' #rits is do a set number of iterations, r-true uses true residual, PETSc default is preconditioned residual
#redist solver
fmmFlag=0
#
if useHex:
    hex=True
    soname="oneD_advection_c0q"+`pDegree_ncls`+"_"+timeIntegration_ncls+"_level_"+`lRefinement`
else:
    soname="oneD_advection_c0p"+`pDegree_ncls`+"_"+timeIntegration_ncls+"_level_"+`lRefinement`

#My Own Coefficients
class MyCoefficients(NCLS.Coefficients):
    def attachModels(self,modelList):
        self.model = modelList[self.modelIndex]
        self.u_dof_old = np.copy(self.model.u[0].dof)
        self.u_dof_old_old = np.copy(self.model.u[0].dof)
        self.q_v = np.zeros((self.model.mesh.nElements_global,self.model.nQuadraturePoints_element,self.model.nSpace_global),'d')+1E10
        self.ebqe_v = np.zeros((self.model.mesh.nExteriorElementBoundaries_global,self.model.nElementBoundaryQuadraturePoints_elementBoundary,self.model.nSpace_global),'d')
        self.model.q[('velocity',0)]=self.q_v
        self.model.ebqe[('velocity',0)]=self.ebqe_v
        if self.RD_modelIndex != None:
            #print self.RD_modelIndex,len(modelList)
            self.rdModel = modelList[self.RD_modelIndex]
        else:
            self.rdModel = self.model

        # Divergence. Assume the velocity is div free
        self.q_div_velocity = np.zeros(self.model.q[('u', 0)].shape,'d')
        self.ebqe_div_velocity = np.zeros(self.model.ebqe[('u', 0)].shape,'d')
    def preStep(self,t,firstStep=False):
        # SAVE OLD SOLUTIONS
        self.u_dof_old_old = np.copy(self.u_dof_old)
        self.u_dof_old = np.copy(self.model.u[0].dof)

        # GET VELOCITY AT DOFs (FOR EDGE BASED METHODS)
        x_dof = self.model.u[0].femSpace.mesh.nodeArray[:,0]
        y_dof = self.model.u[0].femSpace.mesh.nodeArray[:,1]
        self.velx_tn_dof = 0.*y_dof+1
        self.vely_tn_dof = 0.*x_dof

        # GET VELOCITY AT QUADRATURE POINTS (FOR CELL BASE METHODS)
        x = self.model.q['x'][...,0]
        y = self.model.q['x'][...,1]
        x_boundary = self.model.ebqe['x'][...,0]
        y_boundary = self.model.ebqe['x'][...,1]

        self.q_v[...,0]  = 1.
        self.q_v[...,1]  = 0.

        self.ebqe_v[...,0]  = 1.
        self.ebqe_v[...,1]  = 0.
        
        #DIVERGENCE OF VEOCITY
        self.q_div_velocity = 0*x
        self.ebqe_div_velocity[...] = 0*x_boundary

        # CHECK MASS
        if self.checkMass:
            self.m_pre = Norms.scalarDomainIntegral(self.model.q['dV_last'],
                                                    self.model.q[('m',0)],
                                                    self.model.mesh.nElements_owned)
            Profiling.logEvent("Phase  0 mass before NCLS step = %12.5e" % (self.m_pre,),level=2)

        copyInstructions = {}
        return copyInstructions
    def postStep(self,t,firstStep=False):
        if(self.FCT==1 and False):
            self.model.FCTStep()
        if self.checkMass:
            self.m_post = Norms.scalarDomainIntegral(self.model.q['dV'],
                                                     self.model.q[('m',0)],
                                                     self.model.mesh.nElements_owned)
            Profiling.logEvent("Phase  0 mass after NCLS step = %12.5e" % (self.m_post,),level=2)
            #Compare mass before and after step
            #np.testing.assert_almost_equal(self.m_pre,self.m_post, err_msg="Mass before and after step are not almost equal", verbose=True)

        copyInstructions = {}
        return copyInstructions
    def evaluate(self,t,c):
        pass
