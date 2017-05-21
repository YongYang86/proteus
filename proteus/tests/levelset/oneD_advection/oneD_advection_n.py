from proteus import *
from proteus.default_n import *
from oneD_advection_p import *
from oneD_advection import *

multilevelNonlinearSolver  = NLNI
levelNonlinearSolver = Newton
fullNewtonFlag = fullNewton
updateJacobian = False

timeIntegration = NCLS.RKEV # SSP33 #mwf right now need timeIntegration to be SSP33 to run
stepController = Min_dt_controller#Min_dt_RKcontroller#Min_dt_controller #mwf we should probably    
if timeIntegration_ncls == "SSP33": #mwf hack
    timeOrder = 3
    nStagesTime = 3
else:
    timeOrder = 1
    nStagesTime = 1


if cDegree_ncls==0:
    if useHex:
        if pDegree_ncls==1:
            femSpaces = {0:C0_AffineLinearOnCubeWithNodalBasis}
        elif pDegree_ncls==2:
            femSpaces = {0:C0_AffineLagrangeOnCubeWithNodalBasis}
    else:
        if pDegree_ncls==1:
            femSpaces = {0:C0_AffineLinearOnSimplexWithNodalBasis}
        elif pDegree_ncls==2:
            femSpaces = {0:C0_AffineQuadraticOnSimplexWithNodalBasis}
    subgridError = None #Advection_ASGS(coefficients,nd,lag=False)
    shockCapturing = NCLS.ShockCapturing(coefficients,nd,shockCapturingFactor=shockCapturingFactor_ncls,lag=lag_shockCapturing_ncls)
    if parallel or LevelModelType == NCLS.LevelModel:
        numericalFluxType = Advection_DiagonalUpwind_IIPG_exterior

elif cDegree_ncls==-1:
    if pDegree_ncls==0:
        femSpaces = {0:DG_AffineP0_OnSimplexWithMonomialBasis}
    elif pDegree_ncls==1:
        femSpaces = {0:DG_AffineP1_OnSimplexWithMonomialBasis}
    elif pDegree_ncls==2:
        femSpaces = {0:DG_AffineP2_OnSimplexWithMonomialBasis}
    elif pDegree_ncls==3:
        femSpaces = {0:DG_AffineP3_OnSimplexWithMonomialBasis}
    numericalFluxType = Advection_DiagonalUpwind
    limiterType =   {0:TimeIntegration.DGlimiterPkMonomial2d}

if useHex:
    elementQuadrature = CubeGaussQuadrature(nd,quad_order)
    elementBoundaryQuadrature = CubeGaussQuadrature(nd-1,quad_order)
else:
    elementQuadrature = SimplexGaussQuadrature(nd,quad_order)
    elementBoundaryQuadrature = SimplexGaussQuadrature(nd-1,quad_order)

nonlinearSmoother = None#NLGaussSeidel

tolFac = 0.0
linTolFac = tolFac

nl_atol_res = 100*atolVolumeOfFluid
l_atol_res = atolVolumeOfFluid

maxNonlinearIts = 20
maxLineSearches = 0

matrix = SparseMatrix

if parallel:
    multilevelLinearSolver = KSP_petsc4py#PETSc
    levelLinearSolver = KSP_petsc4py#PETSc
    linear_solver_options_prefix = 'ncls_'
    linearSolverConvergenceTest = 'r-true'
else:
    multilevelLinearSolver = LU
    levelLinearSolver = LU

conservativeFlux = {}
if checkMass:
    auxiliaryVariables = [MassOverRegion()]
