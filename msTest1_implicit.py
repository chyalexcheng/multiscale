""" Author: Hongyang Cheng <chyalexcheng@gmail>
    Test #1: 2D Membrane-wrapped granular material
"""
from esys.escript import *
from esys.finley import ReadGmsh
from esys.weipa import saveVTK
from esys.escript.pdetools import Projector
from esys.escript.linearPDEs import LinearPDE,SolverOptions

from msFEM2DImplicit import MultiScale
from saveGauss import saveGauss2D
import time
import shutil

####################
##  key controls  ##
####################

# sample size, 1.2m by 1.2m
dim =2; lx = 1.2; ly = 1.2
# read Gmsh mesh with 6-node triangle element (8 tri6); each element has 3 Gauss points
mshName = 'Test3'
mydomain = ReadGmsh(mshName+'.msh',numDim=dim,integrationOrder=2)
# number of Gauss points
numg = 3*64
packNo=range(0,numg)
# density and damping ratio
rho = 2254.; damp = .2
# number of processes in multiprocessing
nump = 16
# max iteration steps per timestep, safety factor for dt and total number of timesteps 
iterMax = 100; safetyFactor = 0.1; duration = 1.
# import FE-DE boundary node mapping
FEDEmap = numpy.load('FEDEmap'+mshName+'.npy').item()
# import membrane node Ids in exterior DE domain
mIds = numpy.load('mNodesIds'+mshName+'.npy')
# boundary pressure on DE membrane
confining = -0.1*1e4
# open file to write force on the top surface with its length
graphDir = './result/graphs/ms'+mshName+'_Implicit/Iter'
fout=file(graphDir+str(iterMax)+'_safe'+str(safetyFactor)+'.dat','w')

###################
##  model setup  ##
###################

# multiscale model description
prob = MultiScale(domain=mydomain,ng=numg,np=nump,rho=rho,random=False,rtol=1e-2,\
                  verbose=True,mIds=mIds,FEDEmap=FEDEmap,confining=confining)

# nodal coordinate
x = mydomain.getX()
bx = FunctionOnBoundary(mydomain).getX()

# Dirichlet BC positions, fixed four corners
Dbc = whereZero(x[0])*whereZero(x[1])*[1,1] +\
      whereZero(x[0]-lx)*whereZero(x[1])*[1,1] +\
      whereZero(x[0]-lx)*whereZero(x[1]-ly)*[1,1] +\
      whereZero(x[0])*whereZero(x[1]-ly)*[1,1]
# Dirichlet BC values
Dbc_val = whereZero(x[0])*whereZero(x[1])*[0,0] +\
          whereZero(x[0]-lx)*whereZero(x[1])*[0,0] +\
          whereZero(x[0]-lx)*whereZero(x[1]-ly)*[0,0] +\
          whereZero(x[0])*whereZero(x[1]-ly)*[0,0]

# initial stress
stress = prob.getCurrentStress()
proj = Projector(mydomain)
# project Gauss point value to nodal value
sig = proj(stress)
# interpolate to stress at the boundary
sig_bounda = interpolate(sig,FunctionOnBoundary(mydomain)) 
# compute boundary traction by s_ij*n_j
traction = matrix_mult(sig_bounda,mydomain.getNormal())
# get mask for boundary nodes on top surface
topSurf = whereNonNegative(bx[1]-ly)
# traction at top surface
tractTop = traction*topSurf
# resultant force at top
forceTop = integrate(tractTop,where=FunctionOnBoundary(mydomain))
# length of top surface
lengthTop = integrate(topSurf,where=FunctionOnBoundary(mydomain))
# force magnitude
magForceTop = sqrt(forceTop.dot(forceTop))
# write stress on the top
fout.write('0 '+str(magForceTop/lengthTop)+'\n')

######################
##  Initialization  ##
######################

# compute appropriate size of timestep
T = prob.getCurrentTangent()
maxM = max(T[0,0,0,0].toListOfTuples())
PwaveVel = sqrt(maxM/rho)
dt = safetyFactor*inf(mydomain.getSize()/PwaveVel)

# initialize partial difference equation
prob.initialize(specified_u_mask=Dbc, specified_u_val=Dbc_val, dt=dt)

########################################
##  Run simulations for nt timesteps  ##
########################################

# start the simulation
time_start = time.time()
t = 1
nt = int(duration/(safetyFactor*dt))
# directory to save vtk data
vtkDir = "./result/vtk/ms"+mshName
shutil.rmtree(vtkDir)
mkDir(vtkDir)

while t <= nt:
   # update displacement and velocity at (n+1) timesteps
   u, u_t = prob.solve(iter_max=iterMax)
   # write data at selected timesteps
   if t%int(1./safetyFactor) == 0:
      # get stress at (n) timesteps
      stress = prob.getCurrentStress()
      proj = Projector(mydomain)
      # project Gauss point value to nodal value
      sig = proj(stress)
      # interpolate to stress at the boundary
      sig_bounda = interpolate(sig,FunctionOnBoundary(mydomain)) 
      # compute boundary traction by s_ij*n_j
      traction = matrix_mult(sig_bounda,mydomain.getNormal())
      # get mask for boundary nodes on top surface
      topSurf = whereNonNegative(bx[1]-ly)
      # traction at top surface
      tractTop = traction*topSurf
      # resultant force at top
      forceTop = integrate(tractTop,where=FunctionOnBoundary(mydomain))
      # length of top surface
      lengthTop = integrate(topSurf,where=FunctionOnBoundary(mydomain))
      # force magnitude
      magForceTop = sqrt(forceTop.dot(forceTop))
      # write stress on the top
      fout.write(str(t*dt)+' '+str(magForceTop/lengthTop)+'\n')      
      # export FE scene
      saveVTK(vtkDir+"/ms"+mshName+"FE_%d.vtu"%t,u=u,u_t=u_t,sig=sig)
      # export DE scene
      prob.VTKExporter(vtkDir=vtkDir+"/ms"+mshName+"DE",t=t)
   # next iteration
   print "Step NO.%d finished, L2 norm of velocity: %e"%(t,L2(u_t))
   t += 1

time_elapse = time.time() - time_start
fout.write("#Elapsed time in hours: "+str(time_elapse/3600.)+'\n')   
fout.close()
prob.exitSimulation()
