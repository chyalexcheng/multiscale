""" Author: Hongyang Cheng <chyalexcheng@gmail>
    Test #1: 2D Membrane-wrapped granular material
"""
from esys.escript import *
from esys.weipa import saveVTK
from esys.escript.pdetools import Projector
from esys.escript.linearPDEs import LinearPDE,SolverOptions

from msFEM2DExplicit import MultiScale
from saveGauss import saveGauss2D
import time
import shutil

####################
##  key controls  ##
####################

# sample size, 1.2m by 1.2m
dim = 2; lx = 1.2; ly = 1.2
# read Gmsh mesh with 6-node triangle element (8 tri6); each element has 4 Gauss points
mshName = 'Msh2'
# number of Gauss points
numg = 3*32
packNo=range(0,numg)
# density and damping ratio
rho = 2254.; damp = .2
# number of processes in multiprocessing
nump = 32
# safety factor for timestep size and real-time duration of simulation 
safe = 2.0; duration = 25
# import membrane node Ids in exterior DE domain
mIds = numpy.load('mNodesIds'+mshName+'.npy')
# import FE-DE boundary node mapping
FEDENodeMap = numpy.load('FEDENodeMap'+mshName+'.npy').item()
# import FE-DE boundary element mapping
FEDEBoundMap = numpy.load('FEDEBoundMap'+mshName+'.npy').item()
# boundary pressure on DE membrane
conf = 0
# open file to write force on the top surface with its length
graphDir = './result/graphs/msTest1_Explicit/'
fout=file(graphDir+'safe_%1.1f_'%safe+'t_%1.1f_'%duration+mshName+'.dat','w')

###################
##  model setup  ##
###################

# multiscale model description
prob = MultiScale(mshName=mshName,dim=dim,ng=numg,np=nump,rho=rho,mIds=mIds,\
                  FEDENodeMap=FEDENodeMap,FEDEBoundMap=FEDEBoundMap,conf=conf)

# nodal coordinate
dom = prob.getDomain()
x = dom.getX()
bx = FunctionOnBoundary(dom).getX()

# Dirichlet BC positions, fixed four corners
Dbc = whereZero(x[0])*whereZero(x[1])*[1,1] +\
      whereZero(x[0]-lx)*whereZero(x[1])*[1,1] +\
      whereZero(x[1])*[0,1]

# Dirichlet BC values
Dbc_val = whereZero(x[0])*whereZero(x[1])*[0,0] +\
          whereZero(x[0]-lx)*whereZero(x[1])*[0,0] +\
          whereZero(x[1])*[0,0]

######################
##  Initialization  ##
######################

# compute appropriate size of timestep
eigFreq = sqrt(prob.getMaxEigenvalue())
dt = safe*(2./eigFreq)

#~ T = prob.getCurrentTangent()
#~ maxM = max(T[0,0,0,0].toListOfTuples())
#~ PwaveVel = sqrt(maxM/rho)
#~ dt = safe*inf(mydomain.getSize()/PwaveVel)

# initialize partial difference equation and return timestep
prob.initialize(specified_u_mask=Dbc, specified_u_val=Dbc_val, dt=dt)

########################################
##  Run simulations for nt timesteps  ##
########################################

# start the simulation
time_start = time.time()
t = 1
nt = int(duration/dt)
# directory to save vtk data
vtkDir = "./result/vtk/msTest1/explicit"   

while t <= nt:
   # update displacement and velocity at (n+1) timesteps
   u, u_t = prob.solve(damp=damp)
   # write data at selected timesteps
   if t%int(500/safe) == 0:
      # get stress at (n) timesteps
      stress = prob.getCurrentStress()
      dom = prob.getDomain()
      proj = Projector(dom)
      # project Gauss point value to nodal value
      sig = proj(stress)
      # interpolate to stress at the boundary
      sig_bounda = interpolate(sig,FunctionOnBoundary(dom)) 
      # compute boundary traction by s_ij*n_j
      traction = matrix_mult(sig_bounda,dom.getNormal())
      # get mask for boundary nodes on bottom surface
      bottomSurf = whereZero(bx[1])
      # traction at bottom surface
      tractBottom = traction*bottomSurf
      # resultant force at bottom
      forceBottom = integrate(tractBottom,where=FunctionOnBoundary(dom))
      # length of bottom surface
      lengthBottom = integrate(bottomSurf,where=FunctionOnBoundary(dom))
      # force magnitude
      magForceBottom = sqrt(forceBottom.dot(forceBottom))
      # write stress on the bottom
      fout.write(str(t*dt)+' '+str(magForceBottom)+' '+str(lengthBottom)+'\n')      
      # export FE scene
      saveVTK(vtkDir+"/ms"+mshName+"FE_%d.vtu"%t,u=u,u_t=u_t,sig=sig)
      # export DE scene
      prob.VTKExporter(vtkDir=vtkDir+"/ms"+mshName+"DE",t=t)
      print "force at bottom: %e"%magForceBottom
   # next iteration
   print "Step NO.%d finished, L2 norm of velocity at %2.1es: %e"%(t,t*dt,L2(u_t))
   t += 1

#~ prob.getCurrentPacking(pos=(),time=t,prefix='./result/packing/test1/')
time_elapse = time.time() - time_start
fout.write("#Elapsed time in hours: "+str(time_elapse/3600.)+'\n')   
fout.close()
prob.exitSimulation()
