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

####################
##  key controls  ##
####################

# sample size, 1.2m by 1.2m
dim = 2; lx = 1.2; ly = 1.2
# read Gmsh mesh with 6-node triangle element (8 tri6); each element has 4 Gauss points
mshName = 'Msh2'; numOfElements = 32
# number of Gauss points
gp = 1; numg = gp*numOfElements;
packNo=range(0,numg)
# density and damping ratio
rho = 2254.; damp = .2
# number of processes in multiprocessing
nump = 32
# safety factor for timestep size and real-time duration of simulation 
safe = 2.0; duration = 25
# directory for exterior DE scenes and variables
sceneExt ='./DE_exts/Test2/'
# import membrane node Ids in exterior DE domain
mIds = numpy.load(sceneExt+'mNodesIds'+mshName+'.npy')
# import FE-DE boundary node mapping
FEDENodeMap = numpy.load(sceneExt+'FEDENodeMap'+mshName+'.npy').item()
# surcharnge pressure
surcharge=-2.e4
# open file to write force on the top surface with its length
graphDir = './result/graphs/msTest2_Explicit/gp'+str(gp)+'/'
fout=file(graphDir+'safe_%1.1f_'%safe+'t_%1.1f_'%duration+mshName+'.dat','w')

###################
##  model setup  ##
###################

# multiscale model description
prob = MultiScale(mshName=mshName,dim=dim,ng=numg,np=nump,rho=rho,\
						mIds=mIds,FEDENodeMap=FEDENodeMap)

# nodal coordinate
dom = prob.getDomain()
x = dom.getX()
bx = FunctionOnBoundary(dom).getX()

# Dirichlet BC positions, smooth lateral boundary
Dbc = whereZero(x[0])*[1,0] + whereZero(x[0]-lx)*[1,0]

# Dirichlet BC values
Dbc_val = whereZero(x[0])*[0,0] + whereZero(x[0]-lx)*[0,0]

# Neumann BC postions, known pressure on the top
Nbc = whereZero(bx[1]-ly)*[0,surcharge]

######################
##  Initialization  ##
######################

# compute appropriate size of timestep
eigFreq = sqrt(prob.getMaxEigenvalue())
dt = safe*(2./eigFreq)

#~ T = prob.getCurrentTangent()
#~ maxM = max(T[0,0,0,0].toListOfTuples())
#~ PwaveVel = sqrt(maxM/rho)
#~ dt = safe*inf(prob.getDomain().getSize()/PwaveVel)

# initialize partial difference equation and return timestep
prob.initialize(f=Nbc, specified_u_mask=Dbc, specified_u_val=Dbc_val, dt=dt)

########################################
##  Run simulations for nt timesteps  ##
########################################

# start the simulation
time_start = time.time()
t = 1
nt = int(duration/dt)
# directory to export vtk data and packing scenes
Dir = 'msTest2/explicit/gp'+str(gp)+'/'+mshName+'/'
vtkDir = './result/vtk/'+Dir
packDir = './result/packing/'+Dir

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
      topSurf = whereZero(bx[1]-ly)
      # traction at bottom surface
      tractTop = traction*topSurf
      # resultant force at bottom
      forceTop = integrate(tractTop,where=FunctionOnBoundary(dom))
      # length of bottom surface
      lengthTop = integrate(topSurf,where=FunctionOnBoundary(dom))
      # write stress on the bottom
      fout.write(str(t*dt)+' '+str(forceTop[0])+' '+str(forceTop[1])+' '+str(lengthTop)+'\n')      
      # get local void ratio
      vR=prob.getLocalVoidRatio()
      #~ # get local fabric intensity
      #~ fabric=prob.getLocalFabric()
      #~ iso_fabric = trace(fabric)
      #~ aniso_fabric = symmetric(fabric) - iso_fabric*kronecker(prob.getDomain())/dim
      #~ aniso = sqrt(1./2*inner(aniso_fabric,aniso_fabric))
      # get local rotation
      rotation=prob.getLocalAvgRotation()
      # get local shear strain
      strain = prob.getCurrentStrain()
      volume_strain = trace(strain)
      dev_strain = symmetric(strain) - volume_strain*kronecker(prob.getDomain())/dim
      shear = sqrt(2*inner(dev_strain,dev_strain))
      # export FE scene
      #~ saveVTK(vtkDir+"/ms"+mshName+"FE_%d.vtu"%t,u=u,sig=sig)
      saveVTK(vtkDir+"/ms"+mshName+"FE_%d.vtu"%t,u=u,sig=sig,shear=shear,e=vR,rot=rotation)
      # export DE scene
      prob.VTKExporter(vtkDir=vtkDir+"/ms"+mshName+"DE",t=t)
      print "stress ratio at bottom: %e"%(forceTop[0]/forceTop[1])
   # next iteration
   print "Step NO.%d finished, L2 norm of velocity at %2.1es: %e"%(t,t*dt,L2(u_t))
   t += 1

prob.getCurrentPacking(pos=(),time=t,prefix=packDir)
time_elapse = time.time() - time_start
fout.write("#Elapsed time in hours: "+str(time_elapse/3600.)+'\n')   
fout.close()
prob.exitSimulation()
