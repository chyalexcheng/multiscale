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
# read Gmsh mesh with 3-node triangle element; each element has 1 Gauss point
mshName = 'Msh5'; numOfElements = 2*(int(mshName[3])*2)**2
# number of Gauss points
gp = 1; numg = gp*numOfElements;
packNo=range(0,numg)
# density and damping ratio
rho = 2254.; damp = .2
# number of processes in multiprocessing
nump = 32
# safety factor for timestep size and real-time duration of simulation 
safe = 4.0; duration = 15
# directory for exterior DE scenes and variables
sceneExt ='./DE_exts/Test2/'
# import membrane node Ids in exterior DE domain
mIds = numpy.load(sceneExt+'mNodesIds'+mshName+'.npy')
# import FE-DE boundary node mapping
FEDENodeMap = numpy.load(sceneExt+'FEDENodeMap'+mshName+'.npy').item()
# surcharnge pressure
surcharge=-2.e4
# open file to write force on the bottom surface with its length
graphDir = './result/graphs/msTest2_Explicit/gp'+str(gp)+'/'
fout=file(graphDir+'safe_%1.1f_'%safe+'t_%1.1f_'%duration+mshName+'.dat','w')

###################
##  model setup  ##
###################

# multiscale model description
prob = MultiScale(mshName=mshName[:4],dim=dim,ng=numg,np=nump,rho=rho,\
                  mIds=mIds,FEDENodeMap=FEDENodeMap)

# nodal coordinate
dom = prob.getDomain()
x = dom.getX()
bx = FunctionOnBoundary(dom).getX()

# Dirichlet BC positions, smooth lateral boundary
Dbc = whereZero(x[0])*[1,0] + whereZero(x[0]-lx)*[1,0] + \
      whereZero(x[0]-lx)*whereZero(x[1])*[1,1]

# Dirichlet BC values
Dbc_val = whereZero(x[0])*[0,0] + whereZero(x[0]-lx)*[0,0] + \
          whereZero(x[0]-lx)*whereZero(x[1])*[0,0]

# Neumann BC postions, known pressure on the top
Nbc = whereZero(bx[1]-ly)*[0,surcharge]

######################
##  Initialization  ##
######################

# compute appropriate size of timestep
eigFreq = sqrt(prob.getMaxEigenvalue())
dt = safe*(2./eigFreq)

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
gaussDir = './result/gauss/'+Dir
tWrite = int(0.15/dt)

while t <= nt:
   # update displacement and velocity at (n+1) timesteps
   u, u_t = prob.solve(damp=damp)
   # write data at selected timesteps
   if t%tWrite == 0:
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
      botSurf = whereZero(bx[1])
      # traction at bottom surface
      tractBot = traction*botSurf
      # resultant force at bottom
      forceBot = integrate(tractBot,where=FunctionOnBoundary(dom))
      # length of bottom surface
      lengthBot = integrate(botSurf,where=FunctionOnBoundary(dom))
      # write stress on the bottom
      fout.write(str(t*dt)+' '+str(forceBot[0])+' '+str(forceBot[1])+' '+str(lengthBot)+'\n')
      
      # get local void ratio
      vR = prob.getLocalVoidRatio(); vR = proj(vR)
      #~ # get local fabric intensity
      #~ fab = prob.getLocalFabric()
      #~ dev_fab = 4.*(fab-trace(fab)/dim*kronecker(prob.getDomain()))
      #~ anis = sqrt(.5*inner(dev_fab,dev_fab)); anis = proj(anis)
      # get local rotation
      rot = prob.getLocalAvgRotation(); rot = proj(rot)
      # get local shear strain
      strain = prob.getCurrentStrain()
      volume_strain = trace(strain)
      dev_strain = symmetric(strain) - volume_strain*kronecker(prob.getDomain())/dim
      shear = sqrt(2*inner(dev_strain,dev_strain)); shear = proj(shear)

      # export FE scene
      saveVTK(vtkDir+"/ms"+mshName+"FE_%d.vtu"%t,u=u,sig=sig,shear=shear,e=vR,rot=rot)
      #~ saveVTK(vtkDir+"/ms"+mshName+"FE_%d.vtu"%t,u=u,sig=sig,shear=shear,e=vR,rot=rot,anis=anis)
      # export DE scene
      prob.VTKExporter(vtkDir=vtkDir+"/ms"+mshName+"DE",t=t)
      # export local response at Gauss point
      #~ saveGauss2D(gaussDir+"/time_"+str(t)+".dat",strain=strain,stress=stress,fab=fab)
      print "stress ratio at bottom: %e"%(forceBot[0]/forceBot[1])

   # next iteration
   print "Step NO.%d finished, L2 norm of velocity at %2.1es: %e"%(t,t*dt,L2(u_t))
   t += 1

prob.getCurrentPacking(pos=(),time=t,prefix=packDir)
time_elapse = time.time() - time_start
fout.write("#Elapsed time in hours: "+str(time_elapse/3600.)+'\n')   
fout.close()
prob.exitSimulation()
