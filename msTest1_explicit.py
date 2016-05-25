""" Author: Hongyang Cheng <chyalexcheng@gmail>
    Test #1: 2D Membrane-wrapped granular material
"""
from esys.escript import *
from esys.weipa import saveVTK
from esys.finley import ReadGmsh
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
# name of mesh file
mshName = 'MshQuad3_0';
# Mesh with 8-node triangle elements; each element has 4 Gauss point
if mshName[3:7] == 'Quad': numOfElements = (int(mshName[-3])*2)**2
# Mesh with 3-node triangle elements; each element has 1 Gauss point
else:	numOfElements = 2*(int(mshName[-3])*2)**2
# number of Gauss points
gp = 4; numg = gp*numOfElements;
packNo = range(0,numg)
# density and damping ratio
rho = 2254.; damp = .2
# number of processes in multiprocessing
nump = 32
# safety factor for timestep size and real-time duration of simulation 
safe = 4.0; duration = 25
# directory for exterior DE scenes and variables
sceneExt ='./DE_exts/Test1/'
# import node IDs of membrane in exterior DE domain
mIds = numpy.load(sceneExt+'mNodesIds'+mshName+'.npy')
# import FE-DE mapping of boundary node IDs
FEDENodeMap = numpy.load(sceneExt+'FEDENodeMap'+mshName+'.npy').item()
# state filename of initial membrane DE elements
DE_ext = './DE_exts/Test1/DE_ext_'+mshName+'.yade.gz'
#~ # import FE-DE mapping of boundary element IDs (deprecated)
#~ FEDEBoundMap = numpy.load(sceneExt+'FEDEBoundMap'+mshName+'.npy').item()
# file to write force on the bottom
graphDir = './result/graphs/msTest1_Explicit/gp'+str(gp)+'/'
fout=file(graphDir+'safe_%1.1f_'%safe+'t_%1.1f_'%duration+mshName+'.dat','w')

###################
##  model setup  ##
###################

# multiscale model description
dom = ReadGmsh(mshName[:-2]+'.msh',numDim=dim,integrationOrder=2)
prob = MultiScale(domain=dom,dim=dim,ng=numg,np=nump,rho=rho,mIds=mIds,\
                  FEDENodeMap=FEDENodeMap,DE_ext=DE_ext)

# nodal coordinate
dom = prob.getDomain()
x = dom.getX()
bx = FunctionOnBoundary(dom).getX()

# Dirichlet BC positions: four corners fixed, bottom fixed along y
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

# compute appropriate timestep from eigenvalue
eigFreq = sqrt(prob.getMaxEigenvalue())
dt = safe*(2./eigFreq)

# compute appropriate timestep from PWave velocity (deprecated)
#~ T = prob.getCurrentTangent()
#~ maxM = max(T[0,0,0,0].toListOfTuples())
#~ PwaveVel = sqrt(maxM/rho)
#~ dt = safe*inf(prob.getDomain().getSize()/PwaveVel)

# initialize partial differential equation
prob.initialize(specified_u_mask=Dbc, specified_u_val=Dbc_val, dt=dt)

########################################
##  Run simulations for nt timesteps  ##
########################################

# start the simulation
time_start = time.time()
t = 1
nt = int(duration/dt)
# directory to export vtk data and packing scenes
Dir = 'msTest1/explicit/gp'+str(gp)+'/'+mshName+'/'
vtkDir = './result/vtk/'+Dir
packDir = './result/packing/'+Dir
gaussDir = './result/gauss/'+Dir
tWrite = int(0.25/dt)

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
      # get mask for boundary nodes on the bottom
      botSurf = whereZero(bx[1])
      # traction at the bottom
      tractBot = traction*botSurf
      # resultant force at the bottom
      forceBot = integrate(tractBot,where=FunctionOnBoundary(dom))
      # length of the bottom surface
      lengthBot = integrate(botSurf,where=FunctionOnBoundary(dom))
      # force magnitude
      magforceBot = sqrt(forceBot.dot(forceBot))
      # write stress at the bottom surface
      fout.write(str(t*dt)+' '+str(magforceBot)+' '+str(lengthBot)+'\n')
      
      # get local void ratio
      vR = prob.getLocalVoidRatio(); vR = proj(vR)
      # get local fabric intensity
      fab = prob.getLocalFabric()
      dev_fab = 4.*(fab-trace(fab)/dim*kronecker(prob.getDomain()))
      anis = sqrt(.5*inner(dev_fab,dev_fab))
      # set anis to -1 if no contact
      for i in range(numg):
         if math.isnan(anis.getTupleForDataPoint(i)[0]): anis.setValueOfDataPoint(i,-1)
      # get local rotation
      rot = prob.getLocalAvgRotation(); rot = proj(rot)
      # get local shear strain
      strain = prob.getCurrentStrain()
      volume_strain = trace(strain)
      dev_strain = symmetric(strain) - volume_strain*kronecker(prob.getDomain())/dim
      shear = sqrt(2*inner(dev_strain,dev_strain)); shear = proj(shear)
      
      # export FE scene
      saveVTK(vtkDir+"/ms"+mshName+"FE_%d.vtu"%t,u=u,sig=sig,shear=shear,e=vR,rot=rot,anis=anis)
      # export DE scenes
      prob.VTKExporter(vtkDir=vtkDir+"/ms"+mshName+"DE",t=t)
      # export local responses at Gauss points
      saveGauss2D(gaussDir+"/time_"+str(t)+".dat",strain=strain,stress=stress,fab=fab)
      print "force at the bottom: %e"%magforceBot
      
   # next iteration
   print "Step NO.%d finished, L2 norm of velocity at %2.1es: %e"%(t,t*dt,L2(u_t))
   t += 1

prob.getCurrentPacking(pos=(),time=t,prefix=packDir)
time_elapse = time.time() - time_start
fout.write("#Elapsed time in hours: "+str(time_elapse/3600.)+'\n')   
fout.close()
prob.exitSimulation()
