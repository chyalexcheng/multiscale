__author__="Ning Guo, ceguo@connect.ust.hk; Hongyang Cheng, d132535@hiroshima-u.ac.jp"
__supervisor__="Jidong Zhao, jzhao@ust.hk"
__institution__="The Hong Kong University of Science and Technology"

"""
DEM part for Hierachical Multiscale simulation:
   sets up packings representing material points (RVE) at Gauss points
   of the FE domain and returns updated stress
DEM part for Concurrent Multiscale simulation:
   sets up the concurrent DE domain which receives boundary displacement
   from the FE domain and returns back boundary nodal force to FE solver
"""

import time

# import path where yadeimport.py locates
import sys
from os.path import expanduser
# yadeimport.py is generated by `ln yade-versionNo yadeimport.py`
sys.path.append(expanduser('~')+'/myYade/yadeDevVer2/install/bin') 

# import yade modules
from yadeimport import *
from yade import export

### Below is added for surface coupling by Hongyang Cheng ###

# state filename of intial DE packing
DE_int = './DE_exts/Test2/DE_alum.yade.gz'

# load exterior DE domain scene 
def initLoadExt(DE_ext):
   Omega().load(DE_ext)
   return Omega().sceneToString()

# load exterior DE domain scene and return boundary force
def initNbc(scene,conf,mIds,FEDEMap):
   Omega().stringToScene(scene)
   if conf:
      # apply boundary force, if any, on DE membrane nodes
      applyBoundaryForce(conf,mIds)
   # get boundary traction on membrane elements
   FEf = getBoundaryForce(FEDEMap)
   return FEf
  
"""
# load exterior DE domain scene and return boundary traction (deprecated)
def initNbc(scene,conf,mIds,FEDEMap):
   Omega().stringToScene(scene)
   if conf:
      # apply boundary force, if any, on DE membrane nodes
      applyBoundaryForce(conf,mIds)
   # get boundary traction on membrane elements
   FENbc = getBoundaryTraction(FEDEMap)
   return FENbc
"""

# apply displacement on interface DE nodes and return nodal force
def moveInterface_getForce2D(scene,conf,DEdu,dynRelax,mIds,FEDEMap,ns):
   Omega().stringToScene(scene)
   if conf:
      # apply boundary force on membrane nodes
      applyBoundaryForce(conf,mIds)
   # if in dynamic relaxation mode
   if dynRelax:
      # increase damping ratio
      damp = newton.damping
      newton.damping = 0.99
      # stop loading 
      for n in [mIds[0],mIds[-1]]:
         pullSpeed = O.bodies[n].state.vel
         O.bodies[n].shape.color = pullSpeed
         O.bodies[n].state.vel = Vector3(0,0,0)
   # assign interface velocity
   b = Omega().bodies
   for key in DEdu.keys():
      du = utils.Vector3(DEdu[key][0],DEdu[key][1],0)
      v = du/(ns*Omega().dt)
      b[key].state.vel = utils.Vector3(v[0],v[1],0)   
   # run until reach target du
   Omega().run(ns,True)
   # get boundary force on membrane nodes
   FEf = getBoundaryForce(FEDEMap)
   # set loading speed and damping to normal after microscale simulations
   if dynRelax:
      # reset damping ratio
      newton.damping = damp
      # start loading
      for n in [mIds[0],mIds[-1]]:
         pullSpeed = O.bodies[n].shape.color
         O.bodies[n].state.vel = pullSpeed
   # save Omega to scene and reset Omega to reduce leaked memory
   scene = Omega().sceneToString()
   Omega().reset()
   return FEf, scene

"""
# apply displacement on interface DE nodes and return nodal traction
def moveInterface_getTraction2D(scene,conf,DEdu,mIds,FEDEMap,ns,vel=1.e-3,stabLimit=1.e-3):
   Omega().stringToScene(scene)
   if conf:
      # apply boundary force on membrane nodes
      applyBoundaryForce(conf,mIds)
   b = Omega().bodies
   #~ ns = int(max(numpy.max(numpy.abs(DEdu.values()))/vel/Omega().dt,2))
   for key in DEdu.keys():
      d = utils.Vector3(DEdu[key][0],DEdu[key][1],0)
      v = d/(ns*Omega().dt)
      b[key].state.vel = utils.Vector3(v[0],v[1],0)   
   # run until reach target du
   Omega().run(ns,True)
   # calm the system and run the simulation until quasi-static state
   #~ utils.calm()
   #~ for key in DEdu.keys(): b[key].state.vel = utils.Vector3.Zero
   #~ while True:
      #~ Omega().run(100,True); unb = utils.unbalancedForce()
      #~ if unb < stabLimit: break
   # get boundary traction on membrane elements
   FENbc = getBoundaryTraction(FEDEMap)
   # save Omega to scene and reset Omega to reduce leaked memory
   scene = Omega().sceneToString()
   Omega().reset()
   return FENbc, scene
"""

# get boundary force on interface nodes
def getBoundaryForce(FEDENodeMap,width=0.1):
   Omega().run(1,True)
   FEf = dict((k,0) for k in FEDENodeMap.keys())
   for key in FEDENodeMap.keys():
      DEid = FEDENodeMap[key]
      f = Omega().forces.f(DEid).xy()
      FEf[key] = (f[0]/width,f[1]/width)
   #~ print max(FEf.values())
   return FEf
 
# export interior DE scene in vtk format
def exportInt(params):
   scene = params[0]; vtkDir = params[1]; t = params[2]
   Omega().stringToScene(scene)
   pos = Omega().tags['id']
   sIntrs = export.VTKExporter(vtkDir+'int_'+pos)
   sIntrs.exportInteractions(
   numLabel=t,
      what=[('f_n','i.phys.normalForce.norm()'),
            ('f_s','i.phys.shearForce.norm()'),
           ])
  
# export exterior DE scene in vtk format
def exportExt(scene,mIds,vtkDir,t):
   Omega().stringToScene(scene)
   mIntrs=export.VTKExporter(vtkDir+'ext_')
   mIntrs.exportInteractions(
   numLabel=t,
      ids=[(i,j) for i,j in zip(mIds[:-1],mIds[1:])],
      what=[('Tlimit','i.phys.normalForce.norm()/i.phys.normalAdhesion'),
           ])

# apply deformation on 2D DEM packing
def shear2D(param):
   Omega().stringToScene(param[0])
   dstrain = utils.Matrix3(param[1][0],param[1][1],0, param[1][2],param[1][3],0, 0,0,0)
   ns = param[2]
   dynRelax = param[3]
   Omega().cell.velGrad=dstrain/(ns*Omega().dt)
   if dynRelax: 
      damp = Omega().engines[-1].damping
      Omega().engines[-1].damping = 0.99
      Omega().run(ns,True)
      Omega().engines[-1].damping = damp
   else:
      Omega().run(ns,True)
   return Omega().sceneToString()

# return time step of a DE scene
def getScenetDt(scene):
   Omega().stringToScene(scene)
   return Omega().dt

# apply hydrostatic pressure on membrane nodes
def applyBoundaryForce(conf,mIds,width=0.1):
   # run one step only to set permanent force
   Omega().run(1,True)
   for i in mIds:
      intrs = Omega().bodies[i].intrs()
      segts = []; n = utils.Vector3.Zero
      for intr in intrs:
         if isinstance(intr.phys,CohFrictPhys):
            segt = (intr.geom.refR1+intr.geom.refR2-intr.geom.penetrationDepth)*intr.geom.normal
            segts.append(segt); n += segt
      # segment normal
      n = n.normalized()
      dl = sum([.5*segt.dot(n) for segt in segts])
      # cross product of unit z and interaction normal points inward
      n = n.cross(utils.Vector3.UnitZ)
      # get boundary force per node: f = sigma*segt*width (negative sign for compression)
      permDEf = conf*(dl*width)*n
      # add permanent force on membrane nodes
      Omega().forces.addF(i,permDEf,True)


"""
# This is a function to test memory leak
def moveInterface_getForce2D(scene):
   Omega().stringToScene(scene)
   del scene
   Omega().exitNoBacktrace()
   FENbc = 0
   return FENbc, scene[:]

# get boundary traction on interface nodes (deprecated)
def getBoundaryTraction(FEDENodeMap,width=0.1):
   Omega().run(1,True)
   FENbc = dict((k,0) for k in FEDENodeMap.keys())
   for key in FEDENodeMap.keys():
      DEid = FEDENodeMap[key]
      f = Omega().forces.f(DEid).xy()
      intrs = Omega().bodies[DEid].intrs()
      segts = []; n = utils.Vector3.Zero
      for intr in intrs:
         if isinstance(intr.phys,CohFrictPhys):
            segt = (intr.geom.refR1+intr.geom.refR2-intr.geom.penetrationDepth)*intr.geom.normal
            segts.append(segt); n += segt
      # segment normal
      n = n.normalized()
      dl = sum([.5*segt.dot(n) for segt in segts])
      FENbc[key] = (f[0]/(dl*width),f[1]/(dl*width))
   #~ print max(FENbc.values())
   return FENbc


# get boundary traction on interface elements (deprecated)
def getBoundaryTraction(FEDEBoundMap,width=0.1):
   FENbc = dict((k,0) for k in FEDEBoundMap.keys())
   for key in FENbc.keys():
      DEid = FEDEBoundMap[key]
      f = Vector3.Zero; b = Omega().bodies[DEid]
      for intr in b.intrs():
         df = intr.phys.normalForce+intr.phys.shearForce
         if intr.id2 != DEid: df *=-1
         f = f+df
      dl = (b.shape.node1.state.pos-b.shape.node2.state.pos).norm()
      FENbc[key] = (f[0]/(dl*width),f[1]/(dl*width))
   #~ print max(FENbc.values())
   return FENbc
"""

### Below is Ning Guo's original script ###

def initLoad(ID=0): # where ID identifies the Gauss point location
   if 1:
      # All Gauss points import 'DE_int.yade.gz' resulting in a uniform sample (default)
      Omega().load(DE_int)
   else:
      # Otherwise load different packings to generate random field
      # resulting in an inherently heterogeneous sample
      Omega().load(str(ID)+'.yade.gz')
   Omega().tags['id']=str(ID)
   return Omega().sceneToString()
   
def outputPack(param):
   if len(param) != 3:
      raise RuntimeError,"No. of param. should be exactly 3. 0: RVE scene; 1: step count; 2: name prefix"
   Omega().stringToScene(param[0])
   pos = Omega().tags['id']
   Omega().save(param[2]+'packing_'+pos+'_'+str(param[1])+'.yade.gz')

def numOfParticles(scene):
   Omega().stringToScene(scene)
   return len(Omega().bodies) # !!! for spherical particle packing only

"""
# Apply deformation on 2D DEM packing (Ning's original function)
def shear2D(param):
   if len(param) != 2:
      raise RuntimeError,"No. of param. should be exactly 2. 0: RVE scene; 1: strain."
   Omega().stringToScene(param[0])
   # ns >= 2
   ns=int(max(1e5*numpy.max(numpy.abs(param[1])),2))
   dstrain = utils.Matrix3(param[1][0],param[1][1],0, param[1][2],param[1][3],0, 0,0,0)
   Omega().cell.velGrad=dstrain/(ns*Omega().dt)
   Omega().run(ns,True)
   Omega().cell.velGrad=utils.Matrix3.Zero
   return Omega().sceneToString()
"""
   
# Apply deformation on 3D DEM packing
def shear3D(param):
   if len(param) != 2:
      raise RuntimeError,"No. of param. should be exactly 2. 0: RVE scene; 1: strain."
   Omega().stringToScene(param[0])
   # ns >= 2
   ns=int(max(1e5*numpy.max(numpy.abs(param[1])),2))
   dstrain = utils.Matrix3(param[1][0],param[1][1],param[1][2], param[1][3],param[1][4],param[1][5], param[1][6],param[1][7],param[1][8])
   Omega().cell.velGrad=dstrain/(ns*Omega().dt)
   Omega().run(ns,True)
   Omega().cell.velGrad=utils.Matrix3.Zero
   return Omega().sceneToString()

# get updated stress tensor
# utils.getStress() is implemented in Shop.cpp
def getStress2D(scene):
   Omega().stringToScene(scene)
   s = utils.getStress()
   s = .5*(s+s.transpose())
   return [[s[0,0],s[0,1]],[s[1,0],s[1,1]]]
   
# get contact normal based fabric tensor
def getFabric2D(scene):
   Omega().stringToScene(scene)
   f = utils.fabricTensor(splitTensor=False,revertSign=False)[0]
   return [[f[0,0],f[0,1]],[f[1,0],f[1,1]]]

def getFabric3D(scene):
   Omega().stringToScene(scene)
   f = utils.fabricTensor(splitTensor=False,revertSign=False)[0]
   return [[f[0,0],f[0,1],f[0,2]],[f[1,0],f[1,1],f[1,2]],[f[2,0],f[2,1],f[2,2]]]

""" # Used for clumped particle model only
    # get particle orientation based fabric tensor
def getParOriFabric(scene):
   Omega().stringToScene(scene)
   fab = utils.Matrix3.Zero
   numPar = 0
   for b in Omega().bodies:
      if b.isClump:
         numPar += 1
         keys = b.shape.members.keys()
         pos1 = Omega().bodies[keys[0]].state.pos
         pos2 = Omega().bodies[keys[1]].state.pos
         ori = (pos1-pos2).normalized()
         fab += ori.outer(ori)
   fab /= numPar
   return [[fab[0],fab[1]],[fab[3],fab[4]]]
"""

""" # Used for cohesive particle model only
    # get the bond breakage number within a packing
def getDebondingNumber(param):
   Omega().stringToScene(param[0])
   num = 0
   for id1,id2 in param[1]:
      try:
         i = Omega().interactions[id1,id2]
         if not i.isReal():
            num += 1
         elif i.phys.cohesionBroken:
            num += 1
      except IndexError:
         num += 1
   return num
"""

# get updated stress tensor and tangent operator as a tuple
# utils.getStressAndTangent() is implemented in Shop.cpp
def getStressAndTangent2D(scene):
   Omega().stringToScene(scene)
   st = utils.getStressAndTangent(symmetry=True)
   s = st[0]
   s = .5*(s+s.transpose())
   t=[[[[0,0],[0,0]],[[0,0],[0,0]]],[[[0,0],[0,0]],[[0,0],[0,0]]]]
   t[0][0][0][0]=st[1][0,0]
   t[1][1][0][0]=t[0][0][1][1]=st[1][0,1]
   t[0][1][0][0]=t[0][0][0][1]=t[1][0][0][0]=t[0][0][1][0]=st[1][0,5]
   t[1][1][1][1]=st[1][1,1]
   t[1][1][0][1]=t[0][1][1][1]=t[1][1][1][0]=t[1][0][1][1]=st[1][1,5]
   t[0][1][0][1]=t[0][1][1][0]=t[1][0][0][1]=t[1][0][1][0]=st[1][5,5]
   return [[s[0,0],s[0,1]],[s[1,0],s[1,1]]],t

def getStressAndTangent3D(scene):
   Omega().stringToScene(scene)
   st = utils.getStressAndTangent(symmetry=True)
   s = st[0]
   s = .5*(s+s.transpose())
   t = numpy.zeros((3,3,3,3))
   t[0][0][0][0]=st[1][0,0]
   t[0][0][1][1]=t[1][1][0][0]=st[1][0,1]
   t[0][0][2][2]=t[2][2][0][0]=st[1][0,2]
   t[0][0][1][2]=t[0][0][2][1]=t[1][2][0][0]=t[2][1][0][0]=st[1][0,3]
   t[0][0][0][2]=t[0][0][2][0]=t[0][2][0][0]=t[2][0][0][0]=st[1][0,4]
   t[0][0][0][1]=t[0][0][1][0]=t[0][1][0][0]=t[1][0][0][0]=st[1][0,5]
   t[1][1][1][1]=st[1][1,1]
   t[1][1][2][2]=t[2][2][1][1]=st[1][1,2]
   t[1][1][1][2]=t[1][1][2][1]=t[1][2][1][1]=t[2][1][1][1]=st[1][1,3]
   t[1][1][0][2]=t[1][1][2][0]=t[0][2][1][1]=t[2][0][1][1]=st[1][1,4]
   t[1][1][0][1]=t[1][1][1][0]=t[0][1][1][1]=t[1][0][1][1]=st[1][1,5]
   t[2][2][2][2]=st[1][2,2]
   t[2][2][1][2]=t[2][2][2][1]=t[1][2][2][2]=t[2][1][2][2]=st[1][2,3]
   t[2][2][0][2]=t[2][2][2][0]=t[0][2][2][2]=t[2][0][2][2]=st[1][2,4]
   t[2][2][0][1]=t[2][2][1][0]=t[0][1][2][2]=t[1][0][2][2]=st[1][2,5]
   t[1][2][1][2]=t[1][2][2][1]=t[2][1][1][2]=t[2][1][2][1]=st[1][3,3]
   t[1][2][0][2]=t[1][2][2][0]=t[2][1][0][2]=t[2][1][2][0]=t[0][2][1][2]=t[2][0][1][2]=t[0][2][2][1]=t[2][0][2][1]=st[1][3,4]
   t[1][2][0][1]=t[1][2][1][0]=t[2][1][0][1]=t[2][1][1][0]=t[0][1][1][2]=t[1][0][1][2]=t[0][1][2][1]=t[1][0][2][1]=st[1][3,5]
   t[0][2][0][2]=t[2][0][0][2]=t[0][2][2][0]=t[2][0][2][0]=st[1][4,4]
   t[0][2][0][1]=t[0][2][1][0]=t[2][0][0][1]=t[2][0][1][0]=t[0][1][0][2]=t[1][0][0][2]=t[0][1][2][0]=t[1][0][2][0]=st[1][4,5]
   t[0][1][0][1]=t[0][1][1][0]=t[1][0][0][1]=t[1][0][1][0]=st[1][5,5]
   return [[s[0,0],s[0,1],s[0,2]],[s[1,0],s[1,1],s[1,2]],[s[2,0],s[2,1],s[2,2]]],t

# utils.voidratio2D() is implemented in Shop.cpp
# to get plane void ratio (for 2D test only)   
def getVoidRatio2D(scene):
   Omega().stringToScene(scene)
   zSize = Omega().cell.hSize[2,2]
   return utils.voidratio2D(zlen=zSize)

def getEquivalentPorosity(scene):
   Omega().stringToScene(scene)
   zSize = Omega().cell.hSize[2,2]
   e = utils.voidratio2D(zlen=zSize)
   return e/(1.+e)
   
def getVoidRatio3D(scene):
   Omega().stringToScene(scene)
   p = utils.porosity()
   return p/(1.0-p)

# get average rotation of the particles within a packing
def avgRotation2D(scene):
   Omega().stringToScene(scene)
   rot = 0.0
   for b in Omega().bodies:
      rot += b.state.rot()[2]
   rot /= len(Omega().bodies)
   return rot
   
def avgRotation3D(scene):
   Omega().stringToScene(scene)
   rot = utils.Vector3.Zero
   for b in Omega().bodies:
      rot += b.state.rot()
   rot /= len(Omega().bodies)
   return [rot[0],rot[1],rot[2]]
