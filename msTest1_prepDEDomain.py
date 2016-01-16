""" Author: Hongyang Cheng <chyalexcheng@gmail>
    Test #1: prepare membrane strucutre using chained cylinders
"""

# import modules
from yade import qt, plot
import numpy as np

#####################
##  Key controlls  ##
#####################

# mesh file name
mshName = 'Msh2'
# sample size, 1.2m by 1.2m
lx = 1.2; ly = 1.2
# confining pressure
confining = 0
# pulling speed on ends of DE membrane
pullSpeed = -1.e-2
# type of membrane material
GSType = 'PP'
color = [84./255,89./255,109./255]
# global damping
damp = 0.2
width = 0.1
# assumed radius
rGrid = 5.e-3
# discretization per cylinder
L = .3; nL = 0
if not nL: nL = int(L/(2.*rGrid))
# factor for greater GridCo-GridCo stiffness
stif = 1e0
# density
dScaling = 1e0
rho = dScaling*443.96

# tensile modulus of membrane material
if GSType == 'PE':
   maxStrainStress = [0.17035, 37.971e6]
   thick = .251e-3
if GSType == 'PP':
   maxStrainStress = [0.1425, 80.8e6]
   thick = .393e-3
young = maxStrainStress[1]/maxStrainStress[0]
sigTmax = maxStrainStress[1]

## material parameters for internal behavior
# consider 10 cm width to compute 2D modulus and strength
young *= width
sigTmax *= width
# rescale 2D modulus and strength to fit with assumed radius, !!! obmit pi here (code bug)
young *= thick/(pi*rGrid**2)
sigTmax *= thick/(pi*rGrid**2)
# m: membrane, i:interface
E_m = young; v_m = 0.; phi_m = 0.; sigTmax_m = sigTmax; sigSmax_m = sigTmax
E_i = 0.   ; v_i = 0.; phi_i = 0.; sigTmax_i = sigTmax; sigSmax_i = sigTmax

## material parameters for external behavior
# m2i: membrane-interface
#~ E_m2i = stif*young; v_m2i = 0.33; phi_m2i = radians(34)
E_m2i = stif*young; v_m2i = 0.33; phi_m2i = radians(0)

#################
##  Functions  ##
#################

# sort the node ids in anti-clockwise order
def antiClockSort(ids,pos0):
   mn,mx = aabbExtrema()
   coordC = .5*(mn+mx)
   for i in ids:
      b = O.bodies
      coord = b[i].state.pos
      dis = coord-coordC
      phi = acos(dis[0]/dis.norm())
      if dis[1]/dis.norm() < 0:
         phi = 2.*pi-phi
      ids[ids.index(i)] = [i,phi]
      if coord == pos0:
         startID = i
   ids=sorted(ids, key=lambda x:x[1])
   # return only the sorted ids
   for i in xrange(len(ids)):
      ids[i] = ids[i][0]
   # ID start from [0,0,0]
   index = ids.index(startID)
   ids = ids[index:]+ids[:index]
   return ids

# add gridNode body to O.bodies
def gridNode(center,radius,dynamic=None,fixed=False,wire=False,color=None,highlight=False,material=-1):
   b=Body()
   b.shape=GridNode(radius=radius,color=color if color else utils.randomColor(),wire=wire,highlight=highlight)
   #V=(4./3)*math.pi*radius**3   # will be overriden by the connection
   V=0.
   geomInert=(2./5.)*V*radius**2   # will be overriden by the connection
   utils._commonBodySetup(b,V,Vector3(geomInert,geomInert,geomInert),material,pos=center,dynamic=dynamic,fixed=fixed)
   b.aspherical=False
   b.bounded=False
   b.mask=0   #avoid contact detection with the nodes. Manual interaction will be set for them in "gridConnection" below.
   return b

# add gridConnection body and interaction to O.bodies
def gridConnection(id1,id2,radius,wire=False,color=None,highlight=False,material=-1,mask=1,cellDist=None):
   b=Body()
   b.shape=GridConnection(radius=radius,color=color if color else utils.randomColor(),wire=wire,highlight=highlight)
   sph1=O.bodies[id1] ; sph2=O.bodies[id2]
   i=createInteraction(id1,id2)
   nodeMat=sph1.material
   b.shape.node1=sph1 ; b.shape.node2=sph2
   sph1.shape.addConnection(b) ; sph2.shape.addConnection(b)
   if(O.periodic):
      if(cellDist!=None):
         i.cellDist=cellDist
      segt=sph2.state.pos + O.cell.hSize*i.cellDist - sph1.state.pos
   else: segt=sph2.state.pos - sph1.state.pos
   L=segt.norm()
   V=0.5*L*math.pi*radius**2
   geomInert=(2./5.)*V*radius**2
   utils._commonBodySetup(b,V,Vector3(geomInert,geomInert,geomInert),material,pos=sph1.state.pos,dynamic=False,fixed=True)
   sph1.state.mass = sph1.state.mass + V*nodeMat.density
   sph2.state.mass = sph2.state.mass + V*nodeMat.density
   for k in [0,1,2]:
      sph1.state.inertia[k] = sph1.state.inertia[k] + geomInert*nodeMat.density
      sph2.state.inertia[k] = sph2.state.inertia[k] + geomInert*nodeMat.density
   b.aspherical=False
   if O.periodic:
      i.phys.unp= -(sph2.state.pos + O.cell.hSize*i.cellDist - sph1.state.pos).norm() + sph1.shape.radius + sph2.shape.radius
      b.shape.periodic=True
      b.shape.cellDist=i.cellDist
   else:
      i.phys.unp= -(sph2.state.pos - sph1.state.pos).norm() + sph1.shape.radius + sph2.shape.radius   
   i.geom.connectionBody=b
   I=math.pi*(2.*radius)**4/64.
   E=nodeMat.young
   i.phys.kn=E*math.pi*(radius**2)/L
   # use the correct normalAdhesion and shear Adhesion
   i.phys.normalAdhesion = nodeMat.normalCohesion*math.pi*(radius**2)
   i.phys.shearAdhesion = nodeMat.shearCohesion*math.pi*(radius**2)
   # do not use shear and rolling/bending stiffness
   #~ i.phys.kr=E*I/L
   #~ i.phys.ks=12.*E*I/(L**3)
   #~ G=E/(2.*(1+nodeMat.poisson))
   #~ i.phys.ktw=2.*I*G/L
   b.mask=mask
   return b
   
###############
##  Engines  ##
###############

## Engines need to be defined first since the function gridConnection creates the interaction
O.engines=[
   ForceResetter(),
   InsertionSortCollider([
      Bo1_GridConnection_Aabb(),
   ]),
   InteractionLoop(
      [Ig2_GridNode_GridNode_GridNodeGeom6D(),
       Ig2_GridConnection_GridConnection_GridCoGridCoGeom(),
       ],
      [Ip2_CohFrictMat_CohFrictMat_CohFrictPhys(
         setCohesionNow=False,
         setCohesionOnNewContacts=False),
       Ip2_FrictMat_FrictMat_FrictPhys(),
       ],
      [Law2_ScGeom6D_CohFrictPhys_CohesionMoment(),
       Law2_GridCoGridCoGeom_FrictPhys_CundallStrack(),
       ]
   ),
   NewtonIntegrator(gravity=(0,0,0),damping=damp,label='newton')
]

#################
##  Materials  ##
#################

# interface
O.materials.append(CohFrictMat(young=E_i,poisson=v_i,density=rho,frictionAngle=phi_i \
   ,normalCohesion=sigTmax_i,shearCohesion=sigSmax_i,momentRotationLaw=False,label='iMat'))
# membrane
O.materials.append(CohFrictMat(young=E_m,poisson=v_m,density=rho,frictionAngle=phi_m \
   ,normalCohesion=sigTmax_m,shearCohesion=sigSmax_m,momentRotationLaw=False,label='mMat'))
# membrane to interface
O.materials.append(FrictMat(young=E_m2i,poisson=v_m2i,density=rho,frictionAngle=phi_m2i,label='m2iMat'))
# boundary wall
O.materials.append(FrictMat(young=E_m2i,poisson=v_m2i,density=rho,frictionAngle=0,label='walMat'))

###################
##  Build model  ##
###################

# load the boundary nodes from FE domain
bx = np.load('Bx'+mshName+'.npy').item()
bRefIDs = np.load('BRefID'+mshName+'.npy').item()

##~~~~~~~~~~~~~~~~~~~~~~##
##  generate interface  ##
##~~~~~~~~~~~~~~~~~~~~~~##

# create all GridNodes first
nodeIDs = {}
for key in bx.keys():
   x,y = bx[key]
   bx[key] = O.bodies.append(
      gridNode([x,y,0],rGrid,wire=True,fixed=True,material='iMat',color=[1.,1.,0.]))
   if y != 0:
		nodeIDs[key] = bx[key]
# create connection between the nodes
startPos = Vector3.Zero
iNodesIds = antiClockSort(bx.values(),startPos)
# save interface element IDs
iBodiesIds = []
for i,j in zip(iNodesIds,np.append(iNodesIds[1:],iNodesIds[0])):
   iBodiesIds.append(O.bodies.append(gridConnection(i,j,rGrid,wire=True,material='m2iMat',color=[1.,1.,0.])))
# replace FE boundary reference ID with DE interface element IDs (numbering orders are consistent)
dID = min(iBodiesIds)-min(bRefIDs.values())
for key in bRefIDs: bRefIDs[key] += dID

##~~~~~~~~~~~~~~~~~~~~~##
##  generate membrane  ##
##~~~~~~~~~~~~~~~~~~~~~##

# create all GridNodes first
mNodesIds = []
# primary nodes adjacent to interface nodes
for i in iNodesIds[4:]+[iNodesIds[0]]:
   i1,i2 = O.interactions.withBody(i)
   # create membrane nodes without overlap with interface nodes
   n1 = i1.geom.normal; n2 = i2.geom.normal
   alpha = .5*acos(n1.dot(n2))
   pos = O.bodies[i].state.pos
   dPos = 2.*rGrid/cos(alpha)
   m1 = n1.cross(Vector3(0,0,1))
   m2 = n2.cross(Vector3(0,0,1))
   m0 = (m1+m2).normalized()
   pos_new = pos + dPos*m0
   mNodesIds.append(O.bodies.append(
      gridNode(pos_new,rGrid,wire=True,fixed=False,material='mMat',color=[0.,1.,0.])))
# refine membrane elements
for i,j in zip(mNodesIds[:-1],mNodesIds[1:]):
   dL = (O.bodies[j].state.pos-O.bodies[i].state.pos)/nL
   for k in xrange(nL-1):
      pos = O.bodies[i].state.pos+(k+1)*dL
      mNodesIds.append(O.bodies.append(
         gridNode(pos,rGrid,wire=True,fixed=False,material='mMat',color=[0.,1.,0.])))
# create connection between the nodes
startPos = O.bodies[mNodesIds[0]].state.pos
mNodesIds = antiClockSort(mNodesIds,startPos)
for i,j in zip(mNodesIds[:-1],mNodesIds[1:]):
   O.bodies.append(gridConnection(i,j,rGrid,wire=True,material='m2iMat',color=[0.,1.,0.]))
# block z-direction translational and rotational DOFs of membrane nodes
for i in mNodesIds: O.bodies[i].state.blockedDOFs = 'zXYZ'

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
##  estimate size of timestep       ##
##  apply boundary force if any     ##
##  apply boundary velocity if any  ##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##

# estimate timestep size
O.dt = 2*utils.PWaveTimeStep()
print O.dt
np.save('FEDENodeMap'+mshName+'.npy',nodeIDs)
np.save('FEDEBoundMap'+mshName+'.npy',bRefIDs)
np.save('mNodesIds'+mshName+'.npy',mNodesIds)

# apply boundary force
if confining:
   for i in mNodesIds:
      intrs = O.bodies[i].intrs()
      segts = []; n = utils.Vector3.Zero
      for intr in intrs:
         if isinstance(intr.phys,CohFrictPhys):
            segt = (intr.geom.refR1+intr.geom.refR2-intr.geom.penetrationDepth)*intr.geom.normal
            segts.append(segt); n += segt
      # segment normal
      n = n.normalized()
      dl = sum([.5*segt.dot(n) for segt in segts])
      # cross product of unit z and interaction normal points outward
      n = n.cross(utils.Vector3.UnitZ)
      # get boundary force per node: f = sigma*segt*width (negative sign for compression)
      permDEf = confining*(dl*width)*n
      # add permanent force on membrane nodes
      O.forces.addF(i,permDEf,True)

# apply boundary velocity
if pullSpeed:
	endNodes = [mNodesIds[0],mNodesIds[-1]]
	for i in range(len(endNodes)):
		endID = endNodes[i]
		O.bodies[endID].state.blockedDOFs = 'xyzXYZ'
		O.bodies[endID].state.vel = Vector3(0,pullSpeed,0)
		pos0 = O.bodies[endID].state.pos+Vector3(2*rGrid*(-1)**i,0,0)
		pos1 = pos0+Vector3(0,ly+4*rGrid,0)
		id0 = O.bodies.append(
			gridNode(pos0,rGrid,wire=True,fixed=True,material='iMat',color=[0.,0.,1.]))
		id1 = O.bodies.append(
			gridNode(pos1,rGrid,wire=True,fixed=True,material='iMat',color=[0.,0.,1.]))
		O.bodies.append(gridConnection(id0,id1,rGrid,wire=True,material='walMat',color=[0.,0.,1.]))

while 1:
	O.run(100,True)
	if O.forces.f(8).norm() > 1e-2:
			break

# save exterior DE scene
qt.Renderer().bgColor = color
O.save('DE_ext_32threads.yade.gz')
