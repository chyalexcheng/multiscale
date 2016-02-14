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
# initial pressure on membrane from soil
pressure = -2.e4
# pulling speed on ends of DE membrane
pullSpeed = 1.e-2
# type of membrane material
GSType = 'PP'
color = [84./255,89./255,109./255]
# global damping
damp = 0.2
width = 0.1
# assumed radius
rGrid = 5.e-3
# discretization per cylinder (change this when mesh changes)
L = .3; nL = 0
if not nL: nL = int(L/(2.*rGrid))
rGrid = 1.e-3
# factor for greater GridCo-GridCo stiffness
stif = 1e0
# density
dScaling = 1e3
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
E_m2i = stif*young; v_m2i = 0.33; phi_m2i = radians(34)
#~ E_m2i = stif*young; v_m2i = 0.33; phi_m2i = radians(0)

#################
##  Functions  ##
#################

# sort the node ids in anti-clockwise order
def antiClockSort(ids,pos0,coordC):
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
   NewtonIntegrator(gravity=(0,0,0),damping=damp)
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

###################
##  Build model  ##
###################

# load the boundary nodes from FE domain
bx = np.load('Bx'+mshName+'.npy').item()

##~~~~~~~~~~~~~~~~~~~~~~##
##  generate interface  ##
##~~~~~~~~~~~~~~~~~~~~~~##

# create all GridNodes first
nodeIDs = {}
for key in bx.keys():
   x,y = bx[key]
   if y == 0:
      nodeID = O.bodies.append(
         gridNode([x,y,0],rGrid,wire=True,fixed=True,material='iMat',color=[1.,1.,0.]))
      nodeIDs[key] = nodeID
# create connection between the nodes
startPos = Vector3.Zero
coords = np.transpose(bx.values())
coordC = Vector3(np.mean(coords[0]),np.mean(coords[1]),0)
iNodesIds = antiClockSort(nodeIDs.values(),startPos,coordC)
# save interface element IDs
iBodiesIds = []
for i,j in zip(iNodesIds,iNodesIds[1:]):
   iBodiesIds.append(O.bodies.append(gridConnection(i,j,rGrid,wire=True,material='m2iMat',color=[1.,1.,0.])))

##~~~~~~~~~~~~~~~~~~~~~##
##  generate membrane  ##
##~~~~~~~~~~~~~~~~~~~~~##

# create all GridNodes first
mNodesIds = []
# primary nodes adjacent to interface nodes (change this when mesh changes)
for i in iNodesIds:
   pos_new = O.bodies[i].state.pos-Vector3(0,2.*rGrid,0)
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
mNodesIds = antiClockSort(mNodesIds,startPos,coordC)
for i,j in zip(mNodesIds[:-1],mNodesIds[1:]):
   O.bodies.append(gridConnection(i,j,rGrid,wire=True,material='m2iMat',color=[0.,1.,0.]))
# block z-direction translational and rotational DOFs of membrane nodes
for i in mNodesIds: O.bodies[i].state.blockedDOFs = 'xzXYZ'

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
##  estimate size of timestep       ##
##  apply boundary force if any     ##
##  apply boundary velocity if any  ##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##

# estimate timestep size
O.dt = utils.PWaveTimeStep()
print O.dt
np.save('./DE_exts/Test2/FEDENodeMap'+mshName+'.npy',nodeIDs)
np.save('./DE_exts/Test2/mNodesIds'+mshName+'.npy',mNodesIds)

# create bottom boundary platen
ends = [mNodesIds[0],mNodesIds[-1]]
for i in range(2):
	pos = O.bodies[ends[i]].state.pos-Vector3(0,2.*rGrid,0)
	ends[i] = O.bodies.append(
		gridNode(pos,rGrid,wire=True,fixed=True,material='iMat',color=[1.,1.,0.]))
O.bodies.append(gridConnection(ends[0],ends[1],rGrid,wire=True,material='m2iMat',color=[1.,1.,0.]))

# apply initial force on membrane
f = Vector3(0,pressure*lx/len(iBodiesIds)*width,0)
O.engines[-1].damping = 0.99
for i in iNodesIds[1:-1]:
	# add permanent force on interface nodes
	O.forces.addF(i,f,True)
	O.bodies[i].state.blockedDOFs = 'xzXYZ'
for i in [iNodesIds[0],iNodesIds[-1]]:
	# add permanent force on interface nodes
	O.forces.addF(i,.5*f,True)
	O.bodies[i].state.blockedDOFs = 'xzXYZ'
# wait till simulation stablized
while 1:
	calm(); O.run(10000,True)
	if kineticEnergy()<1e-5:
		break

# fix interface nodes
for i in iNodesIds: O.bodies[i].state.blockedDOFs = 'xyzXYZ'
# free x direction DOF of membrane Nodes
for i in mNodesIds: O.bodies[i].state.blockedDOFs = 'yzXYZ'
# reset forces and kinematics of all bodies
O.forces.reset()
for b in O.bodies:
	b.state.vel = Vector3.Zero
	b.state.angVel = Vector3.Zero
	b.state.refPos = b.state.pos
	b.state.refOri = b.state.ori
# set damping to normal level
O.engines[-1].damping = 0.2

# apply pull-out
O.bodies[mNodesIds[-1]].state.vel = Vector3(pullSpeed,0,0)
O.bodies[mNodesIds[-1]].state.blockedDOFs = 'xyzXYZ'

# save exterior DE scene
qt.Renderer().bgColor = color
O.save('./DE_exts/Test2/DE_ext_'+mshName+'.yade.gz')
