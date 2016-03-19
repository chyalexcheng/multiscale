from yade import pack

random = False
pressure = -2.e4
sceneName = 'Test2'
spMat = O.materials.append(FrictMat(young=7.5e8,poisson=.016,frictionAngle=radians(16)))
sp = pack.SpherePack()
if random:
	size = .24
	O.materials[0].frictionAngle = radians(0)
	sp.makeCloud(minCorner=(0,0,.05),maxCorner=(size,size,.05),\
				    psdSizes=[.003,.005,.01],\
					 psdCumm =[.2  ,.85 , 1.],\
					 num=400, periodic=True)
	sp.toSimulation(material=spMat)
	for b in O.bodies:
		if b.shape.radius<=.003/2:
			b.shape.radius=.003/2
		elif .003/2<b.shape.radius<=.005/2:
			b.shape.radius=.005/2
		elif .005/2<b.shape.radius<=.01/2:
			b.shape.radius=.01/2
	O.cell.hSize = Matrix3(size,0,0, 0,size,0, 0,0,.1)
else:
	sp.load('alumRods0.txt')
	sp.toSimulation(material=spMat)
	O.periodic=True
	O.cell.hSize = Matrix3(0.11280356495985018,0,0, 0,0.10837027507783764,0, 0,0,0.1)
print len(O.bodies)

for p in O.bodies:
   p.state.blockedDOFs = 'zXY'
   p.state.mass = 2700 * 0.1 * pi * p.shape.radius**2
   inertia = 0.5 * p.state.mass * p.shape.radius**2
   p.state.inertia = (.5*inertia,.5*inertia,inertia)

O.dt = utils.PWaveTimeStep()
print O.dt

O.engines = [
   ForceResetter(),
   InsertionSortCollider([Bo1_Sphere_Aabb()]),
   InteractionLoop(
      [Ig2_Sphere_Sphere_ScGeom()],
      [Ip2_FrictMat_FrictMat_FrictPhys()],
      [Law2_ScGeom_FrictPhys_CundallStrack()]
   ),
   PeriTriaxController(
      dynCell=True,
      goal=(0,pressure,0),
      stressMask=2,
      relStressTol=.001,
      absStressTol=1.,
      maxUnbalanced=.001,
      maxStrainRate=(.01,.001,.0),
      doneHook='term()',
      label='biax'
   ),
   NewtonIntegrator(damping=.2)
]

def term():
   O.engines = O.engines[:3]+O.engines[4:]
   print getStress()
   print O.cell.hSize
   setContactFriction(radians(16))
   O.cell.trsf=Matrix3.Identity
   O.cell.velGrad=Matrix3.Zero
   for p in O.bodies:
      p.state.vel = Vector3.Zero
      p.state.angVel = Vector3.Zero
      p.state.refPos = p.state.pos
      p.state.refOri = p.state.ori
   O.save('./DE_exts/'+sceneName+'/DE_alum.yade.gz')
   sp.fromSimulation()
   sp.save('./DE_exts/'+sceneName+'/DE_alum.txt')
   O.pause()

O.run()

