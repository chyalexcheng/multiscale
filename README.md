# multiscale
Multiscale modeling approach for soil-geosynthetic interaction
msTest1: shape-forming test
msTest2: pull-out test

# Escript functions
msFEM2DExplicit.py
msFEM2DExplicit.pyc

# Yade functions
simDEM.py
simDEM.pyc

# other functions
saveGauss.py
saveGauss.pyc

# Yade script to build external scene
msTest1_prepDEDomain.py
geosynthetic inclusion is continously created on the left, right and bottom boundaries. Load is applied by pulling both left and right ends downward.
 |--------------|
 |		|
 |		|
 |		|
 |		|
 |		|
\|/	       \|/
 
msTest2_prepDEDomain.py
geosynthetic inclusion is embedded beneath granular soil in a rectangular box. Load is appllied by horizontally pulling the right end of geosynthetic inclusion.

 ----------------------->

# Yade script to build internal scene
prepareAlumRods.py
create DEM particle packings as RVEs, considering material properties of aluminium rods
Initial stress conditions are:
    shape-forming test	: stress-free condition
    pull-out test	: at-rest state, surcharge pressure = 20 kPa

# main script
msTest1_implicit.py
msTest2_explicit.py

# profiler
msTest1.profile

# basic initial packing
alumRods0.txt

# files to build external DE scene
BxMsh1.npy	BxMsh2.npy	BxMsh3.npy	BxMsh4.npy
BRefIDMsh1.npy	BRefIDMsh2.npy	BRefIDMsh3.npy	BRefIDMsh4.npy

# msh files
Msh1.msh	Msh2.msh	Msh3.msh	Msh4.msh
