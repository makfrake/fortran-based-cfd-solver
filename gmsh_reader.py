import gmsh
import numpy as np
import sys


gmsh.initialize()

# Characterich length generation

l  = 1
h  = 0.3
lc = 0.005   # Target mesh size

# Point definition

p1 = gmsh.model.geo.addPoint(0,0,0,lc,1)
p2 = gmsh.model.geo.addPoint(l,0,0,lc,2)
p3 = gmsh.model.geo.addPoint(l,h,0,lc,3)
p4 = gmsh.model.geo.addPoint(0,h,0,lc,4)
p5 = gmsh.model.geo.addPoint(l/2,h/10,0,lc,5)
p6 = gmsh.model.geo.addPoint(l/2-l/8,0,0,lc,6)
p7 = gmsh.model.geo.addPoint(l/2+l/8,0,0,lc,7)

# Lines definitions

l1 = gmsh.model.geo.addLine(p2,p3,1)
l2 = gmsh.model.geo.addLine(p3,p4,2)
l3 = gmsh.model.geo.addLine(p4,p1,3)
s4 = gmsh.model.geo.addSpline([p1,p6,p5,p7,p2],4)

# 2D Surface generation

loop    = gmsh.model.geo.addCurveLoop([l1,l2,l3,s4],1)
surface = gmsh.model.geo.addPlaneSurface([loop],1)

# Structured mesh generation

transf1 = gmsh.model.geo.mesh.setTransfiniteCurve(1, 120, "Progression", 1)
transf2 = gmsh.model.geo.mesh.setTransfiniteCurve(3, 120, "Progression", -1)
transf3 = gmsh.model.geo.mesh.setTransfiniteCurve(2, 400, "Bump", 1)
transf4 = gmsh.model.geo.mesh.setTransfiniteCurve(4, 400, "Bump", 1)

gmsh.model.geo.mesh.setTransfiniteSurface(surface,"Left",[p1,p2,p3,p4])

# Quad mesh
gmsh.model.geo.mesh.setRecombine(1,1)

gmsh.model.geo.synchronize()

# Boundary conditions

gmsh.model.addPhysicalGroup(2,[l1],1)
gmsh.model.addPhysicalGroup(1,[s4,l2],99,"wall")
gmsh.model.addPhysicalGroup(1,[l3],2,"inlet")
gmsh.model.addPhysicalGroup(1,[l1],3,"outlet")

# Save msh file

gmsh.write("bump.msh")