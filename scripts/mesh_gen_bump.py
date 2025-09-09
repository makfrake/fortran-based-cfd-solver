import gmsh
import sys


gmsh.initialize()

# Characteristic length

l  = 1
h  = 0.3
lc = 0.01               # Target mesh size
hp = int(round(l/lc))   # Number of horizontal grid points
vp = int(round(h/lc))   # Number of horizontal grid points

# Point definition

p1 = gmsh.model.geo.addPoint(0,0,0,lc,1)
p2 = gmsh.model.geo.addPoint(l,0,0,lc,2)
p3 = gmsh.model.geo.addPoint(l,h,0,lc,3)
p4 = gmsh.model.geo.addPoint(0,h,0,lc,4)
p5 = gmsh.model.geo.addPoint(l/2,h/10,0,lc,5)
p6 = gmsh.model.geo.addPoint(l/2-l/8,0,0,lc,6)
p7 = gmsh.model.geo.addPoint(l/2+l/8,0,0,lc,7)

# Lines definition

l1 = gmsh.model.geo.addLine(p2,p3,1)
l2 = gmsh.model.geo.addLine(p3,p4,2)
l3 = gmsh.model.geo.addLine(p4,p1,3)
s4 = gmsh.model.geo.addSpline([p1,p6,p5,p7,p2],4)

# 2D Surface generation

loop    = gmsh.model.geo.addCurveLoop([l1,l2,l3,s4],1)
surface = gmsh.model.geo.addPlaneSurface([loop],1)

# Structured mesh generation

transf1 = gmsh.model.geo.mesh.setTransfiniteCurve(l1, vp, "Progression", 1)
transf2 = gmsh.model.geo.mesh.setTransfiniteCurve(l3, vp, "Progression", -1)
transf3 = gmsh.model.geo.mesh.setTransfiniteCurve(l2, hp, "Bump", 1)
transf4 = gmsh.model.geo.mesh.setTransfiniteCurve(s4, hp, "Bump", 1)

gmsh.model.geo.mesh.setTransfiniteSurface(surface,"Left",[p1,p2,p3,p4])

# Quad mesh
gmsh.model.geo.mesh.setRecombine(1,1)

gmsh.model.geo.synchronize()

# Boundary conditions

gmsh.model.addPhysicalGroup(2,[surface],1)
gmsh.model.addPhysicalGroup(1,[s4,l2],99,"wall")
gmsh.model.addPhysicalGroup(1,[l3],2,"inlet")
gmsh.model.addPhysicalGroup(1,[l1],3,"outlet")

# Generate mesh

gmsh.model.mesh.generate(2)

# Save msh file

gmsh.write("bump.msh")

if '-nopopup' not in sys.argv:
    gmsh.fltk.run()

gmsh.finalize()
