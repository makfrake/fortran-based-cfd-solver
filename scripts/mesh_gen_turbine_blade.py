import gmsh
import sys


gmsh.initialize()

# Characteristic length

lc = 1                  # Target mesh size
passo = 0.8333          # Passo della schiera

# Point definition

p1 = gmsh.model.geo.addPoint(-1., 0.5,0, lc)
p2 = gmsh.model.geo.addPoint(0.5, 0.5, 0, lc)
p3 = gmsh.model.geo.addPoint(1., -0.2, 0, lc)
p4 = gmsh.model.geo.addPoint(2.,-0.2, 0, lc)
p5 = gmsh.model.geo.addPoint(2., -0.2-passo, 0, lc)
p6 = gmsh.model.geo.addPoint(1., -0.2-passo, 0, lc)
p7 = gmsh.model.geo.addPoint(0.5, 0.5-passo, 0, lc)
p8 = gmsh.model.geo.addPoint(-1., 0.5-passo, 0, lc)

# Lines definition

l1  = gmsh.model.geo.addLine(p1,p2,1)
l2  = gmsh.model.geo.addLine(p2,p3,2)
l3  = gmsh.model.geo.addLine(p3,p4,3)
l4  = gmsh.model.geo.addLine(p4,p5,4)
l5  = gmsh.model.geo.addLine(p5,p6,5)
l6  = gmsh.model.geo.addLine(p6,p7,6)
l7  = gmsh.model.geo.addLine(p7,p8,7)
l8  = gmsh.model.geo.addLine(p8,p1,8)

# 2D Surface generation

loop    = gmsh.model.geo.addCurveLoop([l1,l2,l3,l4,l5,l6,l7,l8],1)
surface = gmsh.model.geo.addPlaneSurface([loop],1)

# Unstructured mesh generation

field1 = gmsh.model.mesh.field.add("Attractor",1)
gmsh.model.mesh.field.setNumbers(field1, "EdgesList",[l3,l4,l5,l6,l7,l8])

field2 = gmsh.model.mesh.field.add("Threshold", 2)
gmsh.model.mesh.field.setNumber(field2, "IField",1)
gmsh.model.mesh.field.setNumber(field2, "LcMin",lc*0.005)
gmsh.model.mesh.field.setNumber(field2, "IField",lc*0.025)
gmsh.model.mesh.field.setNumber(field2, "DistMin",h/4)
gmsh.model.mesh.field.setNumber(field2, "DistMax",h)

gmsh.model.mesh.field.setAsBackgroundMesh(2)
gmsh.model.geo.mesh.setRecombine(1,1)

gmsh.model.geo.synchronize()

# Boundary conditions

gmsh.model.addPhysicalGroup(2,[surface],1)
gmsh.model.addPhysicalGroup(1,[l1,l2,l3,l4,l5,l7,l8,l10],99,"wall")
gmsh.model.addPhysicalGroup(1,[l11],2,"inlet")
gmsh.model.addPhysicalGroup(1,[l6,l9],3,"outlet")

# Generate mesh

gmsh.model.mesh.generate(2)

# Save msh file

gmsh.write("intake.msh")

if '-nopopup' not in sys.argv:
    gmsh.fltk.run()

gmsh.finalize()