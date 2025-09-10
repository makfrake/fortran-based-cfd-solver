import gmsh
import sys
from turbine_blade_points import s9


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

loop1    = gmsh.model.geo.addCurveLoop([l1,l2,l3,l4,l5,l6,l7,l8],1)
loop2    = gmsh.model.geo.addCurveLoop(s9,2)

surface = gmsh.model.geo.addPlaneSurface([loop1,loop2],1)

# Periodic boundaries



# Unstructured mesh generation

field1 = gmsh.model.mesh.field.add("Attractor",1)
gmsh.model.mesh.field.setNumbers(field1, "EdgesList",s9)

field2 = gmsh.model.mesh.field.add("Threshold", 2)
gmsh.model.mesh.field.setNumber(field2, "IField",1)
gmsh.model.mesh.field.setNumber(field2, "LcMin",0.01)
gmsh.model.mesh.field.setNumber(field2, "IField",0.05)
gmsh.model.mesh.field.setNumber(field2, "DistMin",0.25)
gmsh.model.mesh.field.setNumber(field2, "DistMax",0.5)

field3 = gmsh.model.mesh.field.add("Min")
gmsh.model.mesh.field.setNumbers(field3, "FieldsList",l2)
gmsh.model.mesh.field.setAsBackgroundMesh(3)

field4 = gmsh.model.mesh.field.add("BoundaryLayer",4)
gmsh.model.mesh.field.setNumber(field4, "EdgesList",s9)
gmsh.model.mesh.field.setNumber(field4, "hwall_n",0.001)
gmsh.model.mesh.field.setNumber(field4, "thickness",0.05)
gmsh.model.mesh.field.setNumber(field4, "ratio",1.2)
gmsh.model.mesh.field.setNumber(field4, "Quads",1)
gmsh.model.mesh.field.setAsBoundaryLayer(4)

#gmsh.model.geo.mesh.setRecombine(1,1)

gmsh.model.geo.synchronize()

# Boundary conditions

gmsh.model.addPhysicalGroup(2,[surface],1)
gmsh.model.addPhysicalGroup(1,s9,99,"blade")
gmsh.model.addPhysicalGroup(1,[l5,l6,l7],50,"")
gmsh.model.addPhysicalGroup(1,[l1,l2,l3],40,"")
gmsh.model.addPhysicalGroup(1,[l8],2,"inlet")
gmsh.model.addPhysicalGroup(1,[l4],3,"outlet")

# Generate mesh

gmsh.model.mesh.generate()

# Save msh file

gmsh.write("turbine_blade.msh")

if '-nopopup' not in sys.argv:
    gmsh.fltk.run()

gmsh.finalize()