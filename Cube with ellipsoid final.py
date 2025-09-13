# author Rohit Kumar Yadav
# python script to create representative model for AM MarM 509 heat treated at 1250C

from abaqus import *
from abaqusConstants import *
import __main__

import section
import regionToolset
import displayGroupMdbToolset as dgm
import part
import material
import assembly
import step
import interaction
import load
import mesh
import optimization
import job
import sketch
import visualization
import xyPlot
import displayGroupOdbToolset as dgo
import connectorBehavior
import math as m
import random as rm
import numpy as np
import itertools


# function to find distance
def dist(x1, y1, z1=0, x2=0, y2=0, z2=0):
    d = m.sqrt((x2 - x1) * (x2 - x1) + (y2 - y1) * (y2 - y1) + (z2 - z1) * (z2 - z1))
    return d

#function to find volume of ellipse
def ellipsoidVol(a = 0.0, b = 0.0):
    vol = 4 * m.pi * a * b * b / 3
    return vol

#function to find corners of cube used to find cells of carbides
def cube_corners(x, y, z, a):
    half = a / 2
    offsets = [-half, half]
    corners = [(x + dx, y + dy, z + dz) for dx, dy, dz in itertools.product(offsets, repeat=3)]
    corners = tuple(corners)
    return corners

#domain parameters;

cubeEdge = 5.0  #um     #edge of cube
rqdvf = 0.02     #required volume fraction of carbide/total volume


cubeVol = cubeEdge * cubeEdge * cubeEdge

#creating carbide list
carbidevf = 0.0
noOfCarbides = 0.0
carbideVol = 0.0
x = []      #x-coordinate
y = []      #y-coordinate
z = []      #z-coordinate of ellipsoid center
a = []      #semi-major axis of ellipsoid
b = []      #semi-minor axis of ellipsoid
carbideCellList = ()    #corner regions of ellipsoid to assign section later

while carbidevf <= rqdvf:
    ap = rm.uniform(0.6,1.3)
    bp = rm.uniform(0.25,0.40)
    xp = rm.uniform(-cubeEdge / 2 + bp,cubeEdge / 2 - bp)
    yp = rm.uniform(-cubeEdge / 2 + ap,cubeEdge / 2 - ap)
    zp = rm.uniform(-cubeEdge / 2 + bp,cubeEdge / 2 - bp)
    inter = 0
    for j in range(len(x)):
        dxz = dist(x[j], 0, z[j], xp, 0, zp)
        dy = dist(0, y[j], 0, 0, yp, 0)
        if dxz < bp + b[j] and dy < ap + a[j]:
            inter = 1
            break
    if inter == 0:
        a.append(ap)
        b.append(bp)
        x.append(xp)
        y.append(yp)
        z.append(zp)
        noOfCarbides = noOfCarbides + 1
        print('Carbide = ',noOfCarbides)
        carbideVol = carbideVol + ellipsoidVol(ap,bp)
        print('Center = ', xp,yp,zp)
        carbidevf = carbideVol/cubeVol
        print('a = ',ap,' b = ',bp)
        print('Carbide volume fraction = ',carbidevf)
        carbideCellList+= cube_corners(xp,yp,zp,bp)
    
        
print('Ellipsoid carbide size and location generated')

#creating model
#Mdb() if creating new file
model = mdb.Model(name='Model-1',modelType = STANDARD_EXPLICIT)

#creating materials
model.Material(name = 'Carbide')
model.materials['Carbide'].Elastic(table=((0.365707, 0.25), ))
model.Material(name='Matrix')
model.materials['Matrix'].Elastic(table=((0.266272, 0.3), ))
model.materials['Matrix'].Plastic(table=((
        0.00045, 0.0), (0.000512664760145, 0.000438040461157546), (
        0.000556893977682, 0.00109572525679919), (0.000600686982051, 
        0.00213631117327229), (0.000643655185527, 0.00367752154598405), (
        0.000685256060044, 0.00580997348841567), (0.000724903537752, 
        0.00856869796541231), (0.000762153777714, 0.011921407262596), (
        0.000796764037128, 0.0157990158182437), (0.000828810123759, 
        0.02009939170478), (0.000858444125988, 0.0247392113866344), (
        0.000886026526206, 0.029631785363798), (0.000911823317349, 
        0.034721658924612), (0.000936106869976, 0.0399638076404829), (
        0.000959109337081, 0.0453250684560245), (0.000981026313674, 
        0.0507806594994877), (0.00100202419481, 0.0563116251807552), (
        0.00102223787629, 0.0619037826486917), (0.00104177730817, 
        0.06754633150446), (0.00106073585037, 0.0732306827182445), (
        0.0010792090724, 0.0789488418140858), (0.00109722549955, 
        0.0846979309706963), (0.00111486253249, 0.0904720226851532), (
        0.00113218014588, 0.0962667025476286), (0.0011491883141, 
        0.102080947038558), (0.00116594213948, 0.107910965437832), (
        0.00118247858244, 0.113754258310442), (0.00119880426498, 
        0.119610255359081), (0.00121496604517, 0.125475997255539), (
        0.00123096909518, 0.13135108794343), (0.00124683042637, 
        0.137234437042601), (0.00126257864394, 0.143124317987559), (
        0.00127821293239, 0.14902073467496), (0.00129376530889, 
        0.154921809160459), (0.00130923671821, 0.160827472434916), (
        0.00132463222402, 0.166737415090426), (0.00133997701206, 
        0.172650196366473), (0.00135527169974, 0.17856577873434), (
        0.00137051890012, 0.184484009556526), (0.00138573775067, 
        0.19040381484515), (0.001393337, 0.193364281767816)))
print('Material properties created')

#creating sections

model.HomogeneousSolidSection(name='CarbideSection',material = 'Carbide',thickness=None)
model.HomogeneousSolidSection(name='MatrixSection',material = 'Matrix',thickness=None)

print('Section created')

#creating matrix cube

s = model.ConstrainedSketch(name='__profile__',sheetSize = 200.0)
g,v,d,c = s.geometry, s.vertices, s.dimensions, s.constraints
s.setPrimaryObject(option=STANDALONE)
s.rectangle(point1 = (-0.5*cubeEdge,-0.5*cubeEdge), point2 = (0.5*cubeEdge,0.5*cubeEdge))
model.Part(name='CubePart',dimensionality = THREE_D, type=DEFORMABLE_BODY)
cubePart = model.parts['CubePart']
cubePart.BaseSolidExtrude(sketch=s,depth=cubeEdge)
s.unsetPrimaryObject()
del model.sketches['__profile__']

#assigning section

matrixCell = cubePart.cells.findAt(((0.5*cubeEdge,0.5*cubeEdge,0.5*cubeEdge),))
cubePart.Set(cells=matrixCell,name='Set-Matrix')
region = cubePart.sets['Set-Matrix']
cubePart.SectionAssignment(region=region,sectionName='MatrixSection',offset=0.0,offsetType=MIDDLE_SURFACE,offsetField='',thicknessAssignment=FROM_SECTION)

print('Matrix section assigned')

#creating cube assembly

assembly = model.rootAssembly
assembly.DatumCsysByDefault(CARTESIAN)
assembly.Instance(name='CubeMatrix',part = cubePart, dependent = ON)
assembly.translate(instanceList = ('CubeMatrix',),vector = (0.0,0.0,-0.5*cubeEdge))

#creating carbide ellipsoid

i = 0
while i<noOfCarbides:
    s = mdb.models['Model-1'].ConstrainedSketch(name='__profile__',sheetSize = 5.0)
    g,v,d,c = s.geometry, s.vertices,s.dimensions,s.constraints
    constructionLine = s.ConstructionLine(point1=(0.0,-100.0),point2 = (0.0,100.0))
    constructionLineID = constructionLine.id
    s.FixedConstraint(entity=g[constructionLineID])
    e = s.EllipseByCenterPerimeter((0,0),(b[i],0),(0,a[i]))
    eid = e.id
    l = s.Line((0,a[i]),(0,-a[i]))
    lid = l.id
    s.autoTrimCurve(g[eid],(-b[i],0))
    p = model.Part(name = 'Carbide-{}'.format(i+1),dimensionality = THREE_D,type = DEFORMABLE_BODY)
    p.BaseSolidRevolve(sketch = s,angle = 360.0,flipRevolveDirection = OFF)
    p = model.parts['Carbide-{}'.format(i+1)]
    del model.sketches['__profile__']
    carbideCell = p.cells.findAt(((0.0,0.0,0.0),))
    p.Set(cells=carbideCell,name='Set-Carbide')
    region = p.sets['Set-Carbide']
    p.SectionAssignment(region=region,sectionName='CarbideSection',offset=0.0,offsetType=MIDDLE_SURFACE,offsetField='',thicknessAssignment=FROM_SECTION)
    XYP = p.DatumPlaneByPrincipalPlane(principalPlane = XYPLANE, offset = 0.0)
    XYPid = XYP.id
    XZP = p.DatumPlaneByPrincipalPlane(principalPlane = XZPLANE, offset = 0.0)
    XZPid = XZP.id
    YZP = p.DatumPlaneByPrincipalPlane(principalPlane = YZPLANE, offset = 0.0)
    YZPid = YZP.id
    d = p.datums
    c = p.cells.findAt((0.0,0.0,0.0),)
    p.PartitionCellByDatumPlane(datumPlane = d[XYPid],cells = c)
    c = p.cells.getByBoundingBox(-b[i], -a[i], -b[i], b[i], a[i], b[i])
    p.PartitionCellByDatumPlane(datumPlane=d[XZPid], cells=c)
    c = p.cells.getByBoundingBox(-b[i], -a[i], -b[i], b[i], a[i], b[i])
    p.PartitionCellByDatumPlane(datumPlane=d[YZPid], cells=c)
    print('Created carbide ',i+1,' of ',noOfCarbides)
    #creating assembly
    assembly.Instance(name = 'Carbide-{}'.format(i+1),part = p, dependent = ON)
    assembly.translate(instanceList = ('Carbide-{}'.format(i+1),),vector = (x[i],y[i],z[i]))
    if i == 0:
        selectedInstances = (assembly.instances['Carbide-{}'.format(i+1)],)
    else:
        carbideInstances = assembly.instances['Carbide-{}'.format(i+1)]
        selectedInstances = selectedInstances + (carbideInstances,)
    #increment
    i = i + 1
#loop over
selectedInstances = selectedInstances + (assembly.instances['CubeMatrix'],)
assembly.InstanceFromBooleanMerge(name = 'MicroMarM', instances = selectedInstances, 
                                 keepIntersections = ON,originalInstances=DELETE, domain=GEOMETRY)

print('carbide and matrix are merged\nMicro MarM created')
assembly.regenerate()

i = 0
while i <noOfCarbides:
    del model.parts['Carbide-{}'.format(i+1)]
    i += 1

#creating step
model.StaticStep(name='Step-1', previous='Initial',maxNumInc=1000, initialInc=0.001, minInc=1e-08, maxInc=0.1)
print('Created step-1')

#creating boundary conditions
#encastreBC

marmFace = assembly.instances['MicroMarM-1'].faces
encastreFace = marmFace.findAt(((0.0,-0.5*cubeEdge,0.0),))
region = assembly.Set(faces=encastreFace,name='Set-Encastre')
model.EncastreBC(name='EncastreBC',createStepName='Initial',region=region,localCsys=None)

#displacementBC
displacementFace = marmFace.findAt(((0.0,0.5*cubeEdge,0.0),))
region = assembly.Set(faces=displacementFace,name='Set-Displacement')
model.DisplacementBC(name='DisplacementBC',createStepName='Step-1',region=region,u1=0.0, u2=2.0, u3=0.0, ur1=0.0, ur2=0.0, ur3=0.0, 
        amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', 
        localCsys=None)
        
print('BC created')

