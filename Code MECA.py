from abaqus import *
from abaqus import backwardCompatibility
backwardCompatibility.setValues(reportDeprecated=False)

from abaqusConstants import *
import __main__

import math
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

############## VARIABLES #################

# IMPORT USER VARIABLES #
with open("constant.py", "r") as f:
    exec(f.read())                             # Maximal temperature reachable per increment (K)

# CONSTANTS #
ABSOLUT_ZERO = 0                               # Absolute zero (K)
STEFANBOLTZMANN = 5.67e-8                      # Stefan–Boltzmann constant (W/m2.K4)
R = 8.314                                      # Ideal gas constant

# MESH #
LAYERMESH_SIZE = 15e-06                        # Size of the largest layer mesh (the smallest mesh size is LAYERMESH_SIZE/3)
SUPPORTMESH_SIZE = 0.8e-04                     # Size of the largest support mesh (the smallest mesh size is LAYERMESH_SIZE to have a continuity btw layer&support) (m)

# SIMULATION #
NB_INCREMENT = 100000                           # Max number of increments per step
MIN_INCR = 1e-20                                # Minimal time increment size (s)
INITIAL_INCR = H_SPACING/(LASER_SPEED*5)        # Initial time increment size (s)
MAXTEMP_INCR = 800                              # Maximal temperature reachable per increment (K)

# JOB #
TEMP_FILE = JOB                                 # "D:/Eglantine/Abaqus/"  Location and name of the Thermal job (without ".odb" at the end)
NB_STEPS = (2*NB_LASER)*NB_LAYERS+1 # Nb total of steps (used to define the end of input temperature field)


############## ALGORITHMES ###############

def step_stat(filename, nb_step):
# Code to extract the number of increment per step of the thermal Job & send it to the strain Job ##

    results = [None] * nb_step
    last_step = None
    last_increment = None

    with open(filename, 'r') as f:
        for line in f: # Read line-by-line
            tokens = line.strip().split() # Split line when space
            if len(tokens) < 2:
                continue
            try:
                step = int(tokens[0])
                increment = int(tokens[1])
            except ValueError:
                continue

            if (last_step is not None and step != last_step):
                # if the step change from one line to another, we save the first line data (last increment of the step)
                results[last_step-1] = last_increment #"-1" bc table index starts at 0

            last_step = step
            last_increment = increment

    # for the last step
    if last_step is not None:
        results[last_step-1] = last_increment

    return results

sta_stat = step_stat(TEMP_FILE+".sta",NB_STEPS)      # Table of nb of increments per steps on the Thermal Job extracted from the sta file


def find_top_face(s1):
# Trouver la face superieure de chaque cellule
    tab = []
    for face in s1.faces:
        normal = face.getNormal()
        if normal == (0,1,0) :
            tab.append(face)
    return tab

def find_lower_face(s1):
# Trouver la face inférieure de chaque cellule
    tab = []
    for face in s1.faces:
        normal = face.getNormal()
        if normal == (0,-1,0) :
            tab.append(face)
    return tab

def find_contour_face(s1):
# Trouver les faces extérieures de chaque cellule
    tab_x = []
    tab_z = []
    for face in s1.faces:
        if len(face.getCells()) == 1: # Exterior face
            normal = face.getNormal()
            if normal == (1,0,0) or normal == (-1,0,0) : #x
                tab_x.append(face)
            if normal == (0,0,1) or normal == (0,0,-1) : #z
                tab_z.append(face)
    return tab_x, tab_z


            
###########################################

def Constant():
# Definir les constantes universelles du modele
    mdb.models['Model-1'].setValues(absoluteZero=ABSOLUT_ZERO, stefanBoltzmann=STEFANBOLTZMANN, universalGas=R)

def Material():
# Creer les materiaux du support et layer
    ## Stainless Steel 316L Solid ##
    mdb.models['Model-1'].Material(name=MAT_NAME+'_Solid') 
    mdb.models['Model-1'].materials[MAT_NAME+'_Solid'].Density(table = ((ambient_density, ),) )
    mdb.models['Model-1'].materials[MAT_NAME+'_Solid'].Conductivity(temperatureDependency=ON, table = data_cond)
    mdb.models['Model-1'].materials[MAT_NAME+'_Solid'].SpecificHeat(temperatureDependency=ON, table = data_spe_heat)
    mdb.models['Model-1'].materials[MAT_NAME+'_Solid'].LatentHeat(table=((LATENT_HEAT, SOLIDUS, LIQUIDUS), ))
    mdb.models['Model-1'].materials[MAT_NAME+'_Solid'].Elastic(temperatureDependency=ON, table= data_elastic)
    mdb.models['Model-1'].materials[MAT_NAME+'_Solid'].Expansion(temperatureDependency=ON, table=data_expansion)
    mdb.models['Model-1'].materials[MAT_NAME+'_Solid'].Plastic(temperatureDependency=ON, scaleStress=None, table=data_plastic)
    mdb.models['Model-1'].materials[MAT_NAME+'_Solid'].plastic.AnnealTemperature(table=((SOLIDUS, ), )) # Once temp reach solidus, all plastic strain are forgotten has the solid melt


def Section():
# Creer les sections 'Support' et 'Layer' pour ensuite assigner les materiaux
    mdb.models['Model-1'].HomogeneousSolidSection(name='Support', material=MAT_NAME+'_Solid', thickness=None)
    mdb.models['Model-1'].HomogeneousSolidSection(name='Layer', material=MAT_NAME+'_Solid', thickness=None)


def SupportCreate():
# Creation du support & mesh

    ## Creer le support ##
    s1 = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=0.01)
    s1.sketchOptions.setValues(decimalPlaces=4)
    s1.rectangle(point1=(0, -SUPPORT_THICKNESS), point2=(LAYER_LENGTH, 0))
    p = mdb.models['Model-1'].Part(name='Support', dimensionality=THREE_D, type=DEFORMABLE_BODY)
    p.BaseSolidExtrude(sketch=s1, depth=LAYER_WIDTH)
    s1.unsetPrimaryObject()
    del mdb.models['Model-1'].sketches['__profile__']

    ## Ajouter le materiau au support ##
    # ! Avoir creer la section 'Support' avant
    region = regionToolset.Region(cells=p.cells[:])
    p.SectionAssignment(region=region, sectionName='Support', thicknessAssignment=FROM_SECTION)

    ## Creer le mesh du support ##
    p = mdb.models['Model-1'].parts['Support']
    for edge in p.edges:
        # Tangent to the edge
        verts = edge.getVertices()
        v1 = p.vertices[verts[0]].pointOn # Coord extremité 1 edge
        v2 = p.vertices[verts[1]].pointOn # Coord extremité 2 edge
        p1 = v1[0]
        p2 = v2[0]
        tx = p2[0] - p1[0]
        ty = p2[1] - p1[1]
        tz = p2[2] - p1[2]

        if(tx == 0 and tz == 0): 
            # Si orienté selon y, on retrécit le mesh à la surface
            if(ty > 0):
                # Si edge orienté vers le haut
                p.seedEdgeByBias(biasMethod=SINGLE, end2Edges=part.EdgeArray((edge,)), minSize=LAYERMESH_SIZE, maxSize=SUPPORTMESH_SIZE, constraint=FINER)
            else :
                # Si edge orienté vers le bas
                p.seedEdgeByBias(biasMethod=SINGLE, end1Edges=part.EdgeArray((edge,)), minSize=LAYERMESH_SIZE, maxSize=SUPPORTMESH_SIZE, constraint=FINER)
        else :
            # Sinon mesh normal
            p.seedEdgeBySize(edges=part.EdgeArray((edge,)), size=1.5*LAYERMESH_SIZE, deviationFactor=0.1, constraint=FINER)
    p.generateMesh()


def LayerCreate1():
# Creation d'une couche et des cells

    for layer in range(NB_LAYERS):
        ## Creation de la couche ##
        s1 = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=0.01)
        s1.sketchOptions.setValues(decimalPlaces=4)
        s1.rectangle(point1=(0, layer*LAYER_THICKNESS), point2=(LAYER_LENGTH, LAYER_THICKNESS+layer*LAYER_THICKNESS))
        p = mdb.models['Model-1'].Part(name=f'Layer-{layer}', dimensionality=THREE_D, type=DEFORMABLE_BODY)
        p.BaseSolidExtrude(sketch=s1, depth=LAYER_WIDTH)

        s1.unsetPrimaryObject()
        del mdb.models['Model-1'].sketches['__profile__']

        ## Ajouter le materiau a la couche ##
        # ! Avoir creer la section 'Layer' avant
        region = regionToolset.Region(cells=p.cells[:])
        p.SectionAssignment(region=region, sectionName='Layer', thicknessAssignment=FROM_SECTION)

        ## Decoupage en plans ##    
        ref_face = p.faces

        for i in range(1, NB_LASER):
            offset_val = i * H_SPACING
            # Decoupage x 
            p.DatumPlaneByOffset(plane=ref_face[4], flip=SIDE2, offset=offset_val)
            
        ## Decoupage des plans en cells ##
        datum_planes = [d for d in p.datums.values()]
        for plane in datum_planes:
            pickedCells = p.cells[:]
            p.PartitionCellByDatumPlane(datumPlane=plane, cells=pickedCells)

        ## Mesh ##
        for edge in p.edges:
            # Tangent to the edge
            verts = edge.getVertices()
            v1 = p.vertices[verts[0]].pointOn # Coord extremité 1 edge
            v2 = p.vertices[verts[1]].pointOn # Coord extremité 2 edge
            p1 = v1[0]
            p2 = v2[0]
            tx = p2[0] - p1[0]
            ty = p2[1] - p1[1]
            tz = p2[2] - p1[2]

            if(tx == 0 and tz == 0): 
                # Si orienté selon y, on retrécit le mesh à la surface
                if(ty > 0):
                    # Si edge orienté vers le haut
                    p.seedEdgeByBias(biasMethod=SINGLE, end2Edges=part.EdgeArray((edge,)), minSize=LAYERMESH_SIZE/3, maxSize=LAYERMESH_SIZE, constraint=FINER)
                else :
                    # Si edge orienté vers le bas
                    p.seedEdgeByBias(biasMethod=SINGLE, end1Edges=part.EdgeArray((edge,)), minSize=LAYERMESH_SIZE/3, maxSize=LAYERMESH_SIZE, constraint=FINER)
            elif(tx == 0 and ty == 0):
                # Si orienté selon tz, on reduit le mesh au milieu
                if(tz > 0):
                    p.seedEdgeByBias(biasMethod=DOUBLE, centerEdges=part.EdgeArray((edge,)), minSize=LAYERMESH_SIZE/2, maxSize=LAYERMESH_SIZE, constraint=FINER)
                else :
                    p.seedEdgeByBias(biasMethod=DOUBLE, endEdges=part.EdgeArray((edge,)), minSize=LAYERMESH_SIZE/2, maxSize=LAYERMESH_SIZE, constraint=FINER)
            else :
                # Sinon mesh normal
                p.seedEdgeBySize(edges=part.EdgeArray((edge,)), size=LAYERMESH_SIZE, deviationFactor=0.1, constraint=FINER)
        p.generateMesh()

        # Creation d'un set 'LayerSet' qui contient toute la couche
        p.Set(elements=p.elements, name=f'LayerSet-{layer}')
        p.Set(nodes=p.nodes, name=f'LayerNodes-{layer}')


def Assembly2():
    # Assembler la couche et le support
    a = mdb.models['Model-1'].rootAssembly
    a.DatumCsysByDefault(CARTESIAN)
    allsets = [] #Pour créer un set qui contient tout

    ## Assembler les layers
    for layer in range(NB_LAYERS):
        p = mdb.models['Model-1'].parts[f'Layer-{layer}']
        inst = a.Instance(name=f'Layer-{layer}', part=p, dependent=ON)
        inst_set = inst.sets[f'LayerSet-{layer}']
        allsets.append(inst_set)

    a.SetByMerge(name='AllLayers', sets=allsets)

    ## Assembler le support
    p = mdb.models['Model-1'].parts['Support']
    a.Instance(name='Support-1', part=p, dependent=ON)
    
def Steps3():
# Creer les activations de cells a chaque etape
    a = mdb.models['Model-1'].rootAssembly
    a.regenerate()
    
    ## Initial : Toutes cells eteintes ##
    mdb.models['Model-1'].StaticStep(name='Death', previous='Initial', timePeriod=0.1, initialInc=0.01, minInc=1e-8, maxInc=0.1, nlgeom=ON)
    region = a.sets['AllLayers']
    mdb.models['Model-1'].ModelChange(name='Desactive', createStepName='Death', region=region, regionType=ELEMENTS, activeInStep=False, includeStrain=False) # Desactiver toutes les cells au debut
    
    # Temperature Field #
    mdb.models['Model-1'].Temperature(name='Temperature-Death', createStepName='Death', distributionType=FROM_FILE, fileName=TEMP_FILE+".odb", beginStep=1, beginIncrement=1, endStep=1, endIncrement=sta_stat[0], interpolate=OFF, absoluteExteriorTolerance=0.0, exteriorTolerance=0.05)
    

    for layer in range(NB_LAYERS):
        ## Etape 1 : Active la premiere cell ##
        if layer == 0 :
            mdb.models['Model-1'].StaticStep(name='Layer-0', previous='Death', timePeriod=(LAYER_TIME*(NB_LASER-1))+LASER_TIME+COOLING_TIME, maxNumInc=NB_INCREMENT*NB_LASER, initialInc=INITIAL_INCR, minInc=MIN_INCR, maxInc=LASER_TIME, nlgeom=ON)
        elif layer == NB_LAYERS-1 : # Longer cooling time at the end
            mdb.models['Model-1'].StaticStep(name=f'Layer-{layer}', previous=f'Layer-{layer-1}', timePeriod=(LAYER_TIME*(NB_LASER-1))+LASER_TIME+COOLING_TIME*10, maxNumInc=NB_INCREMENT*NB_LASER, initialInc=INITIAL_INCR, minInc=MIN_INCR, maxInc=LASER_TIME, nlgeom=ON)
        else : 
            mdb.models['Model-1'].StaticStep(name=f'Layer-{layer}', previous=f'Layer-{layer-1}', timePeriod=(LAYER_TIME*(NB_LASER-1))+LASER_TIME+COOLING_TIME, maxNumInc=NB_INCREMENT*NB_LASER, initialInc=INITIAL_INCR, minInc=MIN_INCR, maxInc=LASER_TIME, nlgeom=ON)
        
        region =a.instances[f'Layer-{layer}'].sets[f'LayerSet-{layer}']
        mdb.models['Model-1'].ModelChange(name=f'Activate-{layer}', createStepName=f'Layer-{layer}', region=region, regionType=ELEMENTS, activeInStep=True, includeStrain=False)
        
        # Temperature Field #
        mdb.models['Model-1'].Temperature(name=f'Temperature-{layer}', createStepName=f'Layer-{layer}', distributionType=FROM_FILE, fileName=TEMP_FILE+".odb", beginStep=(2*NB_LASER)*layer+2, beginIncrement=1, endStep=(2*NB_LASER)*(layer+1)+1, endIncrement=sta_stat[(2*NB_LASER)*(layer+1)], interpolate=OFF, absoluteExteriorTolerance=0.0, exteriorTolerance=0.05)


def Surface4():
# Relier les couches entre elles
    a = mdb.models['Model-1'].rootAssembly
    mdb.models['Model-1'].ContactProperty('Surf2Surf')
    mdb.models['Model-1'].interactionProperties['Surf2Surf'].TangentialBehavior(formulation=FRICTIONLESS)
    mdb.models['Model-1'].interactionProperties['Surf2Surf'].NormalBehavior(pressureOverclosure=HARD, allowSeparation=OFF, constraintEnforcementMethod=DEFAULT)
    
    ## Lier face support et face layer-0 ##
    s1 = a.instances['Support-1']
    face1 = find_top_face(s1)
    region1=regionToolset.Region(side1Faces=part.FaceArray(face1))
    s2 = a.instances['Layer-0']
    face2 = find_lower_face(s2)
    region2=regionToolset.Region(side1Faces=part.FaceArray(face2))
    regionDef = s2.sets['LayerNodes-0']
    #mdb.models['Model-1'].SurfaceToSurfaceContactStd(name='Surface-Support', createStepName='Layer-0', main=region1, secondary=region2,  sliding=FINITE, thickness=ON, interactionProperty='Surf2Surf', adjustMethod=OVERCLOSED, initialClearance=OMIT, datumAxis=None,  clearanceRegion=None, tied=OFF)
    mdb.models['Model-1'].Tie(name='Surface-Support', main=region1, secondary=region2, positionToleranceMethod=COMPUTED, adjust=ON, tieRotations=ON, thickness=ON)

    ## Lier les autres layers ##
    for layer in range(NB_LAYERS-1):
        s1 = a.instances[f'Layer-{layer}']
        face1 = find_top_face(s1)
        region1 = a.Surface(side1Faces=part.FaceArray(face1), name=f'TopFace-Layer{layer}')
        s2 = a.instances[f'Layer-{layer+1}']
        face2 = find_lower_face(s2)
        region2=regionToolset.Region(side1Faces=part.FaceArray(face2))
        regionDef = s2.sets[f'LayerNodes-{layer+1}']
        mdb.models['Model-1'].SurfaceToSurfaceContactStd(name=f'Layer{layer}-Layer{layer+1}', createStepName=f'Layer-{layer+1}', main=region1, secondary=region2,  sliding=FINITE, thickness=ON, interactionProperty='Surf2Surf', adjustMethod=OVERCLOSED, initialClearance=OMIT, datumAxis=None,  clearanceRegion=None, tied=OFF)
        #mdb.models['Model-1'].Tie(name=f'Layer{layer}-Layer{layer+1}', main=region1, secondary=region2, positionToleranceMethod=COMPUTED, adjust=ON, tieRotations=ON, thickness=ON)


def Loads6():
# Definir les champs de temperature et gravité 

    ## Gravity ##
    mdb.models['Model-1'].Gravity(name='Gravity', createStepName='Death', comp2=-9.81, distributionType=UNIFORM, field='')
    
    ## Symmetry on the contours ##
    a = mdb.models['Model-1'].rootAssembly
    s = a.instances['Support-1']
    x_faces, z_faces = find_contour_face(s)
    y_face = find_lower_face(s)
    for layer in range(NB_LAYERS):
        s = a.instances[f'Layer-{layer}']
        x_f, z_f = find_contour_face(s)
        x_faces += x_f
        z_faces += z_f

    region_x = regionToolset.Region(faces=part.FaceArray(x_faces))
    region_z = regionToolset.Region(faces=part.FaceArray(z_faces))
    region_y =regionToolset.Region(faces=part.FaceArray(y_face))
    mdb.models['Model-1'].ZsymmBC(name='Symmetry_z', createStepName='Initial', region=region_z, localCsys=None)
    mdb.models['Model-1'].XsymmBC(name='Symmetry_x', createStepName='Initial', region=region_x, localCsys=None)
    mdb.models['Model-1'].DisplacementBC(name='NoDisp_y', createStepName='Initial', region=region_y, u1=UNSET, u2=SET, u3=UNSET, ur1=UNSET, ur2=UNSET, ur3=UNSET, amplitude=UNSET, distributionType=UNIFORM, fieldName='', localCsys=None)


def LayerElementType8():
# Definir le format de notre simulation pour une couche (ici static, mechanical)
    elemType1 = mesh.ElemType(elemCode=C3D8R, elemLibrary=STANDARD, kinematicSplit=AVERAGE_STRAIN, secondOrderAccuracy=OFF, hourglassControl = DEFAULT, distortionControl=DEFAULT)
    elemType2 = mesh.ElemType(elemCode=C3D6, elemLibrary=STANDARD)
    elemType3 = mesh.ElemType(elemCode=C3D4, elemLibrary=STANDARD)
    for layer in range(NB_LAYERS):
        p = mdb.models['Model-1'].parts[f'Layer-{layer}']
        cells = p.cells[:]
        p.setElementType(regions=(cells,), elemTypes=(elemType1, elemType2, elemType3))


def SupportElementType9():
# Definir le format de notre simulation (ici static, mechanical)
    elemType1 = mesh.ElemType(elemCode=C3D8R, elemLibrary=STANDARD, kinematicSplit=AVERAGE_STRAIN, secondOrderAccuracy=OFF, hourglassControl=DEFAULT, distortionControl=DEFAULT)
    elemType2 = mesh.ElemType(elemCode=C3D6, elemLibrary=STANDARD)
    elemType3 = mesh.ElemType(elemCode=C3D4, elemLibrary=STANDARD)
    p = mdb.models['Model-1'].parts['Support']
    cells = p.cells[:]
    p.setElementType(regions=(cells,), elemTypes=(elemType1, elemType2, elemType3))


def FieldOutput10():
# Define Field Outputs & Create Job
    a = mdb.models['Model-1']

    # Field Outputs #
    for req in list(a.fieldOutputRequests.keys()): # delete all previous outputs
        del a.fieldOutputRequests[req]
    a.FieldOutputRequest(name='Thermal-Strain', createStepName='Death', variables=('S', 'THE', 'U', 'ENER', 'ELEDEN', 'NT'))

    # Job #
    a.rootAssembly.regenerate()
    mdb.Job(name=JOB+'-Mech', model='Model-1', description='', type=ANALYSIS, 
        atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
        memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
        explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
        modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
        scratch='', resultsFormat=ODB, numThreadsPerMpiProcess=1, 
        multiprocessingMode=DEFAULT, numCpus=NB_CPU, numDomains=NB_CPU, numGPUs=0)
    del mdb.models['Model-1'].historyOutputRequests['H-Output-1']

    

Constant()
Section()
Material()
SupportCreate()
LayerCreate1()
Assembly2()
Steps3()
Surface4()
Loads6()
LayerElementType8()
SupportElementType9()
FieldOutput10()