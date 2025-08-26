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
    exec(f.read())

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
MAXTEMP_INCR = 800                              # Maximal temperature change per increment (K)

############## ALGORITHMES ###############

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

###########################################

def Constant():
# Definir les constantes universelles du modele
    mdb.models['Model-1'].setValues(absoluteZero=ABSOLUT_ZERO, stefanBoltzmann=STEFANBOLTZMANN, universalGas=R)

def Material():
# Creer les materiaux du support et layer
    ## Stainless Steel 316L Solid ##
    mdb.models['Model-1'].Material(name=MAT_NAME+'_Solid') 
    mdb.models['Model-1'].materials[MAT_NAME+'_Solid'].Density(temperatureDependency=ON, table = data_density)
    mdb.models['Model-1'].materials[MAT_NAME+'_Solid'].Conductivity(temperatureDependency=ON, table = data_cond)
    mdb.models['Model-1'].materials[MAT_NAME+'_Solid'].SpecificHeat(temperatureDependency=ON, table = data_spe_heat)
    mdb.models['Model-1'].materials[MAT_NAME+'_Solid'].LatentHeat(table=((LATENT_HEAT, SOLIDUS, LIQUIDUS), ))
    mdb.models['Model-1'].materials[MAT_NAME+'_Solid'].Elastic(temperatureDependency=ON, table= data_elastic)
    mdb.models['Model-1'].materials[MAT_NAME+'_Solid'].Expansion(temperatureDependency=ON, table=data_expansion)
    
    ## Stainless Steel 316L Powder ##
    mdb.models['Model-1'].Material(name=MAT_NAME+'_Powder') 
    mdb.models['Model-1'].materials[MAT_NAME+'_Powder'].Depvar(n=1) # User-defined material variable
    mdb.models['Model-1'].materials[MAT_NAME+'_Powder'].UserDefinedField() # Equal to Depvar, used for Density dependency on state
    mdb.models['Model-1'].materials[MAT_NAME+'_Powder'].UserMaterial(type=THERMAL, thermalConstants=(LATENT_HEAT, ))
    mdb.models['Model-1'].materials[MAT_NAME+'_Powder'].Density(temperatureDependency=ON, dependencies=1, table = data_density_powder)
    mdb.models['Model-1'].materials[MAT_NAME+'_Powder'].Elastic(temperatureDependency=ON, table=data_elastic)
    mdb.models['Model-1'].materials[MAT_NAME+'_Powder'].Expansion(temperatureDependency=ON, table=data_expansion)


def Section():
# Creer les sections 'Support' et 'Layer' pour ensuite assigner les materiaux
    mdb.models['Model-1'].HomogeneousSolidSection(name='Support', material=MAT_NAME+'_Solid', thickness=None)
    mdb.models['Model-1'].HomogeneousSolidSection(name='Layer', material=MAT_NAME+'_Powder', thickness=None)


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
    mdb.models['Model-1'].HeatTransferStep(name='Death', previous='Initial', timePeriod=0.1, initialInc=0.01, minInc=1e-8, maxInc=0.1, deltmx=350.0) # Creer l'etape initiale
    region = a.sets['AllLayers']
    mdb.models['Model-1'].ModelChange(name='Desactive', createStepName='Death', region=region, regionType=ELEMENTS, activeInStep=False, includeStrain=False) # Desactiver toutes les cells au debut

    for layer in range(NB_LAYERS):
        ## Etape 1 : Active la premiere cell ##
        if layer == 0 :
            mdb.models['Model-1'].HeatTransferStep(name='Layer-0Step-0', previous='Death', timePeriod=LASER_TIME, maxNumInc=NB_INCREMENT, initialInc=INITIAL_INCR, minInc=MIN_INCR, maxInc=5*INITIAL_INCR, deltmx=MAXTEMP_INCR)
            region =a.instances['Layer-0'].sets['LayerSet-0']
            mdb.models['Model-1'].ModelChange(name='Activate-0', createStepName='Layer-0Step-0', region=region, regionType=ELEMENTS, activeInStep=True, includeStrain=False)
        else :
            mdb.models['Model-1'].HeatTransferStep(name=f'Layer-{layer}Step-0', previous=f'Cooling-{layer-1}', timePeriod=LASER_TIME, maxNumInc=NB_INCREMENT, initialInc=INITIAL_INCR, minInc=MIN_INCR, maxInc=5*INITIAL_INCR, deltmx=MAXTEMP_INCR)
            region =a.instances[f'Layer-{layer}'].sets[f'LayerSet-{layer}']
            mdb.models['Model-1'].ModelChange(name=f'Activate-{layer}', createStepName=f'Layer-{layer}Step-0', region=region, regionType=ELEMENTS, activeInStep=True, includeStrain=False)
            
        for i in range(1, NB_LASER):
            mdb.models['Model-1'].HeatTransferStep(name=f'Cooling-{layer}Step-{i-1}', previous=f'Layer-{layer}Step-{i-1}', timePeriod=LAYER_TIME, maxNumInc=NB_INCREMENT, initialInc=INITIAL_INCR, minInc=MIN_INCR, maxInc=LAYER_TIME/2, deltmx=MAXTEMP_INCR)
            mdb.models['Model-1'].HeatTransferStep(name=f'Layer-{layer}Step-{i}', previous=f'Cooling-{layer}Step-{i-1}', timePeriod=LASER_TIME, maxNumInc=NB_INCREMENT, initialInc=INITIAL_INCR, minInc=MIN_INCR, maxInc=5*INITIAL_INCR, deltmx=MAXTEMP_INCR)
            
        ## Cooling entre chaque couche
        if layer == NB_LAYERS-1 : # Longer final cooling
            mdb.models['Model-1'].HeatTransferStep(name=f'Final-Cooling', previous=f'Layer-{NB_LAYERS-1}Step-{NB_LASER-1}', timePeriod=100*COOLING_TIME, maxNumInc=NB_INCREMENT, initialInc=INITIAL_INCR, minInc=MIN_INCR, maxInc=5*COOLING_TIME, deltmx=MAXTEMP_INCR)
        else :
            mdb.models['Model-1'].HeatTransferStep(name=f'Cooling-{layer}', previous=f'Layer-{layer}Step-{NB_LASER-1}', timePeriod=COOLING_TIME, maxNumInc=NB_INCREMENT, initialInc=INITIAL_INCR, minInc=MIN_INCR, maxInc=COOLING_TIME/2, deltmx=MAXTEMP_INCR)
    

def Surface4():
# Relier les couches entre elles
    a = mdb.models['Model-1'].rootAssembly
    
    ## Lier face support et face layer-0 ##
    s1 = a.instances['Support-1']
    face1 = find_top_face(s1)
    region1=regionToolset.Region(side1Faces=part.FaceArray(face1))
    s2 = a.instances['Layer-0']
    face2 = find_lower_face(s2)
    region2=regionToolset.Region(side1Faces=part.FaceArray(face2))
    mdb.models['Model-1'].Tie(name='Surface-Support', main=region1, secondary=region2, positionToleranceMethod=COMPUTED, adjust=ON, tieRotations=ON, thickness=ON)

    ## Lier les autres layers ##
    for layer in range(NB_LAYERS-1):
        s1 = a.instances[f'Layer-{layer}']
        face1 = find_top_face(s1)
        region1 = a.Surface(side1Faces=part.FaceArray(face1), name=f'TopFace-Layer{layer}')
        s2 = a.instances[f'Layer-{layer+1}']
        face2 = find_lower_face(s2)
        region2=regionToolset.Region(side1Faces=part.FaceArray(face2))
        mdb.models['Model-1'].Tie(name=f'Layer{layer}-Layer{layer+1}', main=region1, secondary=region2, positionToleranceMethod=COMPUTED, adjust=ON, tieRotations=ON, thickness=ON)

    ## Créer la topsurface pour la dernière couche ##
    s1 = a.instances[f'Layer-{NB_LAYERS-1}']
    face1 = find_top_face(s1)
    a.Surface(side1Faces=part.FaceArray(face1), name=f'TopFace-Layer{NB_LAYERS-1}')


def HeatTransfer5():
# Definir les differentes rayonnements du support & couches 
    a = mdb.models['Model-1'].rootAssembly
    c1 = a.instances['Support-1'].cells[:]

    ## LAYERS ##
    for layer in range(NB_LAYERS):
        regiontop = a.surfaces[f'TopFace-Layer{layer}']

        # laser  top #
        for i in range(NB_LASER):   
            mdb.models['Model-1'].SurfaceHeatFlux(name=f'Laser-{layer}{i}', createStepName=f'Layer-{layer}Step-{i}', region=regiontop, magnitude=1.0, distributionType=USER_DEFINED)
            if i == NB_LASER-1:
                if(layer < NB_LAYERS-1):
                    mdb.models['Model-1'].loads[f'Laser-{layer}{i}'].deactivate(f'Cooling-{layer}')
            else :
                mdb.models['Model-1'].loads[f'Laser-{layer}{i}'].deactivate(f'Layer-{layer}Step-{i+1}')
        
        # Radiative & Convective top  #
        mdb.models['Model-1'].RadiationToAmbient(ambientTemperature=AMBIENT_TEMP, ambientTemperatureAmp='', createStepName=f'Layer-{layer}Step-0', distributionType=UNIFORM, emissivity=EMISSIVITY, field='', name=f'Radiation-{layer}', radiationType=AMBIENT, surface=regiontop)
        mdb.models['Model-1'].FilmCondition(name=f'Convection-{layer}', createStepName=f'Layer-{layer}Step-0', surface=regiontop, definition=EMBEDDED_COEFF, filmCoeff=FILMCOEFF, filmCoeffAmplitude='', sinkTemperature=AMBIENT_TEMP, sinkAmplitude='', sinkDistributionType=UNIFORM, sinkFieldName='')
        if layer != NB_LAYERS-1:
            mdb.models['Model-1'].interactions[f'Radiation-{layer}'].deactivate(f'Layer-{layer+1}Step-0')
            mdb.models['Model-1'].interactions[f'Convection-{layer}'].deactivate(f'Layer-{layer+1}Step-0')

def AmbientTemp7():
# Initialiser la temperature ambiente dans la phase initiale
    a = mdb.models['Model-1'].rootAssembly

    # Recuperer toutes les regions
    cells_all = []
    faces_all = []
    edges_all = []
    vertices_all = []
    for inst in a.instances.values():
        cells_all += inst.cells[:]
        faces_all += inst.faces[:]
        edges_all += inst.edges[:]
        vertices_all += inst.vertices[:]
    region = regionToolset.Region(vertices=part.VertexArray(vertices_all), edges=part.EdgeArray(edges_all), faces=part.FaceArray(faces_all), cells=part.CellArray(cells_all))
    
    # Appliquer la temperature ambiente
    mdb.models['Model-1'].Temperature(name='InitialTemp', createStepName='Initial', region=region, distributionType=UNIFORM, crossSectionDistribution=CONSTANT_THROUGH_THICKNESS, magnitudes=(800,))


def LayerElementType8():
# Definir le format de notre simulation pour une couche (ici heat transfer)
    elemType1 = mesh.ElemType(elemCode=DC3D8, elemLibrary=STANDARD)
    elemType2 = mesh.ElemType(elemCode=DC3D6, elemLibrary=STANDARD)
    elemType3 = mesh.ElemType(elemCode=DC3D4, elemLibrary=STANDARD)
    for layer in range(NB_LAYERS):
        p = mdb.models['Model-1'].parts[f'Layer-{layer}']
        cells = p.cells[:]
        p.setElementType(regions=(cells,), elemTypes=(elemType1, elemType2, elemType3))


def SupportElementType9():
# Definir le format de notre simulation (ici heat transfer)
    elemType1 = mesh.ElemType(elemCode=DC3D8, elemLibrary=STANDARD)
    elemType2 = mesh.ElemType(elemCode=DC3D6, elemLibrary=STANDARD)
    elemType3 = mesh.ElemType(elemCode=DC3D4, elemLibrary=STANDARD)
    p = mdb.models['Model-1'].parts['Support']
    cells = p.cells[:]
    p.setElementType(regions=(cells,), elemTypes=(elemType1, elemType2, elemType3))

def Job10():
# Define Field Outputs & Create Job
    a = mdb.models['Model-1']

    # Field Outputs #
    for req in list(a.fieldOutputRequests.keys()): # delete all previous outputs
        del a.fieldOutputRequests[req]
    a.FieldOutputRequest(name='Thermal', createStepName='Death', variables=('NT', 'SDV'))
    
    # Job #
    a.rootAssembly.regenerate()
    mdb.Job(name=JOB, model='Model-1', description='', type=ANALYSIS, 
        atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
        memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
        explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
        modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
        scratch='', resultsFormat=ODB, numThreadsPerMpiProcess=1, 
        multiprocessingMode=DEFAULT, numCpus=NB_CPU, numDomains=NB_CPU, numGPUs=0)
    mdb.jobs[JOB].writeInput(consistencyChecking=OFF)


Constant()
Section()
Material()
SupportCreate()
LayerCreate1()
Assembly2()
Steps3()
Surface4()
HeatTransfer5()
AmbientTemp7()
LayerElementType8()
SupportElementType9()
Job10()