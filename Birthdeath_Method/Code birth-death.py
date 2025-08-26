# -*- coding: mbcs -*-
# Do not delete the following import lines
from abaqus import *

from abaqusConstants import *
from abaqus import backwardCompatibility
backwardCompatibility.setValues(reportDeprecated=False)
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

# CONSTANTES #
ABSOLUT_ZERO = -273.15                         # Absolute zero
STEFANBOLTZMANN = 5.67e-8                      # Stefan–Boltzmann constant
R = 8.314                                      # Ideal gas constant

# SIMULATION #
NB_INCREMENT = 1000                            # Max number of increments per step
MIN_INCR = 1e-10                               # Incrément minimum
INITIAL_INCR = 1e-6                            # Increment de temperature initial
MAXTEMP_INCR = 1000                            # Temperature maximale que peut atteindre un increment

# LASER #
NB_LASER = 3                                   # Number of laser passes per side
POWDER_SIZE = 30e-06                           # Size of a powder grain
H_SPACING = 100e-06                            # Spacing between each laser pass

LASER_POWER = 20                              # Laser power
LASER_RADIUS = 40e-06                          # Laser radius
LASER_SPEED = 600e-3                           # Laser speed
LASER_TIME = H_SPACING/LASER_SPEED             # Time period on a cell (size h_spacing*h_spacing)
COOLING_TIME = 8                               # Cooling time between each layer
#laser_formula = (2*LASER_POWER/(math.pi*(LASER_RADIUS**2)))*math.exp(-2*(X**2+Z**2)/(LASER_RADIUS**2))   # Gaussian distribution

# LAYER #
NB_LAYERS = 4 
LAYER_THICKNESS = 30e-06                       # Layer thickness (multiple of powder_size)
LAYER_SIZE = NB_LASER*H_SPACING                # Layer width, can be a direct size but must be multiple of POWDER_SIZE
NB_CELLS = NB_LASER * NB_LASER                 # Total number of cells per layer
LAYER_FILMCOEFF = 242.0                        # Layer surface coefficient
LAYER_EMISSIVITY = 0.75                        # Layer emissivity coefficient

# SUPPORT #
SUPPORT_SIZE = 1e-03                           # Support width
SUPPORT_THICKNESS = 5e-04                      # Support thickness
SUPPORTMESH_SIZE = 5e-04                       # Support mesh size
SUPPORT_FILMCOEFF = 10.0                       # Support surface coefficient
SUPPORT_EMISSIVITY = 0.7                       # Support emissivity coefficient

# EXTERIOR #
AMBIENT_TEMP = 25                              # Ambient temperature in °C


############## ALGORITHMES ###############

def find_top_face(s1):
# Trouver la face superieure de chaque cellule
    for face in s1.faces:
        normal = face.getNormal()
        if normal == (0,1,0) :
            return face

def find_lower_face(s1):
# Trouver la face inférieure de chaque cellule
    for face in s1.faces:
        normal = face.getNormal()
        if normal == (0,-1,0) :
            return face

def find_contour_face(s1):
# Trouver les faces extérieures de chaque cellule
    tab = []
    for face in s1.faces:
        normal = face.getNormal()
        if normal == (1,0,0) or normal == (-1,0,0)  or normal == (0,0,1) or normal == (0,0,-1) :
            tab.append(face)
    return tab


def Cell(layer):
# Permet de creer un tableau avec les cellules ordonnees

    p = mdb.models['Model-1'].parts[f'Layer-{layer}']

    x0 = H_SPACING/2 
    y0 = LAYER_THICKNESS*(layer+1)
    z0 = H_SPACING/2
    tab_cells = [None] * NB_CELLS
    index = 0

    for i in range(NB_LASER//2):
        # Prendre les cellules dans l'ordre pour ensuite pouvoir activer les sets dans le bon ordre
        for j in range(NB_LASER):
                # Gauche a droite
            cell = p.cells.findAt((x0+2*i*H_SPACING,y0,z0+j*H_SPACING))
            tab_cells[index] = cell # Permet d'ordonner les cellules
            index += 1
        for k in range(NB_LASER - 1, -1, -1):
                # Droite a gauche
            cell = p.cells.findAt((x0+(2*i+1)*H_SPACING,y0,z0+k*H_SPACING))
            tab_cells[index] = cell # Permet d'ordonner les cellules
            index += 1
    if(NB_LASER % 2 == 1): #Si le nb de cases par cote n'est pas pair, il faut refaire un passage
        for j in range(NB_LASER):
                # Gauche a droite
            cell = p.cells.findAt((x0+(NB_LASER-1)*H_SPACING,y0,z0+j*H_SPACING))
            tab_cells[index] = cell # Permet d'ordonner les cellules
            index += 1
    
    return tab_cells

############################################

def Constant():
# Definir les constantes universelles du modele
    mdb.models['Model-1'].setValues(absoluteZero=ABSOLUT_ZERO, stefanBoltzmann=STEFANBOLTZMANN, universalGas=R)


def Material():
# Creer les materiaux du support et layer
    ## Stainless Steel 316L ##
    mdb.models['Model-1'].Material(name='SS-316L') 
    mdb.models['Model-1'].materials['SS-316L'].Density(table=((7966.0, ), ))
    mdb.models['Model-1'].materials['SS-316L'].Conductivity(temperatureDependency=ON, table=((14.12, 20), (15.26, 100), (16.69, 200), (18.11, 300), (19.54, 400), (20.96, 500), (22.38, 600), (23.81, 700), (25.23, 800), (26.66, 900), (28.08, 1000), (29.50, 1100), (30.93, 1200), (32.35, 1300), (33.78, 1400)))
    mdb.models['Model-1'].materials['SS-316L'].Expansion(temperatureDependency=ON, table=((14.56, 20), (15.39, 100), (16.21, 200), (16.86, 300), (17.37, 400), (17.78, 500), (18.12, 600), (18.43, 700), (18.72, 800), (18.99, 900), (19.27, 1000), (19.53, 1100), (19.79, 1200), (20.02, 1300), (20.21, 1400)))
    mdb.models['Model-1'].materials['SS-316L'].SpecificHeat(temperatureDependency=ON, table=((0.492, 20), (0.502, 100), (0.514, 200), (0.526, 300), (0.538, 400), (0.550, 500), (0.562, 600), (0.575, 700), (0.587, 800), (0.599, 900), (0.611, 1000), (0.623, 1100), (0.635, 1200), (0.647, 1300), (0.659, 1400)))
    
    ## Layer (Inconel 718) ##
    mdb.models['Model-1'].Material(name='Inconel 718')
    mdb.models['Model-1'].materials['Inconel 718'].Density(table=((8190.0, ), ))
    mdb.models['Model-1'].materials['Inconel 718'].Conductivity(
        temperatureDependency=ON, table=((11.4, 21.0), (12.5, 100.0), (14.0, 
        300.0), (15.5, 500.0), (21.5, 700.0), (31.3, 1350.0)))
    mdb.models['Model-1'].materials['Inconel 718'].SpecificHeat(
        temperatureDependency=ON, table=((427.0, 21.0), (441.0, 100.0), (481.0, 
        300.0), (521.0, 500.0), (601.0, 700.0), (691.0, 1350.0)))


def Section():
# Creer les sections 'Support' et 'Layer' pour ensuite assigner les materiaux
    mdb.models['Model-1'].HomogeneousSolidSection(name='Support', material='SS-316L', thickness=None)
    mdb.models['Model-1'].HomogeneousSolidSection(name='Layer', material='SS-316L', thickness=None)

    
def SupportCreate():
# Creation du support & mesh

    ## Creer le support ##
    indent = (LAYER_SIZE-SUPPORT_SIZE)/2 #Pour bien placer le support au milieu des layers des le debut
    s1 = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=0.01)
    s1.sketchOptions.setValues(decimalPlaces=4)
    s1.rectangle(point1=(indent, -SUPPORT_THICKNESS), point2=(SUPPORT_SIZE+indent, 0))
    p = mdb.models['Model-1'].Part(name='Support', dimensionality=THREE_D, type=DEFORMABLE_BODY)
    p.BaseSolidExtrude(sketch=s1, depth=SUPPORT_SIZE)
    s1.unsetPrimaryObject()
    del mdb.models['Model-1'].sketches['__profile__']

    ## Ajouter le materiau au support ##
    # ! Avoir creer la section 'Support' avant
    region = regionToolset.Region(cells=p.cells[:])
    p.SectionAssignment(region=region, sectionName='Support', offset=0.0, offsetType=MIDDLE_SURFACE, offsetField='', thicknessAssignment=FROM_SECTION)

    ## Creer le mesh du support ##
    p = mdb.models['Model-1'].parts['Support']
    p.seedPart(size=SUPPORTMESH_SIZE, deviationFactor=0.1, minSizeFactor=0.1)
    p.generateMesh()


def LayerCreate1():
# Creation d'une couche et des cells
    for layer in range(NB_LAYERS):
        ## Creation de la couche ##
        s1 = mdb.models['Model-1'].ConstrainedSketch(name='__profile__', sheetSize=0.01)
        s1.sketchOptions.setValues(decimalPlaces=4)
        s1.rectangle(point1=(0, layer*LAYER_THICKNESS), point2=(LAYER_SIZE, LAYER_THICKNESS+layer*LAYER_THICKNESS))
        p = mdb.models['Model-1'].Part(name=f'Layer-{layer}', dimensionality=THREE_D, type=DEFORMABLE_BODY)
        p.BaseSolidExtrude(sketch=s1, depth=LAYER_SIZE)

        s1.unsetPrimaryObject()
        del mdb.models['Model-1'].sketches['__profile__']

        ## Ajouter le materiau a la couche ##
        # ! Avoir creer la section 'Layer' avant
        region = regionToolset.Region(cells=p.cells[:])
        p.SectionAssignment(region=region, sectionName='Layer', offset=0.0, offsetType=MIDDLE_SURFACE, offsetField='', thicknessAssignment=FROM_SECTION)

        ## Decoupage en plans ##    
        ref_face = p.faces

        for i in range(1, NB_LASER):
            offset_val = i * H_SPACING
            # Decoupage x 
            p.DatumPlaneByOffset(plane=ref_face[4], flip=SIDE2, offset=offset_val)
            # Decoupage z
            p.DatumPlaneByOffset(plane=ref_face[0], flip=SIDE2, offset=offset_val)
        
        ## Decoupage des plans en cells ##
        datum_planes = [d for d in p.datums.values()]
        for plane in datum_planes:
            pickedCells = p.cells[:]
            p.PartitionCellByDatumPlane(datumPlane=plane, cells=pickedCells)

        ## Mesh ##
        p.seedPart(size=POWDER_SIZE, deviationFactor=0.1, minSizeFactor=0.1)
        p.generateMesh()

        ## Creation des sets ##
        tab_cells = Cell(layer)

        for i in range(NB_CELLS):
            cell = tab_cells[i]

            # Creation des sets par cell
            set_name = f'Layer-{layer}Set-{i}'
            p.Set(elements=cell.getElements(), name=set_name)

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
    a.translate(instanceList=('Support-1', ), vector=(0.0, 0.0, (LAYER_SIZE-SUPPORT_SIZE)/2)) # Centrer la layer
    

def Steps3():
# Creer les activations de cells a chaque etape
    a = mdb.models['Model-1'].rootAssembly
    a.regenerate()
    
    ## Initial : Toutes cells eteintes ##
    mdb.models['Model-1'].HeatTransferStep(name='Death', previous='Initial', timePeriod=0.1, initialInc=0.01, minInc=1e-8, maxInc=0.1, deltmx=350.0) # Creer l'etape initiale
    region =a.sets['AllLayers'] 
    mdb.models['Model-1'].ModelChange(name='Desactive', createStepName='Death', region=region, regionType=ELEMENTS, activeInStep=False, includeStrain=False) # Desactiver toutes les cells au debut

    for layer in range(NB_LAYERS):
        ## Etape 1 : Active la premiere cell ##
        if layer == 0 :
            mdb.models['Model-1'].HeatTransferStep(name='Layer-0Step-0', previous='Death', timePeriod=LASER_TIME, maxNumInc=NB_INCREMENT, initialInc=INITIAL_INCR, minInc=MIN_INCR, maxInc=LASER_TIME/2, deltmx=MAXTEMP_INCR)
            region =a.instances['Layer-0'].sets['Layer-0Set-0']
            mdb.models['Model-1'].ModelChange(name='Activate-00', createStepName='Layer-0Step-0', region=region, regionType=ELEMENTS, activeInStep=True, includeStrain=True)
        else :
            mdb.models['Model-1'].HeatTransferStep(name=f'Layer-{layer}Step-0', previous=f'Cooling-{layer-1}', timePeriod=LASER_TIME, maxNumInc=NB_INCREMENT, initialInc=INITIAL_INCR, minInc=MIN_INCR, maxInc=LASER_TIME/2, deltmx=MAXTEMP_INCR)
            region =a.instances[f'Layer-{layer}'].sets[f'Layer-{layer}Set-0']
            mdb.models['Model-1'].ModelChange(name=f'Activate-{layer}0', createStepName=f'Layer-{layer}Step-0', region=region, regionType=ELEMENTS, activeInStep=True, includeStrain=True)
        
        ## Etapes n : Activer les cells une a une ##
        for i in range(1, NB_CELLS):
            mdb.models['Model-1'].HeatTransferStep(name=f'Layer-{layer}Step-{i}', previous=f'Layer-{layer}Step-{i-1}', timePeriod=LASER_TIME, maxNumInc=NB_INCREMENT, initialInc=INITIAL_INCR, minInc=MIN_INCR, maxInc=LASER_TIME/2, deltmx=MAXTEMP_INCR)
            region =a.instances[f'Layer-{layer}'].sets[f'Layer-{layer}Set-{i}']
            mdb.models['Model-1'].ModelChange(name=f'Activate-{layer}{i}', createStepName=f'Layer-{layer}Step-{i}', region=region, regionType=ELEMENTS, activeInStep=True, includeStrain=True)

        ## Cooling entre chaque couche
        mdb.models['Model-1'].HeatTransferStep(name=f'Cooling-{layer}', previous=f'Layer-{layer}Step-{NB_CELLS-1}', timePeriod=COOLING_TIME, maxNumInc=NB_INCREMENT, initialInc=INITIAL_INCR, minInc=MIN_INCR, maxInc=COOLING_TIME/2, deltmx=MAXTEMP_INCR)


def Surface4():
# Relier les couches entre elles
    a = mdb.models['Model-1'].rootAssembly
    
    ## Lier face support et face layer-0 ##
    s1 = a.instances['Support-1']
    face1 = find_top_face(s1)
    region1=regionToolset.Region(side1Faces=part.FaceArray((face1,)))
    s2 = a.instances['Layer-0']
    face2 = find_lower_face(s2)
    region2=regionToolset.Region(side1Faces=part.FaceArray((face2,)))
    mdb.models['Model-1'].Tie(name='Surface-Support', main=region1, secondary=region2, positionToleranceMethod=COMPUTED, adjust=ON, tieRotations=ON, thickness=ON)

    ## Lier les autres layers ##
    for layer in range(NB_LAYERS-1):
        s1 = a.instances[f'Layer-{layer}']
        face1 = find_top_face(s1)
        region1=regionToolset.Region(side1Faces=part.FaceArray((face1,)))
        s2 = a.instances[f'Layer-{layer+1}']
        face2 = find_lower_face(s2)
        region2=regionToolset.Region(side1Faces=part.FaceArray((face2,)))
        mdb.models['Model-1'].Tie(name=f'Layer{layer}-Layer{layer+1}', main=region1, secondary=region2, positionToleranceMethod=COMPUTED, adjust=ON, tieRotations=ON, thickness=ON)


def HeatTransfer5():
    # Definir les differentes rayonnements du support & couches 
    a = mdb.models['Model-1'].rootAssembly

    ## SUPPORT ##
    support = a.instances['Support-1'].faces
    side1Faces1 = support.getSequenceFromMask(mask=('[#37 ]', ), )
    region=regionToolset.Region(side1Faces=side1Faces1)
    
    ## Film Condition ##
    mdb.models['Model-1'].FilmCondition(name='SupportSurfaceFilm', createStepName='Death', surface=region, definition=EMBEDDED_COEFF, filmCoeff=SUPPORT_FILMCOEFF, filmCoeffAmplitude='', sinkTemperature=AMBIENT_TEMP, sinkAmplitude='', sinkDistributionType=UNIFORM, sinkFieldName='')
    ## Radiation ##
    mdb.models['Model-1'].RadiationToAmbient(name='SupportSurfaceRadiation', createStepName='Death', surface=region, radiationType=AMBIENT, distributionType=UNIFORM, field='', emissivity=SUPPORT_EMISSIVITY, ambientTemperature=AMBIENT_TEMP, ambientTemperatureAmp='')

    ## LAYERS ##

    for layer in range(NB_LAYERS):
        all_faces = a.instances[f'Layer-{layer}']
        contour_faces = find_contour_face(all_faces)
        regioncontour = regionToolset.Region(side1Faces=part.FaceArray(contour_faces))
        top_face = find_top_face(all_faces)
        regiontop = regionToolset.Region(side1Faces=part.FaceArray((top_face,)))

        ## Film Condition ##
        mdb.models['Model-1'].FilmCondition(name=f'Layer-{layer}ContourFilm', createStepName=f'Layer-{layer}Step-0', surface=regioncontour, definition=EMBEDDED_COEFF, filmCoeff=LAYER_FILMCOEFF, filmCoeffAmplitude='', sinkTemperature=AMBIENT_TEMP, sinkAmplitude='', sinkDistributionType=UNIFORM, sinkFieldName='')
        mdb.models['Model-1'].FilmCondition(name=f'Layer-{layer}TopFilm', createStepName=f'Layer-{layer}Step-0', surface=regiontop, definition=EMBEDDED_COEFF, filmCoeff=LAYER_FILMCOEFF, filmCoeffAmplitude='', sinkTemperature=AMBIENT_TEMP, sinkAmplitude='', sinkDistributionType=UNIFORM, sinkFieldName='')

        ## Radiation ##
        mdb.models['Model-1'].RadiationToAmbient(name=f'Layer-{layer}ContourRadiation', createStepName=f'Layer-{layer}Step-0', surface=regioncontour, radiationType=AMBIENT, distributionType=UNIFORM, field='', emissivity=LAYER_EMISSIVITY, ambientTemperature=AMBIENT_TEMP, ambientTemperatureAmp='')
        mdb.models['Model-1'].RadiationToAmbient(name=f'Layer-{layer}TopRadiation', createStepName=f'Layer-{layer}Step-0', surface=regiontop, radiationType=AMBIENT, distributionType=UNIFORM, field='', emissivity=LAYER_EMISSIVITY, ambientTemperature=AMBIENT_TEMP, ambientTemperatureAmp='')


def Load6():
# Creer un repere au centre sur la face du haut de chaque cell
    def find_top_face2(cell_faces,layer):
    # Trouver la face superieure de chaque cellule
        s1 = a.instances[f'Layer-{layer}']
        for face_id in cell_faces:
            normal = s1.faces[face_id].getNormal()
            if normal == (0,1,0) :
                return s1.faces[face_id]

    for layer in range(NB_LAYERS):
        tab_cells = Cell(layer) 
        a = mdb.models['Model-1'].rootAssembly
        a.regenerate()

        for k in range(NB_CELLS):
            face_cell = tab_cells[k].getFaces() # Recupere toutes les faces d'une cellule
            face = find_top_face2(face_cell, layer) # Selectionne face du dessus
            ## Trouver le centre de la face ##
            center = face.getCentroid()
            ## Creer le repere ##
            cys = a.DatumCsysByThreePoints(name=f'Coord-{layer}-{k}', coordSysType=CARTESIAN, origin=center[0], line1=(1.0, 0.0, 0.0), line2=(0.0, 1.0, 0.0))
            ## Creer la repartition de chaleur du laser en ce repere ##
            mdb.models['Model-1'].ExpressionField(name=f'Layer-{layer}Laser-{k}', localCsys=a.datums[cys.id], description='Laser Gaussien', expression=f'(2*{LASER_POWER}/(pi*({LASER_RADIUS}**2)))*exp((-2)*(X**2+Z**2)/({LASER_RADIUS}**2))')
            
            ## Appliquer le Load Heat Flux ##
            region = regionToolset.Region(side1Faces=part.FaceArray((face,)))
            mdb.models['Model-1'].SurfaceHeatFlux(name=f'Layer-{layer}Heat-{k}', createStepName=f'Layer-{layer}Step-{k}', region=region, magnitude=1.0, distributionType=FIELD, field=f'Layer-{layer}Laser-{k}')

            ## Desactiver la charge dans la step suivante sauf pour la derniere cellule ##
            if k == NB_CELLS-1:
                mdb.models['Model-1'].loads[f'Layer-{layer}Heat-{k}'].deactivate(f'Cooling-{layer}')
            else : 
                mdb.models['Model-1'].loads[f'Layer-{layer}Heat-{k}'].deactivate(f'Layer-{layer}Step-{k+1}')


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
    mdb.models['Model-1'].Temperature(name='AmbientTemp', createStepName='Initial', region=region, distributionType=UNIFORM, crossSectionDistribution=CONSTANT_THROUGH_THICKNESS, magnitudes=(AMBIENT_TEMP,))


def LayerElementType8():
    # Definir le format de notre simulation pour une couche (ici heat transfer)
    elemType1 = mesh.ElemType(elemCode=DCC3D8, elemLibrary=STANDARD)
    elemType2 = mesh.ElemType(elemCode=DC3D6, elemLibrary=STANDARD)
    elemType3 = mesh.ElemType(elemCode=DC3D4, elemLibrary=STANDARD)
    for layer in range(NB_LAYERS):
        p = mdb.models['Model-1'].parts[f'Layer-{layer}']
        cells = p.cells[:]
        p.setElementType(regions=(cells,), elemTypes=(elemType1, elemType2, elemType3))


def SupportElementType9():
    # Definir le format de notre simulation (ici heat transfer)
    elemType1 = mesh.ElemType(elemCode=DCC3D8, elemLibrary=STANDARD)
    elemType2 = mesh.ElemType(elemCode=DC3D6, elemLibrary=STANDARD)
    elemType3 = mesh.ElemType(elemCode=DC3D4, elemLibrary=STANDARD)
    p = mdb.models['Model-1'].parts['Support']
    cells = p.cells[:]
    p.setElementType(regions=(cells,), elemTypes=(elemType1, elemType2, elemType3))



Constant()
Section()
Material()
SupportCreate()
LayerCreate1()
Assembly2()
Steps3()
Surface4()
HeatTransfer5()
Load6()
AmbientTemp7()
LayerElementType8()
SupportElementType9()