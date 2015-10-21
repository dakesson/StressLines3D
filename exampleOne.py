from femLibrary.pycalfem_GeoData import *
from femLibrary.pycalfem_mesh import *
from femLibrary.pycalfem import *
from femLibrary.pycalfem_utils import *

from vtk import *
import numpy as np

np.set_printoptions(precision=2)

class VTKWriter:

    def __init__(self):
        self.renWin = vtk.vtkRenderWindow()
        self.renWin.SetSize(600, 600)
     
        self.iRen = vtk.vtkRenderWindowInteractor()
        self.iRen.SetRenderWindow(self.renWin)
    
    def createInput(self, res=6):
        self.res = res
        g = GeoData()

        #Add Points:
        #warning no point nr 1
        g.addPoint([0, 0, 0], 0, marker = 1)
        g.addPoint([10, 0, 0], 2, marker = 1)
        g.addPoint([10, 10, 0], 3, marker = 2)
        g.addPoint([0, 10, 0], 4, marker = 2)
        g.addPoint([0, 0, 10], 5, marker = 1)
        g.addPoint([10, 0, 10], 6, marker = 1)
        g.addPoint([10, 10, 10], 7, marker = 2)
        g.addPoint([0, 10, 10], 8, marker = 2)
        
        #Add splines:
        elOnCurveVar = self.res
        g.addSpline([0, 2], 0, elOnCurve = elOnCurveVar)
        g.addSpline([2, 3], 1, elOnCurve = elOnCurveVar)
        g.addSpline([3, 4], 2, elOnCurve = elOnCurveVar)
        g.addSpline([4, 0], 3, elOnCurve = elOnCurveVar)
        g.addSpline([0, 5], 4, elOnCurve = elOnCurveVar)
        g.addSpline([2, 6], 5, elOnCurve = elOnCurveVar)
        g.addSpline([3, 7], 6, elOnCurve = elOnCurveVar)
        g.addSpline([4, 8], 7, elOnCurve = elOnCurveVar)
        g.addSpline([5, 6], 8, elOnCurve = elOnCurveVar)
        g.addSpline([6, 7], 9, elOnCurve = elOnCurveVar)
        g.addSpline([7, 8], 10, elOnCurve = elOnCurveVar)
        g.addSpline([8, 5], 11, elOnCurve = elOnCurveVar)
        
        #Add surfaces:
        g.addStructuredSurface([0, 1, 2, 3], 0)
        g.addStructuredSurface([8, 9, 10, 11], 1)
        g.addStructuredSurface([0, 4, 8, 5], 2) #Floow
        g.addStructuredSurface([1, 5, 9, 6], 3)
        g.addStructuredSurface([2, 6, 10, 7], 4, marker=2)
        g.addStructuredSurface([3, 4, 11, 7], 5)
        
        #Add Volume:
        #addStructuredVolume() takes three args. The first is a list of surface IDs (structured surfaces).
        # The surfaces should make a hexahedron (i.e. 6 surfaces). Other kinds of structured volumes than hexahedra will
        # not work for hexahedral elements, which is the only type of 3D element that CALFEM handles.
        #The two optional parameters are the volume ID and volume marker.
        g.addStructuredVolume([0,1,2,3,4,5], 0, marker=2)
        
        elType = 5 #Element type 5 is hexahedron. (See user manual for more element types)
        self.dofsPerNode= 3 #Degrees of freedom per node.
        
        mesher = GmshMesher(geoData = g,
                            gmshExecPath = "Gmsh.app/Contents/MacOS/gmsh", #Path to gmsh.exe. If None then the system PATH variable is queried. Both relative and absolute paths work, like "gmsh\gmsh.exe".
                            elType = elType, 
                            dofsPerNode= self.dofsPerNode)
        
        #Mesh the geometry:
        E=200e9
        v=0.3
        ptype = 4 #for 3D
        
        self.ep=ptype
        self.D=hooke(ptype,E,v)        
        
        self.coords, self.edof, self.dofs, self.bdofs, self.elementmarkers = mesher.create()
        self.ex, self.ey, self.ez = coordxtr(self.edof, self.coords, self.dofs)
        
        self.numEl = size(self.ex, 0)
        self.nodesPerEl = size(self.ex,1)
        self.numNodes = size(self.coords,0)
        self.dofPerNode = 3
        
    def solveProblem(self):        
        self.nDofs = self.edof.max() 
        self.K = zeros([self.nDofs,self.nDofs])
        self.f = zeros([self.nDofs,1])

        #Add gravity
        eq = array([0,0,0])
       
        for eltopo, elx, ely, elz, elMarker in zip(self.edof, self.ex, self.ey,self.ez, self.elementmarkers):
            #Calc element stiffness matrix: Conductivity matrix D is taken 
            # from Ddict and depends on which region (which marker) the element is in.
            ke, fe = soli8e(elx,ely,elz,self.ep,self.D,eq)
            assem(eltopo, self.K, ke, self.f, fe)
        
        print "Solving equation system..."

        self.bc = array([],'i')
        bcVal = array([],'i')


        self.bc, self.bcVal = applybc(self.bdofs, self.bc, bcVal, 1, 0.0, 0)
#        applyforce(self.bdofs, self.f, 2, 10e5, 2)

#        self.addAdjecentBC()
#        print self.bc

        #Add force
        self.f[subtract(self.bdofs[2][1::3],1)] = -1000
#        self.f[1::3] = 1000


        self.a,self.r = solveq(self.K,self.f,self.bc,self.bcVal)
        self.ed = extractEldisp(self.edof, self.a)
        
    def addAdjecentBC(self):
        for bc in self.bc:
            nodeNr, direction = dofToNode(bc, self.dofsPerNode)
            nodeCoord = self.coords[nodeNr,:]
            print nodeCoord            

            for i in range(self.numNodes):
                coord = self.coords[i,:]
                distance = linalg.norm(coord-nodeCoord)


                if abs(distance) < 1.2 and distance != 0:
                    
                    self.bc = np.append(self.bc,[(i)*self.dofsPerNode+1])
                    self.bc = np.append(self.bc,[(i)*self.dofsPerNode+2])
                    self.bc = np.append(self.bc,[(i)*self.dofsPerNode+3])
                    
                    self.bcVal = np.append(self.bcVal,[0])        
                    self.bcVal = np.append(self.bcVal,[0])        
                    self.bcVal = np.append(self.bcVal,[0])                            

    def makeHexahedron(self):
    
        uGrid = vtk.vtkUnstructuredGrid()

    
    
        #Setup points
        points = vtk.vtkPoints()    
        scalars = vtk.vtkFloatArray()
                        
        for i in range(self.numNodes):
            points.InsertNextPoint(self.coords[i,0],self.coords[i,1],self.coords[i,2])
            
        uGrid.SetPoints(points)
        
        self.gaussPoints = vtk.vtkPoints()

    
        #Setup topology
        self.es = zeros([self.numEl,6])        
        self.globalEsSmooth = zeros([self.numEl,self.nodesPerEl,6])
        
        for i in range(self.numEl):
            hex_ = vtk.vtkHexahedron()

            #Compute stress for every element
            elEs, eci = soli8s(self.ex[i,:],self.ey[i,:],self.ez[i,:],self.ep,self.D,asmatrix(self.ed[i,:]))
            self.globalEsSmooth[i,:,:] = elEs
            self.es[i,:] = [average(elEs[:,0]), average(elEs[:,1]) ,average(elEs[:,2]), average(elEs[:,3]), average(elEs[:,4]), average(elEs[:,5])]           
           
            mises = (0.5*((average(elEs[:,0])- average(elEs[:,1]))**2 + (average(elEs[:,1]) - average(elEs[:,2]))**2 + (average(elEs[:,0]) - average(elEs[:,2]))**2 + 6*(average(elEs[:,4])**2 + average(elEs[:,5])**2 + average(elEs[:,3])**2)))
            scalars.InsertNextValue(mises)       
          
            #Getting the node number from edof
            for j in range(0, self.nodesPerEl):    
                hex_.GetPointIds().SetId(j, self.edof[i,(self.dofPerNode-1+j*self.dofPerNode)]/self.dofPerNode-1)
                #self.gaussPoints.InsertNextPoint(self.gpCoord[:,j])                
            
            uGrid.InsertNextCell(hex_.GetCellType(), hex_.GetPointIds())
            uGrid.GetCellData().SetScalars(scalars) 
            
        return uGrid

    def displayBodies(self):

        self.grid = self.makeHexahedron()

        a,b = self.grid.GetScalarRange()
        xmi, xma, ymi, yma, zmi, zma = self.grid.GetBounds()
        print("Scalar range:" + str(a) + "   " + str(b))

        lut=vtkColorTransferFunction()
        lut.AddRGBPoint(a,0,0,1)
        lut.AddRGBPoint(b,1,0,0)        
        
        ctf = vtkLookupTable()
        ctf.SetHueRange(0.667,1.0 )
        ctf.SetValueRange(1.0, 1.0)
        ctf.SetSaturationRange(1.0, 1.0)
        ctf.SetTableRange(a,b)

     
        # Create and link the mappers actors and renderers together.
        volumeMapper = vtk.vtkDataSetMapper()
        volumeMapper.SetInput(self.grid)
        volumeMapper.SetLookupTable(lut)
  
        volumeActor = vtk.vtkActor()
        volumeActor.SetMapper(volumeMapper)
        volumeActor.GetProperty().SetOpacity(0.3)  
        
        self.ren = vtkRenderer()
        self.ren.SetBackground(1.0, 1.0, 1.0)        
        self.ren.AddViewProp(volumeActor)
        
    def drawForceArrows(self):
        
        #Force arrow
        self.forcePoints = vtk.vtkPoints()
        self.forceVectors = vtk.vtkDoubleArray()
        self.forceVectors.SetNumberOfComponents(3)
        self.forceScalars = vtk.vtkFloatArray()
        
                #Set up force arrows
        for forceScalar, dofNr in zip(self.f, range(size(self.f,0))):
            if forceScalar != 0:

                nodeNr, direction = dofToNode(dofNr, self.dofsPerNode)
                self.forcePoints.InsertNextPoint(add(self.coords[nodeNr-1,:],[0,1.0,0]))
                
                #print "Dof nr: " + str(dofNr) + "  Node nr: " + str(nodeNr) + "Direction: " + str(direction)
                
                if direction == 0:
                    self.forceVectors.InsertNextTuple3(1,0,0)
                elif direction == 1:
                    self.forceVectors.InsertNextTuple3(0,1,0)
                elif direction == 2:
                    self.forceVectors.InsertNextTuple3(0,0,1)                
            
                self.forceScalars.InsertNextValue(forceScalar)
        
        
        dataForce = vtk.vtkPolyData()
        dataForce.SetPoints(self.forcePoints)
        dataForce.GetPointData().SetScalars(self.forceScalars) 
        dataForce.GetPointData().SetVectors(self.forceVectors)
        
        c,d = dataForce.GetScalarRange()

        forceArrow = vtkArrowSource()
        forceArrow.SetTipRadius(0.25)
        forceArrow.SetShaftRadius(0.15)
        
        forceArrowGlyph = vtkGlyph3D()
        forceArrowGlyph.SetInput(dataForce)
        forceArrowGlyph.SetSource(forceArrow.GetOutput())
        forceArrowGlyph.SetScaleFactor(1.0/max([abs(c),abs(d)]))
        
        forceArrowMapper = vtkPolyDataMapper()
        forceArrowMapper.SetInputConnection(forceArrowGlyph.GetOutputPort())
        forceArrowActor = vtkActor()
        forceArrowActor.SetMapper(forceArrowMapper) 
        forceArrowActor.GetProperty().SetColor(0.0,1.0,0.0)
        
        self.ren.AddActor(forceArrowActor)
        
    def drawBoundary(self):
        #Bc points
        self.bcPoints = vtk.vtkPoints()
        
        #Calculate the node nr from dof number to get coordinates
        for bcPoint in self.bc[self.dofsPerNode-1::self.dofsPerNode]/self.dofsPerNode-1:
            self.bcPoints.InsertNextPoint(self.coords[bcPoint,:])
            
        dataBc = vtk.vtkPolyData()
        dataBc.SetPoints(self.bcPoints)
        
        #Draw BC balls
        bcball = vtk.vtkSphereSource()
        bcball.SetRadius(0.2)
        bcball.SetThetaResolution(8)
        bcball.SetPhiResolution(8)
        
        bcballGlyph = vtk.vtkGlyph3D()
        bcballGlyph.SetInput(dataBc)
        bcballGlyph.SetSourceConnection(bcball.GetOutputPort())
        bcballGlyph.SetScaleModeToScaleByScalar()
        bcballGlyph.SetScaleModeToDataScalingOff()

        bcballMapper = vtkPolyDataMapper()
        bcballMapper.SetInputConnection(bcballGlyph.GetOutputPort())
        
        bcballActor = vtkActor()
        bcballActor.SetMapper(bcballMapper)         
        bcballActor.GetProperty().SetColor(0.0,1.0,0.0)   

        self.ren.AddActor(bcballActor)     
        
    def drawGaussPoints(self):
        dataGp = vtk.vtkPolyData()
        dataGp.SetPoints(self.gaussPoints)
        
        #Draw BC balls
        gpball = vtk.vtkSphereSource()
        gpball.SetRadius(0.01)
        gpball.SetThetaResolution(8)
        gpball.SetPhiResolution(8)
        
        gpballGlyph = vtk.vtkGlyph3D()
        gpballGlyph.SetInput(dataGp)
        gpballGlyph.SetSourceConnection(gpball.GetOutputPort())
        gpballGlyph.SetScaleModeToScaleByScalar()
        gpballGlyph.SetScaleModeToDataScalingOff()

        gpballMapper = vtkPolyDataMapper()
        gpballMapper.SetInputConnection(gpballGlyph.GetOutputPort())
        
        gpballActor = vtkActor()
        gpballActor.SetMapper(gpballMapper)         
        gpballActor.GetProperty().SetColor(1.0,1.0,1.0)   

        self.ren.AddActor(gpballActor)    
        
    def cutPlane(self):
        #create a plane to cut,here it cuts in the XZ direction (xz normal=(1,0,0);XY =(0,0,1),YZ =(0,1,0)
        plane=vtk.vtkPlane()
        plane.SetOrigin(0.5,0,0)
        plane.SetNormal(1,0,0)
         
        #create cutter
        cutter=vtk.vtkCutter()
        cutter.SetCutFunction(plane)
        cutter.SetInput(self.grid)
        cutter.Update()
        cutterMapper=vtk.vtkPolyDataMapper()
        cutterMapper.SetInputConnection( cutter.GetOutputPort())
         
        #create plane actor
        planeActor=vtk.vtkActor()
        planeActor.GetProperty().SetColor(1.0,0.5,0)
        planeActor.GetProperty().SetLineWidth(2)
        planeActor.SetMapper(cutterMapper)
        
        self.ren.AddActor(planeActor)
        

    def drawFlowLines(self):
        
        print "Drawing flow lines"
        xRes = self.res*1
        yRes = self.res*1
        zRes = self.res*1

        self.nrFlowLines = 6        
        
        xmin,xmax, ymin,ymax, zmin,zmax = self.grid.GetBounds()
        
        self.xStepSize = (xmax - xmin) / (xRes)
        self.yStepSize = (ymax - ymin) / (yRes)
        self.zStepSize = (zmax - zmin) / (zRes)

        self.xRange = linspace(xmin + self.xStepSize/2, xmax - self.xStepSize/2, xRes)
        self.yRange = linspace(ymin + self.yStepSize/2, ymax - self.yStepSize/2, yRes)
        self.zRange = linspace(zmin + self.zStepSize/2, zmax - self.zStepSize/2, zRes)

        
#        maxForce = self.forceMatrix.max()
        
        allLines = []
        allForceInLines = []        
        
        for xlineNr in range(1,self.nrFlowLines+1):
            for zlineNr in range(1,self.nrFlowLines+1):
                
                forceLine = []
                forceInLine = 0
                
                #Make guess of x
                x = ((max(self.xRange) + self.xStepSize - min(self.xRange)) / (self.nrFlowLines+1)) * (xlineNr)
                #Loop over Y
                for y in self.yRange:                    
                    
                    #Find z first
                    z, zForce = self.findPositionIndex(2,x,y,0,zlineNr)
                    #Find x
                    x, xForce = self.findPositionIndex(0,0,y,z,xlineNr)  

                    #Find z again   
                    oldPos = [x,y,z]

                    i = 1
                    ok = True                    
                    while True:
                        z, zForce = self.findPositionIndex(2,x,y,0,zlineNr)
                        x, xForce = self.findPositionIndex(0,0,y,z,xlineNr)  
                        
                        dist = distance([x,y,z], oldPos)
                        if (dist < 0.1):
                            break
                        
                        if i > 10:
                            print "---DID NOT CONVERGE----"
                            print i
                            print dist
                            ok = False
                            break
                        
                        oldPos = [x,y,z]
                        i = i + 1
                        
                    if ok:
                        forceLine.append([x,y,z])
                        forceInLine = forceInLine + xForce + zForce
                        
                allLines.append(forceLine)
                allForceInLines.append(forceInLine)


        self.forceScalarRange = [min(allForceInLines), max(allForceInLines)]
        
        for forceLineLoop, forceInLineLoop in zip(allLines, allForceInLines):
            self.drawLine(forceLineLoop, forceInLineLoop)
                    
                
        
        
    def findPositionIndex(self,direction,x,y,z,lineNr):
        # Direction x = 0, y = 1, z = 2
        # Returns the index
    
        forceSum = 0.0
        numberPoints = 0
        forceVector = []
        stepSize = 0
        
        if direction == 0:
            for x in self.xRange:
                forceVector.append(self.zForceForCoord(x,y,z))
                coordRange = self.xRange
                
        elif direction == 1:
            for y in self.yRange:
                forceVector.append(self.zForceForCoord(x,y,z))
                coordRange = self.yRange
                
        elif direction == 2:
            for z in self.zRange:
                forceVector.append(self.zForceForCoord(x,y,z))
                coordRange = self.zRange

        forceSum  =sum(forceVector)
        numberPoints = len(forceVector)
                    
        #Find the coord for flowline
        accumulatedForce = 0.0
        forceStep = (forceSum / (self.nrFlowLines + 1))
        targetAccumulatedForce = forceStep * lineNr

        #Find z

        for index in range(numberPoints):
            accumulatedForce += forceVector[index]
            
#            print "accum: " + str(accumulatedForce) + " target: " + str(targetAccumulatedForce)
            
            #Its abs on all forces
            if accumulatedForce > targetAccumulatedForce:
                #print "y:" + str(y) + " accum: " + str(accumulatedForce) + " target: " + str(targetAccumulatedForce) + "line nr: " + str(lineNr) + "index: " + str(index)
                forceOver = accumulatedForce - targetAccumulatedForce
                relative = forceOver / forceVector[index]                    
                
                stepSize = coordRange[len(coordRange)-1] - coordRange[len(coordRange)-2]
                coord = coordRange[index] - stepSize * relative + stepSize/2

                return coord, forceStep

        print "----------------------"
        print "Failed to find position"
        print direction
        print accumulatedForce
        print targetAccumulatedForce
        print "----------------------"  
        
    def zForceForCoord(self,x,y,z):
        
            cellId = self.findCellFor([x,y,z])
            
            esSmooth = self.globalEsSmooth[cellId,:,:]
            zForce = abs(esSmooth[:,1])/self.yStepSize + abs(esSmooth[:,3])/self.xStepSize + 10*abs(esSmooth[:,5])/self.zStepSize                     

            zPointForce = soli8sPointInterpolate(zForce, self.ex[cellId,:], self.ey[cellId,:], self.ez[cellId,:], [x,y,z])                
#            print "----force"
#            print zPointForce
#            print "---slut force"
            if zPointForce < 0:
                print "--------------"
                print zPointForce
                print self.ex[cellId,:]
                print self.ey[cellId,:]
                print self.ez[cellId,:]
                print [x,y,z]
                print "--------------"
            
            return zPointForce
    
    
    
    def drawLine(self, lineIn, forceInLine):

        threshold = self.forceScalarRange[0] + (self.forceScalarRange[1] - self.forceScalarRange[0]) * 0.20
        
        if forceInLine < threshold:
            return
            
        # Create a vtkPoints object and store the points in it        
        points = vtk.vtkPoints()
        
        for point in lineIn:
            points.InsertNextPoint(point)
         
        # Create a cell array to store the lines in and add the lines to it
        lines = vtk.vtkCellArray()
        
        col = vtkDoubleArray()
        col.SetNumberOfComponents(1)
         
        for i in range(len(lineIn)-1):
          line = vtk.vtkLine()
          line.GetPointIds().SetId(0,i)
          line.GetPointIds().SetId(1,i+1)
          lines.InsertNextCell(line)

          
        for i in range(len(lineIn)):
          col.InsertNextTupleValue([forceInLine])

        lut=vtkColorTransferFunction()
        lut.AddRGBPoint(self.forceScalarRange[0],0,0,1)
        lut.AddRGBPoint(self.forceScalarRange[1],1,0,0)        
         
        # Create a polydata to store everything in
        linesPolyData = vtk.vtkPolyData()
         
        # Add the points to the dataset
        linesPolyData.SetPoints(points)
         
        # Add the lines to the dataset
        linesPolyData.SetLines(lines)
        linesPolyData.GetPointData().SetScalars(col)
        
        tubeFilter = vtkTubeFilter()
        tubeFilter.SetInput(linesPolyData)
        tubeFilter.SetRadius(0.15)
#        tubeFilter.SetVaryRadiusToVaryRadiusByAbsoluteScalar()       
        tubeFilter.SetNumberOfSides(20)
        tubeFilter.Update()         
         
        # Setup actor and mapper
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputConnection(tubeFilter.GetOutputPort())
        mapper.SetLookupTable(lut)        
         
        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
        actor.GetProperty().SetColor(1.0,0.0,0.0)        
        self.ren.AddActor(actor)        
        
    def findCellFor(self,point):
        locator = vtk.vtkCellLocator()
        locator.SetDataSet(self.grid)
        locator.BuildLocator()
        cellId = locator.FindCell(point)
        return cellId        
        
        
    def drawAxes(self):
        axes = vtk.vtkAxesActor()
        self.ren.AddActor(axes)

        
    def showWindow(self):
        self.renWin.AddRenderer(self.ren)
        self.iRen.Initialize()
        self.renWin.Render()
        self.iRen.Start()
        
        
        
def distance(point1, point2):
    return sqrt((point2[0] - point1[0])**2 + (point2[1] - point1[1])**2 + (point2[2] - point1[2])**2)
        
def dofToNode(dofNr, dofsPerNode):
    direction = dofNr%dofsPerNode
    nodeNr = (dofNr-direction)/dofsPerNode
    if direction > 0:
        nodeNr = nodeNr + 1
    return nodeNr, direction
    
def drange(start, stop, step):
    r = start
    while r < stop:
        yield r
        r += step    
        
if __name__ == '__main__':
    vtkWriter = VTKWriter()
    vtkWriter.createInput()
    vtkWriter.solveProblem()
    vtkWriter.displayBodies()
#    vtkWriter.drawFlowLines()
    vtkWriter.drawForceArrows()
    vtkWriter.drawBoundary()
#    vtkWriter.drawGaussPoints()
#    vtkWriter.drawAxes()
#    vtkWriter.cutPlane()
    vtkWriter.drawFlowLines()
    
    vtkWriter.showWindow()
    print vtk.vtkVersion.GetVTKSourceVersion()
