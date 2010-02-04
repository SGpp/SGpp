import vtk
import numpy
from vtk.util.colors import tomato, banana
from storage import Storage


math = vtk.vtkMath()
storage = Storage()

nr =  2

points = vtk.vtkPoints()
points.Allocate(nr,nr)

vectors = vtk.vtkFloatArray()
vectors.SetNumberOfComponents(3)

#values = storage.extractGrid(100, [0.0,0.5])

#for i in range(1, len(values), 2):
 #   points.InsertNextPoint((i-1)%nr, (i-1) - (i-1)%nr * nr , 0)
  #  vectors.InsertNextTuple3(values[i-1], values[i],0) 


points.InsertNextPoint(0, 0, 0)
points.InsertNextPoint(1, 0, 0)
points.InsertNextPoint(0, 1, 0)
points.InsertNextPoint(1, 1, 0)
vectors.InsertNextTuple3(5, 5,0) 
vectors.InsertNextTuple3(5, -5,0) 
vectors.InsertNextTuple3(-5, 5,0) 
vectors.InsertNextTuple3(-4, -7,0) 

print vectors

#grid  = vtk.vtkUniformGrid()
grid = vtk.vtkStructuredGrid()
grid.SetDimensions([nr, nr, 0])
grid.SetPoints(points)
#grid.SetOrigin(0,0,0)
#grid.SetSpacing(0.01,0.01,0.0)
#grid.SetExtent(0,99,0,99,0,0)
grid.GetPointData().SetVectors(vectors)

print grid

append = vtk.vtkAppendFilter()
append.AddInput(grid)

# Use sphere as glyph source.
balls = vtk.vtkArrowSource()
balls.SetTipLength(0.0035)
balls.SetTipRadius(0.05)
balls.SetTipResolution(6)
balls.SetShaftResolution(6)
balls.SetShaftRadius(0.0003)

glyph = vtk.vtkGlyph3D()
glyph.SetInputConnection(append.GetOutputPort())
glyph.SetSourceConnection(balls.GetOutputPort())
glyph.SetVectorModeToUseVector()
glyph.OrientOn()
glyph.SetScaleModeToScaleByVector()
glyph.SetScaleFactor(3)

glyphMapper = vtk.vtkPolyDataMapper()
glyphMapper.SetInputConnection(glyph.GetOutputPort())

glyph = vtk.vtkActor()
glyph.SetMapper(glyphMapper)
#glyph.GetProperty().SetDiffuseColor(tomato)
#glyph.GetProperty().SetSpecular(.3)
#glyph.GetProperty().SetSpecularPower(30)

# bounding box
#outline = vtk.vtkOutlineFilter()
#outline.SetInput( grid2.GetOutput() )
#outlineMapper = vtk.vtkPolyDataMapper()
#outlineMapper.SetInput( outline.GetOutput() )
#outlineActor = vtk.vtkActor()
#outlineActor.SetMapper( outlineMapper )
#outlineActor.GetProperty().SetColor(0.0,0.0,1.0)

ren = vtk.vtkRenderer()
ren.SetBackground(.8, .8, .8)
renWin = vtk.vtkRenderWindow()
renWin.SetSize(800, 800)
renWin.AddRenderer( ren )

iren = vtk.vtkRenderWindowInteractor()
iren.SetRenderWindow(renWin)

#ren.AddActor(outlineActor)
ren.AddActor(glyph)

renWin.Render()

iren.Initialize()
iren.Start()
