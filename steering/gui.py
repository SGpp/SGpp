import vtk
from Tkinter import *
from storage import Storage

from vtk import *
from vtk.tk.vtkTkRenderWidget import vtkTkRenderWidget

import time

class SampleViewer:

  storage = None
  TK = None
  currParam = None
  renderWidget = None
  dimensionSize = None
  grid = None 
  ren = None
  renWin = None
  once = None
  values = None
  seeds = None
  stream = None
  created = False

  def __init__ ( self ):
    self.storage = Storage()
    self.currParam = [0.5, 0.5]
    self.dimensionSize = 10
    self.once = 0 
    self.grid = vtk.vtkStructuredGrid()
    self.createGUI()

  
  def createGUI(self): 
    # Make a root window
    self.TK = Tk() 

    #info = Label(self.TK, text = "Storage Info: " + str(self.storage.grid.getStorage().size()) + " points.")
    info = Label(self.TK, text = "Storage Info: " + " points.")
    info.grid(row=0,column=2)  

    # Add a vtkTkRenderWidget
    self.renderWidget = vtkTkRenderWidget(self.TK,width=1200,height=800)
    self.renderWidget.grid(row=1, column=2)

    # Get the render window from the widget
    scale = Scale ( self.TK, label='T', orient='vertical', from_=1, to=100, command=self.changeTime)
    scale.grid(row=1,column=0,sticky=N+S)
    scale = Scale ( self.TK, label='Re', orient='vertical', from_=1, to=100, command=self.changeRe)
    scale.grid(row=1,column=1,sticky=N+S)

    var = IntVar()
    check = Checkbutton(self.TK, text="Colored Vectors", variable=var)
    check.grid(row=3,column=2)    

    self.renWin = self.renderWidget.GetRenderWindow()

    # Next, do the VTK stuff
    ren = vtkRenderer()
    ren.SetBackground(.8, .8, .8)
    self.renWin.AddRenderer(ren)

    self.loadGrid()
    self.created = False
  
    append = vtk.vtkAppendFilter()
    append.AddInput(self.grid)

    balls = vtk.vtkArrowSource()

    glyph = vtk.vtkGlyph3D()
    glyph.SetInputConnection(append.GetOutputPort())
    glyph.SetSourceConnection(balls.GetOutputPort())
    glyph.SetVectorModeToUseVector()
    glyph.SetScaleModeToScaleByVector()
    glyph.SetScaleFactor(10)
    glyph.SetColorModeToColorByVector()

    glyphMapper = vtk.vtkPolyDataMapper()
    glyphMapper.SetInputConnection(glyph.GetOutputPort())

    glyph = vtk.vtkActor()
    glyph.SetMapper(glyphMapper)

    outline = vtk.vtkOutlineFilter()
    outline.SetInputConnection(append.GetOutputPort())
    outlineMapper = vtk.vtkPolyDataMapper()
    outlineMapper.SetInput( outline.GetOutput() )
    outlineActor = vtk.vtkActor()
    outlineActor.SetMapper( outlineMapper )
    outlineActor.GetProperty().SetColor(0.0,0.0,1.0)

    seeds = vtk.vtkLineSource()
    seeds.SetPoint1(5, 5, 5)
    seeds.SetPoint2(7, 7, 7)
    seeds.SetResolution(1)

    integ = vtk.vtkRungeKutta4()
    sl = vtk.vtkStreamLine()
    sl.SetInput(self.grid)
    sl.SetStartPosition(5, 5, 5)
    sl.SetMaximumPropagationTime(500)
    sl.SetStepLength(0.0005)
    sl.SetIntegrationStepLength(0.00005)
    sl.SetIntegrationDirectionToIntegrateBothDirections()
    sl.SetIntegrator(integ)
 
    tube = vtk.vtkTubeFilter()
    tube.SetInputConnection(sl.GetOutputPort())
    tube.SetRadius(0.1)

    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputConnection(tube.GetOutputPort())

    actor = vtk.vtkActor()
    actor.SetMapper(mapper)

    ren.AddActor(outlineActor)
    ren.AddActor(glyph)
    #ren.AddActor(actor)

    # start up the event loop
    self.TK.mainloop()
  
  def loadGrid(self):
    print "Loading new grid."

    nr = self.dimensionSize

    if not self.created:
      points = vtk.vtkPoints()
      points.Allocate(nr,nr)

    vectors = vtk.vtkFloatArray()
    vectors.SetNumberOfComponents(3)

    x = time.time()
    values = self.storage.extractGrid(nr, [0.5,0.5])#self.currParam)
    print "Extract time ", time.time() - x 

    dim = 3
  
    for i in range(0, dim * (nr**3), dim):
      if not self.created:
        points.InsertNextPoint( (i/(3 * nr)) % nr , (i/3) % (nr), i / (3 * nr*nr) )
      vectors.InsertNextTuple3(values[i], values[i+1], values[i+2]) 
  
    if not self.created:
      self.grid.SetPoints(points)

    self.grid.GetPointData().SetVectors(vectors)

    return

  def changeTime(self, value):
    self.currParam[0] = float(value)/100.0
    x = time.time()
    self.loadGrid()
    print "Load Grid Time ChangeT ", time.time() - x
    self.renWin.Render()
  
  def changeRe(self, value):
    self.currParam[1] = float(value)/100.0
    x = time.time()
    self.loadGrid()
    print "Load Grid Time ChangeRe ", time.time() - x
    self.renWin.Render()


if __name__ == '__main__':
    S = SampleViewer()
