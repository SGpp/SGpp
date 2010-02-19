import vtk
from Tkinter import *
from storage import Storage

from vtk import *
from vtk.tk.vtkTkRenderWidget import vtkTkRenderWidget


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
    self.dimensionSize = 100
    self.once = 0 
    self.grid = vtk.vtkStructuredGrid()
    self.createGUI()

  
  def createGUI(self): 
    # Make a root window
    self.TK = Tk() 

    # Add a vtkTkRenderWidget
    self.renderWidget = vtkTkRenderWidget(self.TK,width=800,height=800)
    self.renderWidget.grid(row=0, column=2)

    # Get the render window from the widget
    scale = Scale ( self.TK, label='T', orient='vertical', from_=1, to=100, command=self.changeTime)
    scale.grid(row=0,column=0, sticky=N+S)
    scale = Scale ( self.TK, label='Re', orient='vertical', from_=1, to=100, command=self.changeRe)
    scale.grid(row=0,column=1, sticky=N+S)

    self.renWin = self.renderWidget.GetRenderWindow()

    # Next, do the VTK stuff
    ren = vtkRenderer()
    ren.SetBackground(.8, .8, .8)
    self.renWin.AddRenderer(ren)

    self.loadGrid()
    self.created = True
  
    append = vtk.vtkAppendFilter()
    append.AddInput(self.grid)

    balls = vtk.vtkArrowSource()

    glyph = vtk.vtkGlyph3D()
    glyph.SetInputConnection(append.GetOutputPort())
    glyph.SetSourceConnection(balls.GetOutputPort())
    glyph.SetVectorModeToUseVector()
    glyph.SetScaleModeToScaleByVector()
    glyph.SetScaleFactor(10)

    glyphMapper = vtk.vtkPolyDataMapper()
    glyphMapper.SetInputConnection(glyph.GetOutputPort())

    glyph = vtk.vtkActor()
    glyph.SetMapper(glyphMapper)
  
    rake = vtk.vtkLineSource()
    rake.SetPoint1(15, -5, 0)
    rake.SetPoint2(15, 5, 0)
    rake.SetResolution(5)
    rakeMapper = vtk.vtkPolyDataMapper()
    rakeMapper.SetInputConnection(rake.GetOutputPort())
    rakeActor = vtk.vtkActor()

    rk4 = vtk.vtkRungeKutta4()
    streamer = vtk.vtkStreamTracer()
    streamer.SetInputConnection(append.GetOutputPort())
    #streamer.SetSource(rake)
    streamer.SetStartPosition(50, 50, 0)
    #streamer.SetMaximumPropagationTime(100)
    #streamer.SetIntegrationStepLength(.2)
    #streamer.SetStepLength(.001)
    #streamer.SetNumberOfThreads(1)
    #streamer.SetIntegrationDirectionToForward()
    #streamer.VorticityOn()
    streamer.SetIntegrator(rk4)
    #rf = vtk.vtkRibbonFilter()
    #rf.SetInputConnection(streamer.GetOutputPort())
    #rf.SetWidth(0.1)
    #rf.SetWidthFactor(5)
    streamMapper = vtk.vtkPolyDataMapper()
    streamMapper.SetInputConnection(streamer.GetOutputPort())
    streamMapper.SetScalarRange(append.GetOutput().GetScalarRange())
    streamline = vtk.vtkActor()
    streamline.SetMapper(streamMapper)
    streamline.VisibilityOn()

    outline = vtk.vtkOutlineFilter()
    outline.SetInputConnection( append.GetOutputPort() )
    outlineMapper = vtk.vtkPolyDataMapper()
    outlineMapper.SetInput( outline.GetOutput() )
    outlineActor = vtk.vtkActor()
    outlineActor.SetMapper( outlineMapper )
    outlineActor.GetProperty().SetColor(0.0,0.0,1.0)

    ren.AddActor(outlineActor)
    ren.AddActor(glyph)
    ren.AddActor(rakeActor)
    ren.AddActor(streamline)
 
    #lineWidget.SetInteractor(self.renWin.GetInteractor()) 
    #lineWidget.AddObserver("StartInteractionEvent", self.generateStreamlines)
    #lineWidget.AddObserver("InteractionEvent", self.generateStreamlines)
    #lineWidget.AddObserver("EndInteractionEvent", self.generateStreamlines)
 
    #button = Button(text="Quit",command=quit)
    #button.pack(expand='true',fill='x')
 
    # start up the event loop
    self.TK.mainloop()

  def constructActors(self):
    return 
   
  def loadGrid(self):
    nr = self.dimensionSize

    if not self.created:
      points = vtk.vtkPoints()
      points.Allocate(nr,nr)

    vectors = vtk.vtkFloatArray()
    vectors.SetNumberOfComponents(3)

    values = self.storage.extractGrid(self.dimensionSize, self.currParam)

    for i in range(0, 20000, 2):
   
      if not self.created:
        points.InsertNextPoint( i / (2 * nr) , (i/2) % nr, 0)
  
      vectors.InsertNextTuple3(values[i], values[i+1],0) 
      #vectors.InsertNextTuple3(i%3, 1, 0) 
  
    if not self.created:
      self.grid.SetPoints(points)

    self.grid.GetPointData().SetVectors(vectors)

    return

  def generateStreamlines(obj, event):
    print "bla bla"
    obj.GetPolyData(seeds)
    self.stream.VisibilityOn()
    self.renWin.Render()

  def changeTime(self, value):
    self.currParam[0] = float(value)/100.0
    self.loadGrid()
    self.renWin.Render()
  
  def changeRe(self, value):
    self.currParam[1] = float(value)/100.0
    self.loadGrid()
    self.renWin.Render()


if __name__ == '__main__':
    S = SampleViewer()
