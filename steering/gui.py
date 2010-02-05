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
  
    append = vtk.vtkAppendFilter()
    append.AddInput(self.grid)

    balls = vtk.vtkArrowSource()

    glyph = vtk.vtkGlyph3D()
    glyph.SetInputConnection(append.GetOutputPort())
    glyph.SetSourceConnection(balls.GetOutputPort())
    glyph.SetVectorModeToUseVector()
    glyph.SetScaleModeToScaleByVector()
    glyph.SetScaleFactor(50)

    glyphMapper = vtk.vtkPolyDataMapper()
    glyphMapper.SetInputConnection(glyph.GetOutputPort())

    glyph = vtk.vtkActor()
    glyph.SetMapper(glyphMapper)

    # bounding box
    outline = vtk.vtkOutlineFilter()
    outline.SetInputConnection( append.GetOutputPort() )
    outlineMapper = vtk.vtkPolyDataMapper()
    outlineMapper.SetInput( outline.GetOutput() )
    outlineActor = vtk.vtkActor()
    outlineActor.SetMapper( outlineMapper )
    outlineActor.GetProperty().SetColor(0.0,0.0,1.0)

    ren.AddActor(outlineActor)
    ren.AddActor(glyph)

    #button = Button(text="Quit",command=quit)
    #button.pack(expand='true',fill='x')
 
    # start up the event loop
    self.TK.mainloop()

  def constructActors(self):
    return 
   
  def loadGrid(self):
    nr = self.dimensionSize

    points = vtk.vtkPoints()
    points.Allocate(nr,nr)

    vectors = vtk.vtkFloatArray()
    vectors.SetNumberOfComponents(3)

    values = self.storage.extractGrid(100, self.currParam)

    print len(values)

    for i in range(0, len(values), 2):
      points.InsertNextPoint( i / (2 * nr) , (i/2) % nr, 0)
      vectors.InsertNextTuple3(-1 * values[i], values[i+1],0) 

    self.grid.SetPoints(points)
    self.grid.GetPointData().SetVectors(vectors)

    return

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
