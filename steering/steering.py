from pysgpp import *
import re
import tools
import pickle
from client import SocketClient
from storage import Storage

class Steering:
  
  socketClient = None
  storage = None

  def __init__(self):
    return
 
  def run(self):
    print "Create Storage."
    self.storage = Storage()

    print "Start Client."
    self.socketClient = SocketClient(self.storage)
 
    print "Enter Communication Loop."
    while 1:
      self.socketClient.doCommunication()
    #  params = self.socketClient.recv_end()
    #  print params
      #data = self.extractGrid(100, 0.5)
      #self.socketClient.sendCurrentGrid(data)
      #self.socketClient.sendTerminationMessage()   

    return
# end Steering class

### Set up and start of the program
r = Steering()
r.run()
