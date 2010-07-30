import socket
import sys
from xml.dom.minidom import parseString, getDOMImplementation
from storage import Storage


class SocketClient:
  
  storage = None

  HOST = 'atsccs24.informatik.tu-muenchen.de'  # The remote host
  SERVER_PORT = 9000         # The same port as used by the server
  GUI_PORT = 14007 
 
  serverSocket = None
  guiSocket = None
  
  MESS_SCOPE_QUERY  = "scopeQuery"
  MESS_SCOPE_ANSWER = "scopeAnswer"
  MESS_DATA_QUERY   = "dataQuery"
  MESS_DATA_ANSWER  = "dataAnswer"
 
  EOF = "eof"

  
  def __init__(self, storage):
    self.connect()
    self.storage = storage

  def doCommunication(self):
    message = self.recv_end()
    #message = getDOMImplementation().createDocument(None, self.MESS_SCOPE_QUERY , None).toxml() 
    print "Received: " + message
    xml = parseString(message)    
      
    if xml.documentElement.tagName == self.MESS_SCOPE_QUERY:
      self.respondToScopeQuery()
    elif xml.documentElement.tagName == self.MESS_DATA_QUERY:
      self.respondToDataQuery(xml)
    else:
      print "ERROR: Unrecognized Message Tag."

    return


  def respondToScopeQuery(self):
    answer = getDOMImplementation().createDocument(None, self.MESS_SCOPE_ANSWER, None)
    
    scopeID = answer.createElement("scopeID")
    scopeID.setAttribute("value", "0")   

    dimension = answer.createElement("dimension")
    dimension.setAttribute("value", "2")   

    bottomRefPoint = answer.createElement("bottomReferencePoint")
    bottomRefPoint.setAttribute("x0", "0")   
    bottomRefPoint.setAttribute("x1", "0")   
    
    boundingBox = answer.createElement("boundingBox")
    boundingBox.setAttribute("h0", "100")   
    boundingBox.setAttribute("h1", "100")   
    
    answer.documentElement.appendChild(scopeID)
    answer.documentElement.appendChild(dimension)
    answer.documentElement.appendChild(bottomRefPoint)
    answer.documentElement.appendChild(boundingBox)
    
    self.sendServer(answer.documentElement.toprettyxml() + "\n eof")    
    self.sendGUI(answer.documentElement.toprettyxml() + "\n eof")    
 
  def respondToDataQuery(self, dataQuery):
    
    answer = getDOMImplementation().createDocument(None, self.MESS_DATA_ANSWER, None)
	
    scopeID = answer.createElement("scopeID")
    scopeID.setAttribute("value", "0")   

    dimension = answer.createElement("meshDimensions")
    dimension.setAttribute("h0", "64")   
    dimension.setAttribute("h1", "64")   
    
    bottomRefPoint = answer.createElement("bottomReferencePoint")
    bottomRefPoint.setAttribute("x0", "0")   
    bottomRefPoint.setAttribute("x1", "0")   
   
    boundingBox = answer.createElement("boundingBox")
    boundingBox.setAttribute("h0", "1")   
    boundingBox.setAttribute("h1", "1") 

    data = answer.createElement("data")
    
    values = self.storage.extractGrid(64, [0.0,0.5])   

    myfile = open("myfile", "w")
 
    #for i in range(len(values)):
    #  myfile.write(str(values[i]) + " " + str(values[i+1]) + " 0\n" )
    #  i +=2     

    for i in range(len(values)):
      subdata = answer.createElement("subdata")
      subdata.setAttribute("value", str(values[i]))
      data.appendChild(subdata)
    
    myfile.close()

    answer.documentElement.appendChild(scopeID)
    answer.documentElement.appendChild(dimension)
    answer.documentElement.appendChild(bottomRefPoint)
    answer.documentElement.appendChild(boundingBox)
    answer.documentElement.appendChild(data)

    self.sendServer(answer.documentElement.toprettyxml() + "\n eof")    


  def sendServer(self, dataString):
    totalsent = 0
    while totalsent < len(dataString):
      sent = self.serverSocket.send(dataString[totalsent:])
      if sent == 0:
        raise RuntimeError, "socket connection broken"
      totalsent = totalsent + sent


  def sendGUI(self, dataString):
    totalsent = 0
    while totalsent < len(dataString):
      sent = self.guiSocket.send(dataString[totalsent:])
      if sent == 0:
        raise RuntimeError, "socket connection broken"
      totalsent = totalsent + sent

 
  def recv_end(self):
    End = self.EOF
    total_data=[]
    data=''
    got_end=False
    while True:
      data=self.serverSocket.recv(1024)
      if not data: break
      if End in data:
        total_data.append(data[:data.find(End)])
        got_end=True
        break
      total_data.append(data)
       
      if len(total_data)>1:
        #check if end_of_data was split
        last_pair=total_data[-2]+total_data[-1]
        if End in last_pair:
          total_data[-2]=last_pair[:last_pair.find(End)]
          total_data.pop()
          got_end=True
          break
    return ''.join(total_data)


  def connect(self):
    self.serverSocket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    self.guiSocket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    
    self.serverSocket.connect((self.HOST, self.SERVER_PORT))
    self.guiSocket.connect((self.HOST, self.GUI_PORT))
  
#  def disconnect(self):
 #   if self.sendSocket is None:
  #    return
   # else: 
    #  self.sendSocket.close()
