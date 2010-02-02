import socket
import sys
from xml.dom.minidom import parseString, getDOMImplementation
from storage import Storage


class SocketClient:
  
  storage = None

  HOST = 'atsccs24'    # The remote host
  PORT = 9000         # The same port as used by the server
  
  sendSocket = None
  receiveSocket = None
  
  MESS_SCOPE_QUERY  = "scopeQuery"
  MESS_SCOPE_ANSWER = "scopeAnswer"
  MESS_DATA_QUERY   = "dataQuery"
  MESS_DATA_ANSWER  = "dataAnswer"
 
  EOF = "EOF"

  
  def __init__(self, storage):
    #self.connect()
    self.storage = storage


  def doCommunication(self):
    #message = recv_end()
    message = getDOMImplementation().createDocument(None, self.MESS_DATA_QUERY , None).toxml() 
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
    dimension.setAttribute("value", "0")   

    bottomRefPoint = answer.createElement("bottomReferencePoint")
    bottomRefPoint.setAttribute("x1", "0")   

    boundingBox = answer.createElement("boundingBox")
    boundingBox.setAttribute("h1", "0")   

    answer.documentElement.appendChild(scopeID)
    answer.documentElement.appendChild(dimension)
    answer.documentElement.appendChild(bottomRefPoint)
    answer.documentElement.appendChild(boundingBox)
    
    print answer.toprettyxml()   
    #self.send(answer.toxml())    

 
  def respondToDataQuery(self, dataQuery):
    
    answer = getDOMImplementation().createDocument(None, self.MESS_DATA_ANSWER, None)
	
    scopeID = answer.createElement("scopeID")
    scopeID.setAttribute("value", "0")   

    dimension = answer.createElement("dimension")
    dimension.setAttribute("value", "0")   

    bottomRefPoint = answer.createElement("bottomReferencePoint")
    bottomRefPoint.setAttribute("x1", "0")   

    boundingBox = answer.createElement("boundingBox")
    boundingBox.setAttribute("h1", "0")   

    data = answer.createElement("data")
    
    values = self.storage.extractGrid(100, [0.5,0.5])   

    for i in range(len(values)):
      subdata = answer.createElement("subdata")
      subdata.setAttribute("value", str(values[i]))
      data.appendChild(subdata)
    
    answer.documentElement.appendChild(scopeID)
    answer.documentElement.appendChild(dimension)
    answer.documentElement.appendChild(bottomRefPoint)
    answer.documentElement.appendChild(boundingBox)
    answer.documentElement.appendChild(data)

    print answer.toprettyxml()   
    #self.send(answer.toxml())    


  def send(self, dataString):
    totalsent = 0
    while totalsent < len(dataString):
      sent = self.sendSocket.send(dataString[totalsent:])
      if sent == 0:
        raise RuntimeError, "socket connection broken"
      totalsent = totalsent + sent

  def receiveParametersMessage(self):
    xmlString = self.recv_end() 

 
  def recv_end(self):
    End = self.EOF
    total_data=[]
    data=''
    got_end=False
    while True:
      data=self.receiveSocket.recv(1024)
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
    self.sendSocket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    self.receiveSocket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    
    self.sendSocket.connect((self.HOST, self.PORT))
    self.receiveSocket.connect((self.HOST, self.PORT + 1))

  
  def disconnect(self):
    if self.sendSocket is None:
      return
    else: 
      self.sendSocket.close()
