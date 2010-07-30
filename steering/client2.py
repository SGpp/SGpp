import socket
import sys


class SocketClient:
  
  HOST = 'localhost'    # The remote host
  PORT = 9000         # The same port as used by the server
  
  sendSocket = None
  receiveSocket = None
  
  EOF = "EOF"

  header = """
  <?xml version="1.0" encoding="UTF-8" standalone="yes" ?>
  <!DOCTYPE boost_serialization>
  <boost_serialization signature="serialization::archive" version="5">
  <cell class_id="0" tracking_level="0" version="0">
     <type>12</type>
     <data class_id="9" tracking_level="1" version="0" object_id="_0">
       <parent class_id="10" tracking_level="1" version="0" object_id="_1"></parent>
       <re>2</re>
       <time>0.10000000000000001</time> 
        <velocities> """

  footer = """
     </velocities>
    </data>
   </cell>
   EOF
   """

  termination = """
  <?xml version="1.0" encoding="UTF-8" standalone="yes" ?>
  <!DOCTYPE boost_serialization>
  <boost_serialization signature="serialization::archive" version="5">
  <cell class_id="0" tracking_level="0" version="0">
  <type>9</type>
  <data class_id="1" tracking_level="1" version="0" object_id="_0">
      <parent class_id="9" tracking_level="1" version="0" object_id="_1"></parent>
  </data>
  </cell>
  EOF
  """  
  
  def __init__(self):
    self.connect()

  def sendCurrentGrid(self, values):
    serializedData = self.header + "\n<count>" + str(len(values)) + "</count>\n"
    serializedData = serializedData + "<item_version>0</item_version>\n"
 
    #for i in xrange(len(values)):
      #serializedData = serializedData + "<item>" + str(values[i]) + "</item>\n"
  
    #serializedData = serializedData + self.footer
    text = values[2:]
    self.send(values + "\n")

  def sendTerminationMessage(self):
    self.send(self.termination)
    return

  def send(self, dataString):
    totalsent = 0
    while totalsent < len(dataString):
      print "sending"
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
    #self.receiveSocket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    
    self.sendSocket.connect((self.HOST, self.PORT))
    #self.receiveSocket.connect((self.HOST, self.PORT + 1))
  
  def disconnect(self):
    if self.sendSocket is None:
      return
    else: 
      self.sendSocket.close()

 

