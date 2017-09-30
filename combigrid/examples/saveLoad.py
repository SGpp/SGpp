import os
import numpy as np


def save_array(array,filename,dirname):


    if not os.path.exists(dirname):
        os.makedirs(dirname)
    np.savetxt(dirname+"/"+filename,array)

def save_namelist(list,filename,dirname):

    if not os.path.exists(dirname):
        os.makedirs(dirname)

    with open(dirname+"/"+filename, 'w') as file_handler:
        for item in list:
            file_handler.write(str(item)+"\n")

        file_handler.close()

def load_array(filename,dirname):

   return np.loadtxt(dirname+"/"+filename)


def load_all(directory,ending):
    for filename in os.listdir(directory):
        if filename.endswith(ending):
            print(filename)

def load_list(filename,dirname):
    with open(dirname + "/" + filename, 'r') as file_handler:

        return(file_handler.read().splitlines(False))







"""
x=np.zeros((10,10))
y=np.array([3,2,23,5,767,45])

save_array(x,"testiooo",dirname="f1")
save_array(y,"testiooo2",dirname="f1")

load_all("f1")

"""

save_namelist(["test","blub"],"labels.txt","testiii")
