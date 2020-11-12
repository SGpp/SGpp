# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

import random
import numpy as np

def normiere_liste(liste, dimensions):
   B = np.array(liste)
   for dimension in range(0, dimensions):
      teil = B[:, dimension]
      minimum = min(teil)
      maximum = max(teil)
      for i in range(0, len(teil)):
         teil[i] = (teil[i]-minimum) / (maximum-minimum)
         teil[i] = (teil[i]+0.1)*0.8
      B[:, dimension] = teil
   return B

def generate_dataset(dimensions, clusters, setsize, abweichung, rauschensize):
   #generate centers
   centers = []
   for cluster in range(0, clusters):
      abstand = 0
      center = []
      counter = 0
      continue_search = 1
      while continue_search:
         center[:] = []
         if counter == 50:
            abweichung = abweichung * 0.9
            counter = 0
         for dimension in range(0, dimensions):
            center.append(random.random())
         #check point
         continue_search = 0
         for c in centers:
            abstand = 0
            for dimension in range(0, dimensions):
               abstand = abstand + (center[dimension] - c[dimension])**2
            abstand = abstand ** 0.5
            if abstand < 32 * abweichung:
               continue_search = 1
               break
         counter = counter +1
         if not centers:
            break
      centers.append(center)
   #generate datapoints
   dataset = []
   cluster_ret = []
   currentsize = setsize
   for cluster in range(0, clusters):
      if currentsize > 0:
         for i in range(0, min(currentsize, setsize / clusters)):
            point = []
            for dimension in range(0, dimensions):
               point.append(random.gauss(centers[cluster][dimension], abweichung))
            dataset.append(point)
            cluster_ret.append(cluster+1)
            currentsize = currentsize - 1
         #currentsize = currentsize - setsize / clusters
   while currentsize is not 0:
      point = []
      for dimension in range(0, dimensions):
         point.append(random.gauss(centers[cluster][dimension], abweichung))
      dataset.append(point)
      currentsize = currentsize - 1
      cluster_ret.append(clusters)
   for i in range(0, rauschensize):
      point = []
      for dimension in range(0, dimensions):
         point.append(random.uniform(0.05,0.95))
      dataset.append(point)
      cluster_ret.append(-1)
   dataset = normiere_liste(dataset, dimensions)
   cluster_npret = np.array(cluster_ret)
   # print dataset
   return dataset, cluster_npret

# dimensions, clusters, setsize, abweichung, rauschensize
clusters=10
start_size=20000
dim=10
dataset_size=start_size
for i in range(0, 10):
   N = pow(2, i)
   print("creating" + "datasets/weak_N" + str(N) + "_c" + str(clusters) + "_size" + str(dataset_size) + "_dim" + str(dim) + ".arff")
   dataset1, Y1 = generate_dataset(dim, clusters, dataset_size, 0.1, 0)
   np.savetxt("datasets/weak_N" + str(N) + "_c" + str(clusters) + "_size" + str(dataset_size) + "_dim" + str(dim) + ".arff", dataset1)
   np.savetxt("datasets/weak_N" + str(N) + "_c" + str(clusters) + "_size" + str(dataset_size) + "_dim" + str(dim) + "_erg.txt", Y1)
   dataset_size*=2
