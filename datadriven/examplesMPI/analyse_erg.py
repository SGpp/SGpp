from builtins import str
import numpy as np
import sys

#file = raw_input("Filename: ")
#referenz = raw_input("Reference filename: ")
if len(sys.argv) != 3:
    print("Missing filenames for clustering results and optimal results!")
    print("USAGE: python2 analyse_erg.py <clustering result file> <reference file with (optimal) results>")
    sys.exit()
file = sys.argv[1]
referenz = sys.argv[2]
file1 = np.loadtxt(file)
file_ref = np.loadtxt(referenz)

found_clusters = []
clustercount = []
for i in file1:
    if i not in found_clusters:
        found_clusters.append(i)
        clustercount.append(1)
    else:
        counter = 0
        for c in found_clusters:
            if int(c) is int(i):
                clustercount[counter]=clustercount[counter]+1
                break
            counter=counter+1
print("Clusteranalyse: ")
counter = 0
confirmed_hits = []
gesamthits=0
for c in found_clusters:
    if int(c) is 0:
        print("Entfernte Datenpunkte: " + str(clustercount[counter]))
    else:
        print(str(found_clusters[counter]) + " mit " + str(clustercount[counter]))
    counter = counter + 1
    confirmed_hits.append(0)
counter=0
for i in file1:
    tmp=float(file_ref[counter])
    tocomp = int(tmp)
    if tocomp is -1:
        tocomp = 0
    if int(i) is tocomp:
        gesamthits = gesamthits +1
        clustercounter = 0
        for c in found_clusters:
            if int(c) is int(i):
                confirmed_hits[clustercounter]=confirmed_hits[clustercounter]+1
                break
            clustercounter = clustercounter + 1
    counter = counter + 1
print("Trefferanalyse:")
print("Gesamte Treffer: " + str(gesamthits))
counter = 0
for c in found_clusters:
    if int(c) is 0:
        print("Rauschen (zu entfernende Punkte) mit " + str(confirmed_hits[counter])+" Treffer")
    else:
        print(str(found_clusters[counter]) + " mit " + str(confirmed_hits[counter])+" Treffer")
    counter = counter + 1
