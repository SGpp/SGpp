from math import *


def get1Dstretching(level,lbdry,rbdry): 
	nr_pts = level**2 + 1
	stretchingList=[]
	diff_ = (rbdry-lbdry)/2.0
	sum_ = (rbdry + lbdry)/2.0
	factor = pi/(nr_pts - 1.0)
	stretchingList.append(lbdry)
	for i in range(1,nr_pts-1,2):
		stretchingList.append(sum_ - diff_*cos(i*factor))
		stretchingList.append(sum_ - diff_*cos((i+1)*factor))
	stretchingList.append(rbdry) 
	# now return the result 
	return stretchingList		


def get_coeffs(level):
	omegas = [];
	N = 2**level	
	#straight-forward N^2 complexity implenentation 
	#of coeffs computation with clenshaw curtis... 
	for n in range(0,N+1):
		w_n = 2 
		factor = 2*pi*n/N
		for k in range(1,N/2):
			w_n +=4/(1-4*(k**2))*cos(factor*k)
		w_n +=2/(1+N**2)
		w_n /= N #finally divide the result by N;
		omegas.append(w_n) 

	return omegas 


def integrate_clenshaw(f,level,a,b):
	I = 0
	N = 2**level + 1
	x = get1Dstretching(level,a,b)
	omegas = get_coeffs(level)
	print "coeffs of size: ", len(omegas) 
	print 	
	print "points of size: ", len(x) 
	print 

	for n in range(0,N):
		I += f(x[n])*omegas[n]	






