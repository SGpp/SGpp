#!/usr/bin/python

# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at 
# sgpp.sparsegrids.org

##
# Helper tool, normalize dataset using PCA algorithm
# Data normalisation using principal components
# 


#from create_rotation_matrix import *
from numpy import *
import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg.linalg import eig, inv
import logging
import logging.config


## normalization on [0,1] interval
# @param mat: matrix points row-wise
def norm(mat):
	NDIM = shape(mat)[1] # number of dimensions
	for i in xrange(NDIM):
		m = min(mat[:,i])
		M = max(mat[:,i])
		mat[:,i] = (mat[:,i] - m)/(M-m)
	return mat


## remove outliers, where the points deviate on more than koef times standard
# deviation from the mean
# @param mat: matrix points row-wise
# @param koef: koefficient to determine the outliers
# @param target: list with target values for the points
def remove_outliers(mat, koef, target=None):
	#mat = mat.T
	(size, dim) = shape(mat)
	m = mean(mat,0)
	s = std(mat,0)
	result = m
	target_result = []
	for i in xrange(size):
		point = mat[i,:]
		if sum(point > m - koef*s) == dim and sum(point < m + koef*s) == dim:
			#target = targets[i]
			#result = concatenate((result, concatenate((point,[[target]]),axis=1)), axis=0)
			result = concatenate((result, point), axis=0)
			if target != None:
				target_result.append(target[i])
	result = matrix(result)
	return (result[1:, :], target_result)
	

## leave only the data points with coordinated greater than zero
# is convinient for some problems, i.e. photometric redshift
def clean_data(data):
	data = data[np.min(data[:,1:6],axis=1)>0,]
	return data
	

def create_logger():
	# create logger
	logger = logging.getLogger()
	logger.setLevel(logging.ERROR)
	logger.addHandler(logging.StreamHandler())
	return logger
	
	
if __name__ == '__main__':
	logger = create_logger()
	import optparse
	parser = optparse.OptionParser()
	parser.add_option("-d","--directory", action="store", type="string", 
			  dest="csv_dir", help="The input/output directory")
	parser.add_option("-i","--file-in", action="store", type="string", 
		dest="file_in", help="Input file")
	parser.add_option("-o","--file-out", action="store", type="string", 
        	dest="file_out", help="Output file")
	parser.add_option("-k","--koef", action="store", type="int", 
        	dest="koef", help="remove outliers, where the points deviate on " 
        	+ "more than --koef times standard deviation")
	parser.add_option("-v","--verbose", action="store_true", default=False, 
	      dest="verbose", help="Turns the verbose mode on.")
	parser.add_option("-c", "--target-column", dest="target_column", 
					  action="store", type="int", help = "Number of the column with target data")
	# parse options
	(options,args)=parser.parse_args()

	if options.csv_dir == None:
		parser.error("Parameter --directory is required")
	if options.file_in == None:
		parser.error("Parameter --file-in is required")
	#	if options.names=None:
	#        	parser.error("Parameter --name is required") 

	if options.verbose:
		logger.setLevel(logging.INFO)

	# specify the input/output directory
	#	csv_dir = '/home/perun/Arbeit/redshift/photometric_redshift/datasets/original_csv/'
	# specify the filenames
	#	files = ['DR7_cleaned_galaxies_MGS.csv']

		 
	filename = options.file_in
	logger.info( "Processing file %s" % filename )

	#read data from the file
	logger.info("Loading the data...")
	data = genfromtxt(options.csv_dir + filename, skiprows=1, delimiter=',')
	
	#clean data
	logger.info("Cleaning the data...")
	data = clean_data(data)
	if options.target_column != None and options.target_column < shape(data)[1]:
		target = data[:, options.target_column]
		if options.target_column != shape(data)[1]-1:
			data = concatenate((data[:, 0:options.target_column],
						    data[:, options.target_column+1:-1]), axis=0)
		else:
			data = data[:, 0:options.target_column]
	#subselect data dimensions of interest
	means = mean(data,axis=0)
	interest_data = data - means
	
	# transform data using PCA components
	logger.info("Transforming the data...")
	interest_data = matrix(interest_data).T # point are written columnwise
	C = cov(interest_data) # covariance matrix
	(u,V) = eig(C) # eigenvector of C in V columnwise
	invV = inv(V)
	interest_data_transformed = dot(invV, interest_data).T

	# uncomment if you want to use only the points within specified range
	logger.info("Normalizing the data...")
	if options.koef:
		(interest_data_transformed, target) = remove_outliers(
			interest_data_transformed, options.koef, target)
	interest_data_transformed = norm(interest_data_transformed)
	if options.target_column != None: 
		interest_data_transformed = concatenate((interest_data_transformed,
							 matrix(target).T), axis=1)

	# write data
	logger.info("Writing output...")
	if options.file_out:
		savetxt(options.file_out, interest_data_transformed, 
			delimiter=',')
	else:
		savetxt(filename.replace(".csv","_normalized.csv"), 
			interest_data_transformed, delimiter=',')
	logger.info("Complete!")





