import scipy.io
import numpy as np
import cv2
from random import randint
from math import sqrt

#Settings

#Number of Images per Class
num_train = 800000
num_test = 200

#Number of real output Images per Class
out_train = 20

#Size of the image
h = w = 16
t = 2

#output filename
base = 'tricircle_0.25-0.75_%ix%i'%(h,w)
end = '.arff'

train = base + '_%iimages.train'%num_train
test = base + '.test'


#write header for both files
for filename in [test, train]:
	with open(filename,'w') as f:
		f.write('''@RELATION Triangles and circles. 0 for cicles. 1 for triangles.\n\n''')

		for xi in range(w):
			for yi in range(h):
				f.write('''@ATTRIBUTE x%iy%i NUMERIC\n''' % (xi,yi))
		f.write('''@ATTRIBUTE class NUMERIC\n@DATA\n''')


cos = sqrt(3)/2.
sin = 1/2.

for i in range(num_test+num_train):
	#start with two empty images
	circle = np.zeros((h*t,w*t), np.uint8)
	triangle = np.zeros((h*t,w*t), np.uint8)


	#generate random circle
	radius = randint(3,t*min(h,w)/3)
	x = randint(1+radius, t*h-1-radius)
	y = randint(1+radius, t*h-1-radius)

	#draw circle
	cv2.circle(circle, (x,y), radius, 255, 2)


	#random triangle
	#3 random points?
	#cord = np.random.randint(1,t*min(h,w), (3,2), np.uint8)
	#convert to tuple
	#cordt = [(c[0],c[1]) for c in cord]

	cordt = [(x,y+radius), (x+cos*radius, y-sin*radius), (x-cos*radius, y-sin*radius)]
	cordt = [(int(c[0]), int(c[1])) for c in cordt]

	#draw triangle
	cv2.line(triangle, cordt[0], cordt[1], 255, 2)
	cv2.line(triangle, cordt[1], cordt[2], 255, 2)
	cv2.line(triangle, cordt[2], cordt[0], 255, 2)

	#scale image to w*h
	circle = cv2.resize(circle, (h,w), interpolation=cv2.INTER_AREA)
	triangle = cv2.resize(triangle, (h,w), interpolation=cv2.INTER_AREA)

	if i < num_test:
		cv2.imwrite('test_img/cir%i.png'%i, circle)
		cv2.imwrite('test_img/tri%i.png'%i, triangle)
	elif i < num_test + out_train:
		cv2.imwrite('train_img/cir%i.png'%i, circle)
		cv2.imwrite('train_img/tri%i.png'%i, triangle)

	#convert image to [0.1, 0.9]
	images = [circle, triangle]

	for j,img in enumerate(images):
		img = img/255.
		img *= 0.5
		img += 0.25
		fname = (train,test)[i<num_test]
		with open (fname, 'a') as f:
			for val in img.reshape(w*h):
				f.write('%f, '%val)
			f.write('%i\n'%j)


