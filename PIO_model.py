#!/usr/bin/env python

################################################
###  Michael Howland, mhowland@stanford.edu  ###
###  CTR Summer Program 2016                 ###
###  Predictive Inner Outer Model, Marusic   ###
###  2011 in WMLES                           ###
################################################

# Import libraries
import numpy as np
import scipy.io as sio
import struct

# Parameters
retau = 2000
plane_num = 1
d=3.5
folder = ''
folder = 'data/' + 'Retau' + str(retau) + '/' + 'plane_' + str(plane_num) + '/' + str(d) + 'd/'
#print folder

# Some operation functions
def rms(val):
	rms_val = np.sqrt(np.mean(np.square(val)))
	return rms_val
	

# Load calibration data
nt = 10
for t in range(1, nt+1):
	nz = 4608 #number of spanwise locations
	folder_bin = folder + 'cali_mat' + '_t' + str(t) + '.bin'
	#print folder_bin
	f = open(folder_bin, 'rb')
	data = np.fromfile(f, dtype=np.float64)
	#print len(data)
	nx = float(len(data))/(float(nz)+1)
	#print nx
	s = (nx, 1)
	tplus = np.zeros(s, dtype=np.float64)
	s = (nx, nz) 
	tau_w_star = np.zeros(s, dtype=np.float64)
	tau_rms = np.zeros(nz, dtype=np.float64)
	# Unpack the binary data
	for j in range(0, nz+1):
		if j == 0:
			tplus = data[0 : nx - 1]
		else:
			tau_w_star[:,j-1] = data[(j)*nx : (j+1)*nx]
			tau_rms[j-1] = rms(data[(j)*nx : (j+1)*nx])
			#tau_rms[j-1] = np.sqrt(np.mean(np.square(tau_w_star[:,j-1])))
	print "Time step plane = %s" % t
	print "Average tau_rms = %s" % np.mean(tau_rms)
	print "STD tau_rms = %s" % np.std(tau_rms)



