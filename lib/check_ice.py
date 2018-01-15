import numpy as np
import os
import sys 
import mrcfile
from scipy import fftpack
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import pylab as py
import radialProfile
import pyfits

#====================
def checkmics(miclist,apix): 

	'''
	This will open micrographs to check for non-vitreous ice, where micrographs will be discarded if
	there are reflections of non-vitreous ice above a given threshold
	'''
	#Start bad list
	badlist=[]

	#Loop over all micrographs
	for mic in miclist: 

		#Open MRC file: 
		micfile=mrcfile.open(mic)
		
		#Check that single image
		if micfile.is_single_image() is False: 
			print 'Error: Input micrograph %s has more than one file (is it a movie?). Skipping micrograph'
			badlist.append(micfile)
			continue

		#Get dimensions: np.shape(micfile.data)

		#Calc. FT of input micrograph
		F1=fftpack.fft2(micfile.data)
		F2 = fftpack.fftshift( F1 )
		psd2D = np.abs( F2 )**2
		print np.shape(psd2D)
		psd1D = radialProfile.azimuthalAverage(psd2D)
		print np.shape(psd1D)	
		#Can be used to plot with pyfits
		py.figure(1)
		py.clf()
		py.imshow( np.log10( micfile.data ), cmap=py.cm.Greys)
 		py.figure(2)
		py.clf()
		py.imshow( np.log10( psd2D ))
 		py.figure(3)
		py.clf()
		py.semilogy( psd1D )
		py.xlabel('Spatial Frequency')
		py.ylabel('Power Spectrum')
 		py.show()
		

		#Take circular average within resolution shell 3.6-3.8A, relative intensity to 3.1A

		
	return badlist
