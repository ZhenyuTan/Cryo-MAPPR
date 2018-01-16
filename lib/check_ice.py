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

		#Get dimensions: 
		x,y=np.shape(micfile.data)
	
		#If not square, make it a square with smallest dimension
		if x != y: 
			newdim=min(x,y)
			micfile=micfile.data[0:newdim,0:newdim] #np.reshape(micfile.data,(newdim,newdim))
		if x == y: 
			newdim=min(x,y)
			micfile=micfile.data
		
		#Calc. FT of input micrograph
		F1=fftpack.fft2(micfile)
		F2 = fftpack.fftshift( F1 )
		psd2D = np.abs( F2 )**2
		psd1D = radialProfile.azimuthalAverage(psd2D)

		#Get max value b/w 3.6 and 3.8 
		#To find x value entry for a given resolution: 
                # X= (Full XDIM)*(APix)/(Resolution shell wanted)
                # X= (3710 * 0.91)/(3.7) = 912
		Xtoget=round((newdim*apix)/(3.7))
		startX=Xtoget-5
		maxValue=0
		finishX=Xtoget+5
		while startX <= finishX:
			checkVal=psd1D[int(startX)]
			if checkVal>maxValue: 
				maxValue=checkVal
			startX=startX+1
	
		#Get adjacent value for non vitreous ice
		Xtoget=round((newdim*apix)/(4))
		nullValue=psd1D[int(Xtoget)]
		VitreousCheck=maxValue/nullValue
		if VitreousCheck > 0.995:
			badlist.append(mic) 
				
		#Can be used to plot with pyfits
		'''
		import pyfits
		micfile=mrcfile.open(mic)
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
		'''
		
	return badlist
