from PIL import Image
import mahotas as mh
import os
import sys 
import mrcfile
from scipy import fftpack
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import pylab as py

#====================
def checkmics(miclist,apix): 
	'''
	This will open micrographs to check for non-vitreous ice, where micrographs will be discarded if
	there are reflections of non-vitreous ice above a given threshold
	'''
	#Start bad list
	badlist=[]
	counter=0

	#Get micrograph directory
        if len(miclist)>0:
                micdirectory=miclist[0].split('/')
                del micdirectory[-1]
                micdirectory='/'.join(micdirectory)

        if len(badlist)>0:
                micdirectory=badlist[0].split('/')
                del micdirectory[-1]
                micdirectory='/'.join(micdirectory)

	#Open output file
	outfile=open('%s/vitreousIceCheck_output.txt' %(micdirectory),'w')

        #Write info line: 
        outfile.write('#Columns: micrograph name, ratio of 3.7A ring to 4A, Good/bad\n')

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
		psd1D = azimuthalAverage(psd2D)

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
			outfile.write('%s\t%f\t%s\n' %(mic,VitreousCheck,'Bad'))
		if VitreousCheck <= 0.995: 
			outfile.write('%s\t%f\t%s\n' %(mic,VitreousCheck,'Good'))	
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
	#Remove entries from good list
	goodlist=remove_entries_from_good_list(miclist,badlist)
	return goodlist,badlist

#===================
def remove_entries_from_good_list(ingoodlist,badlist): 
	outgoodlist=[]
	for entry in ingoodlist: 
		if entry in badlist: 
			continue
		outgoodlist.append(entry)
	return outgoodlist

#==================
def checkStats(inlist,badlist): 
	'''
	Get average & std dev for each mic, discard if number of pixels outside of 4 sigma is > ____%
	'''

	for mic in inlist: 
		micfile=mrcfile.open(mic)
		avg=np.average(micfile.data)
		stdDev=np.std(micfile.data)

		#print mic
		#print avg
		#print stdDev
		#print np.shape(micfile.data[(np.where(micfile.data>(avg+(3*stdDev))))])

	return inlist,badlist

#===================
def findIce(goodlist,badlist,apix,diameter,percentIceAllowed): 
	'''Loop over icefinder with a set of micrographs, returning good and bad lists
	(goodlist,badlist,params['apix'],params['diam'],percentIceAllowed)'''

	#Set up output lists: 
	miclistNoIce=[]
	miclistWithIce=[]
	percentAreaList=[]

	#Get micrograph directory
	if len(goodlist)>0: 
		micdirectory=goodlist[0].split('/')
		del micdirectory[-1]
                micdirectory='/'.join(micdirectory)

	if len(badlist)>0: 
		micdirectory=badlist[0].split('/')
		del micdirectory[-1]
		micdirectory='/'.join(micdirectory)

	outfile=open('%s/findIce_output.txt' %(micdirectory),'w')

	#Write info line: 
	outfile.write('#Columns: micrograph name, percentage area ice, Good/bad\n')

	#Hardcoded ice parameters (for now)
	particles_size_threshold=diameter*diameter #Suggested 128*128
	instensity_threshold=1.4
	gaussian_filter_sigma=25

	for mic in goodlist: 

		#Execute icefinder, returns percent covered in ice
		iceArea=icefinder(mic,particles_size_threshold, instensity_threshold,gaussian_filter_sigma)
		
		#Add to list	
		percentAreaList.append(iceArea)

		#Choose mic destination into good or bad list depending on % area covered
		if iceArea > percentIceAllowed: 
			miclistWithIce.append(mic)
			outfile.write('%s\t%f\t%s\n' %(mic,iceArea,'Bad'))
		if iceArea <= percentIceAllowed: 
			miclistNoIce.append(mic)
			outfile.write('%s\t%f\t%s\n' %(mic,iceArea,'Good'))
	
	#Output results
	bins = np.linspace(0, 10,100)
	plt.hist(percentAreaList,bins)
    	plt.title("Ratio_Histogram")
    	plt.xlabel("Value")
    	plt.ylabel("Frequency")
	plt.savefig('%s/percentIce_Histogram.png' %(micdirectory)) #plt.show()

	outfile.close()

	return miclistNoIce,miclistWithIce
			
#===================
def icefinder(inputMRC, particles_size_threshold, instensity_threshold,gaussian_filter_sigma):
    mrcOpen=mrcfile.open(inputMRC)
    im = mrcOpen.data #Replaced with mrcfile reading operation: io.ImageCollection(path) #read images from the path
    sz=np.size(im) # image size 
    inverted_im = np.power(2,16)-im # invert contrast 
    im_filtered = mh.gaussian_filter(inverted_im, gaussian_filter_sigma )# gaussian filter 
    im_ins_thresholded = (im_filtered> im_filtered.mean()+instensity_threshold*np.sqrt(im_filtered.var())) # binary the image based on intensity
    labeled, n_nucleus  = mh.label(im_ins_thresholded) #label the patches which has the value 1 
    sizes = mh.labeled.labeled_size(labeled) #calulate the size of each patch,simple count the piex numbers
    too_small = np.where(sizes < particles_size_threshold)
    labeled = mh.labeled.remove_regions(labeled, too_small) #remove the patches which is too small 
    labeled = mh.labeled.remove_bordering(labeled) # remove the patches which attach the edge 
    relabeled, n_left = mh.labeled.relabel(labeled) #relabel 
    sizes_cleared = mh.labeled.labeled_size(relabeled)
    ratio=(np.sum(sizes_cleared[1:])/sz)*100
    return ratio

#====================
def azimuthalAverage(image, center=None):
    """
    Calculate the azimuthally averaged radial profile.

    image - The 2D image
    center - The [x,y] pixel coordinates used as the center. The default is 
             None, which then uses the center of the image (including 
             fracitonal pixels).
    Source: http://www.astrobetter.com/wiki/python_radial_profiles 
    Repo: https://github.com/keflavich/image_tools/blob/master/image_tools/radialprofile.py
    """
    # Calculate the indices from the image
    y, x = np.indices(image.shape)

    if not center:
        center = np.array([(x.max()-x.min())/2.0, (x.max()-x.min())/2.0])

    r = np.hypot(x - center[0], y - center[1])

    # Get sorted radii
    ind = np.argsort(r.flat)
    r_sorted = r.flat[ind]
    i_sorted = image.flat[ind]

    # Get the integer part of the radii (bin size = 1)
    r_int = r_sorted.astype(int)

    # Find all pixels that fall within each radial bin.
    deltar = r_int[1:] - r_int[:-1]  # Assumes all radii represented
    rind = np.where(deltar)[0]       # location of changed radius
    nr = rind[1:] - rind[:-1]        # number of radius bin
    
    # Cumulative sum to figure out sums for each radius bin
    csim = np.cumsum(i_sorted, dtype=float)
    tbin = csim[rind[1:]] - csim[rind[:-1]]

    radial_prof = tbin / nr

    return radial_prof
