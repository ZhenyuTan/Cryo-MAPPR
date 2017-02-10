#!/usr/bin/env python
import math
import time
import linecache
import subprocess
import glob
import os
import shutil
import math
import sys 
import optparse
from optparse import SUPPRESS_HELP
import random
#=========================
def setupParserOptions():
        arglist = []
        for arg in sys.argv:
                arglist.append( arg )
	parser = optparse.OptionParser()
        progname = os.path.basename(arglist[0])
	usage=progname+"\n\nThis program will automatically find and pick particles in your dataset, outputing high resolution particles along with a report of your dataset.\n\nRequired Inputs:\n  -h, --help    show this help message and exit\n  --dir=FILE    Directory containing micrographs\n  --MW=INT      Approximate molecular weight of sample (kDa)\n  --kev=INT     Accelerating voltage in keV (Default=200)\n  --cs=FLOAT    Spherical aberration (mm) (Default=2.0)\n  --apix=FLOAT  Micrograph pixel size in Angstroms/pixel (Default=1.73)\n"
	parser.add_option("--dir",dest="dir",type="string",metavar="FILE",
                    help="Directory containing micrographs")
	parser.add_option("--MW",dest="MW",type="int",metavar="INT",
                    help="Approximate molecular weight of sample (kDa)")
	parser.add_option("--kev",dest="kev",type="int",metavar="INT",default=200,
                    help="Accelerating voltage in keV (Default=200)")
	parser.add_option("--cs",dest="cs",type="int",metavar="FLOAT",default=2.0,
                    help="Spherical aberration (mm) (Default=2.0)")
	parser.add_option("--apix",dest="apix",type="float",metavar="FLOAT",default=1.73,
                    help="Micrograph pixel size in Angstroms/pixel (Default=1.73)")
	parser.add_option("--out",dest="outdir",type="string",metavar="FILE",default='AutoPick',
                    help="Name of output directory (absolute path) for particle picking outputs to be located (Default= Located in input directory as 'AutoPick'")
	parser.add_option("--negstain",dest="negstain",action="store_true",default=False,
                    help="Flag to indicate dataset is from negative stain")
	parser.add_option("--phasePlate",dest="phasePlate",action="store_true",default=False,
                    help="Flag to indicate dataset is from Volta Phase Plate")
	parser.add_option("--continue", action="store_true",dest="continue",default=False,
                    help="Continue previous run")
	parser.add_option("--partLimit",dest="partLimit",type="int",metavar="INT",default=0,
                    help="Particle limit required in order to start initial 2D classification (Default: Negative Stain=5000, Cryo-EM=20000)")
	parser.add_option("--invertOverride",dest="invert",action="store_true",default=False,
                    help="Override automatic contrast inversion of images to create white particles on black background ")
	parser.add_option("--startingDiameter",dest="dogStartDiam",type="int",metavar="INT",default=0,
                    help="Initial starting diameter for DoG Picker search (Angstroms) (Default=MW/2)")
	parser.add_option("--finalDiameter",dest="dogFinDiam",type="int",metavar="INT",default=0,
                    help="Final diameter for DoG Picker search (Angstroms) (Default=MW)")
	parser.add_option("--dogThresh",dest="dogThresh",type="float",metavar="FLOAT",default=0.3,
                    help="DoG Picker Threshold (Default=0.3)")
	parser.add_option("--binning",dest="binning",type="int",metavar="INT",default=0,
                    help="Binning factor for micrographs to be used for all processing (Default=Binning factor that yields pixel size between 5 - 8 Angstroms)")
	parser.add_option("--resmin",dest="resmin",type="int",metavar="INT",default=30,
                    help="Resolution cutoff for class averages to be considered 'good' (Default=30 Angstroms")
	parser.add_option("--microresmin",dest="microresmin",type="int",metavar="INT",default=10,
                    help="Resolution cutoff for micrographs to be considered 'good' (Default=10 Angstroms")
	parser.add_option("-d", action="store_true",dest="debug",default=False,
                    help="debug")
        options,args = parser.parse_args()
        if len(args) > 0:
                parser.error("Unknown commandline options: " +str(args))
        if len(sys.argv) <=3:
		print "usage: " + usage
                print "Please run '" + progname + " -h' for full option listing\n"
		sys.exit()
        params={}
        for i in parser.option_list:
                if isinstance(i.dest,str):
                        params[i.dest] = getattr(options,i.dest)
        return params

#=============================
def checkConflicts(params):

	try:
    		import numpy
	except ImportError:
    		raise ImportError,"numpy is required to run this program."
		sys.exit()
	try:
                import scipy
        except ImportError:
                raise ImportError,"scipy is required to run this program."
                sys.exit()
	try:
                import PIL
        except ImportError:
                raise ImportError,"Python Image Library (PIL, pillow, etc.)  is required to run this program."
                sys.exit()

	gctf=subprocess.Popen("which Gctf", shell=True, stdout=subprocess.PIPE).stdout.read().strip()
	if not gctf: 
		print '\nError: Gctf command not found. Please install and include in PATH as "Gctf"\n'
		sys.exit()
	gauto=subprocess.Popen("which Gautomatch", shell=True, stdout=subprocess.PIPE).stdout.read().strip()
        if not gauto:
                print '\nError: Gautomatch command not found. Please install and include in PATH as "Gautomatch"\n'
		sys.exit()
	spider=subprocess.Popen("which spider", shell=True, stdout=subprocess.PIPE).stdout.read().strip()
        if not spider:
                print '\nError: SPIDER not found. Please install and include in PATH as "spider"\n'
		sys.exit()
	eman2=subprocess.Popen("which e2proc2d.py", shell=True, stdout=subprocess.PIPE).stdout.read().strip()
        if not eman2:
                print '\nError: EMAN2 not found. Please install and include in PATH"\n'
		sys.exit()
	dogpick=subprocess.Popen("which ApDogPicker.py", shell=True, stdout=subprocess.PIPE).stdout.read().strip()
        if not eman2:
                print '\nError: DoG Picker not found. Please install from http://emg.nysbc.org/redmine/projects/software/wiki/DoGpicker and include in PATH"\n'
                sys.exit()
	

#===================
def testPicking(newworkingdir,initial_micros,startingMicNum,negativeStain,apix,mpi,j_threads):

        if os.path.exists(newworkingdir): 
		shutil.rmtree(newworkingdir)
	os.makedirs(newworkingdir)
        stacks,micros,boxsize,workingfolders,diameter=Gauss_pick_micrographs(newworkingdir,initial_micros,startingMicNum,apix,binning,multiplier,starting_diameter,final_diameter,xdim,overlap_multiplier,boxsizemultiplier,startingthresh,finalthresh,in_mics,mpi,j_threads)

	return float(subprocess.Popen('cat %s/*/Import/*.box | wc -l' %(newworkingdir),shell=True, stdout=subprocess.PIPE).stdout.read().strip()),micros

#====================
def pick_run_2D(newworkingdir,initial_micros,startingMicNum,negativeStain,apix,in_mics,mpi,j_threads,gpus,reslim):

	#Create new directory
	os.makedirs(newworkingdir)
	os.makedirs('%s/gctf' %(newworkingdir))

	#Pick from micrographs using gaussian picking over a range of diameters 
	output_stacks,micros,finalbox,dogdirs,diameter=Gauss_pick_micrographs(newworkingdir,initial_micros,startingMicNum,apix,binning,multiplier,starting_diameter,final_diameter,xdim,overlap_multiplier,boxsizemultiplier,startingthresh,finalthresh,in_mics,mpi,j_threads)

	#Symbolically link micrographs to the directory
	outdir=newworkingdir
	workingMicList=[]
        counter=0
        maxnum=len(micros)
        cs=params['cs']
        if maxnum > len(micros):
                maxnum=len(micros)
        while counter< maxnum:
                gpucounter=1
                while gpucounter <= gpus:
                        if counter+gpucounter < maxnum: 
				if not os.path.exists('%s/gctf/process%i' %(outdir,gpucounter)):
                                	os.makedirs('%s/gctf/process%i' %(outdir,gpucounter))
				os.symlink(micros[counter+gpucounter],'%s/gctf/process%i/%s' %(outdir,gpucounter,micros[counter+gpucounter].split('/')[-1]))
                        	workingMicList.append('%s/gctf/process%i/%s' %(outdir,gpucounter,micros[counter+gpucounter].split('/')[-1]))
                        gpucounter=gpucounter+1
                counter=counter+gpus
        gpucounter=1
        wait=False
        while gpucounter <= gpus:
                microstarfile=estimate_CTF_gctf('%s/gctf/process%i' %(outdir,gpucounter),apix,kev,cs,phase,wait,gpucounter-1)
                gpucounter=gpucounter+1

        #Wait for Gctf to finish the combine star files into folder [outdir]/gctf/
        for mic in workingMicList:
                isdone=0
                while isdone == 0:
                        if os.path.exists('%s.ctf' %(mic[:-4])):
                                isdone=1
	starlist=[]
        gpucounter=1
        while gpucounter <=  gpus:
		isdone=0
		while isdone == 0: 
			if os.path.exists('%s/gctf/process%i/micrographs_all_gctf.star' %(outdir,gpucounter)): 
				starlist.append('%s/gctf/process%i/micrographs_all_gctf.star' %(outdir,gpucounter))
                		isdone=1
		gpucounter=gpucounter+1
        time.sleep(5)
	combineStarFiles(starlist,'%s/gctf/micrographs_all_gctf.star' %(outdir))
	'''
        counter=0
        maxnum=len(micros)
        while counter< maxnum:
                os.symlink(micros[counter],'%s/gctf/%s' %(newworkingdir,micros[counter].split('/')[-1]))
                counter=counter+1

	#Run ctf estimation on selected micrographs
        microstarfile=estimate_CTF_gctf('%s/gctf/' %(newworkingdir),apix,kev,cs,phase,False,0)
	#Merge picked particles into a single STAR file: merged_star.star
	mergedStarFile='%s/merged_star.star' %(newworkingdir)
	mergedwrite=open(mergedStarFile,'a')
	mergedwrite.write('\ndata_\n\nloop_\n_rlnCoordinateX #1\n_rlnCoordinateY #2\n_rlnImageName #3\n_rlnMicrographName #4\n_rlnVoltage #5\n_rlnDefocusU #6\n_rlnDefocusV #7\n_rlnDefocusAngle #8\n_rlnSphericalAberration #9\n_rlnCtfBfactor #10\n_rlnCtfScalefactor #11\n_rlnPhaseShift #12\n_rlnAmplitudeContrast #13\n_rlnMagnification #14\n_rlnDetectorPixelSize #15\n_rlnCtfFigureOfMerit #16\n')
	mergedwrite.close()

	#For each diameter, convert CTF and extract particles. Then place into the merged star file
	for stack in output_stacks:
		newstar=transferStar('%s/gctf/micrographs_all_gctf.star' %(newworkingdir),'%s/%s/micrographs_all_gctf.star' %(newworkingdir,stack),apix,apix*binning,'%s/%s' %(newworkingdir,stack),reslim)
                diameter=extractRelion(newstar,apix*binning,finalbox,mpi*j_threads,'%s/%s/Import/' %(newworkingdir,stack),'.box','%s/%s/' %(newworkingdir,stack))
        	newStarFile('%s/%s/Import/Extract/particles.star' %(newworkingdir,stack),mergedStarFile,'')

	#Get number of particles to calculate number of classes for 2D classification
	numParticles=float(subprocess.Popen('cat %s | grep Extract | wc -l' %(mergedStarFile),shell=True, stdout=subprocess.PIPE).stdout.read().strip())
	numClasses=round(numParticles/200)
	if numClasses <5:
        	numClasses=5
	if numClasses > 200:
        	numClasses=200

	#Specify only flip phases for negative stain
	if negativeStain == 0: 
		ctfcorrection='--ctf'
	if negativeStain == 1: #YES
		ctfcorrection='--ctf --only_flip_phases'
	
	#Run Relion 2D classification	
	runRelion2D('%s/' %(newworkingdir),mergedStarFile.split('/')[-1],apix*binning,mpi,numClasses,round(finalbox*.45)*apix*binning,'--dont_check_norm',j_threads,gpu,ctfcorrection,True)
	
	#Format output files: Remove pixel values from any empty classes
	cmd='mv %s/Class2D/run_it025_classes.mrcs %s/Class2D/run_it025_classes_orig.mrcs' %(newworkingdir,newworkingdir)
	subprocess.Popen(cmd,shell=True).wait()
	cmd='relion_image_handler --i %s/Class2D/run_it025_classes_orig.mrcs --o noNan --remove_nan --replace_nan 0' %(newworkingdir)
	subprocess.Popen(cmd,shell=True).wait()
	cmd='mv %s/Class2D/run_it025_classes_orig_noNan.mrcs %s/Class2D/run_it025_classes.mrcs' %(newworkingdir,newworkingdir)
	subprocess.Popen(cmd,shell=True).wait()
	#Return output stack of averages
	return '%s/Class2D/run_it025_classes.mrcs' %(newworkingdir),diameter,finalbox
	'''
	return dogdirs,diameter,finalbox

#====================
def newStarFile(instar,appendstar,newdir): 

	ino=open(instar,'r')
	outo=open(appendstar,'a')
	for line in ino: 
		if len(line)<40: 
			continue
		l=line.split()
		particle=l[2]
		newparticle=particle.split('@')[0]+'@'+newdir+particle.split('@')[-1]
		l[2]=newparticle
		newline='\t'.join(l)+'\n'
		outo.write(newline)					
	ino.close()
	outo.close()

#====================
def find_diam(info):

        peak=0
        coord=0
	info_open=open(info,'r')

        for line in info_open:
                if len(line) == 0:
                        continue
                if line[1] == ';':
                        continue
                intensity=float(line.split()[2])
                if intensity > peak:
                        peak=intensity
                        coord=float(line.split()[3])
	info_open.close()
	info_open=open(info,'r')
	zero=-5
        for line in info_open:
                if len(line) == 0:
                        continue
                if line[1] == ';':
                        continue
		if float(line.split()[3]) > coord: 
			if float(line.split()[2]) <=0.1:
				if zero == -5: 
					zero=float(line.split()[3])
					crossing=float(line.split()[3])
	info_open.close()
	info_open=open(info,'r')
	badflag=0
	for line in info_open:
                if len(line) == 0:
                        continue
                if line[1] == ';':
                        continue
		if float(line.split()[3]) > crossing: 
			if float(line.split()[2]) > 0.05: 
				badflag=1
        return peak,coord,zero*2,badflag

#====================
def get_templates(stack,numimages,mindiam,maxdiam):

	#Convert stack to spider: 	
	cmd='proc2d %s %s norm=0,1 inplace' %(stack,stack)
        subprocess.Popen(cmd,shell=True).wait()
	cmd='proc2d %s %s.spi spiderswap' %(stack,stack[:-4])
	subprocess.Popen(cmd,shell=True).wait()
	cmd='do lb1 [avg]=1,%i\n' %(numimages)
	cmd+='RO I\n'
	cmd+='%s@{*********[avg]}\n' %(stack[:-4])
	cmd+='%s_roi@{**********[avg]}\n' %(stack[:-4])
	cmd+='lb1\n'
	runSpider(cmd)
	o1=open('%s_outInfo.txt' %(stack[:-4]),'w')
	sellist='%s_selList.txt' %(stack[:-4])
	o2=open(sellist,'w')

	counter=1
	while counter <=numimages: 
		cmd='RO SD\n'
		cmd+='%s@%i\n' %(stack[:-4],counter)
		cmd+='%s_rosd\n' %(stack[:-4])
		cmd+='%s_rosd_info\n' %(stack[:-4])
		runSpider(cmd)
		peak,coord,diam,badflag=find_diam('%s_rosd_info.spi' %(stack[:-4]))
		keep=0
		if diam < maxdiam: 
			if diam > mindiam: 
				if peak > 0.5: 
					if badflag == 0: 
						o2.write('%i\n' %(counter-1))
						keep=1
		o1.write('%i\t%f\t%f\t%f\t%i\t%i\n' %(counter,peak,coord,diam,badflag,keep))
		os.remove('%s_rosd_info.spi' %(stack[:-4]))
		os.remove('%s_rosd.spi' %(stack[:-4]))
		counter=counter+1

	cmd='e2proc2d.py %s.img %s_sel.img --select %s' %(stack[:-4],stack[:-4],sellist)
	print cmd 
	subprocess.Popen(cmd,shell=True).wait()

	return '%s_sel.img' %(stack[:-4])

#=====================================
def runRelion2D(directory,inputfile,apix,mpi,numClasses,diamIn,extraOption,j_threads,gpu,ctf,wait): 

	os.makedirs('%s/Class2D/' %(directory))
	#Write relion extraction command: 
	if diamIn == 1:
		diam=float(directory.split('diam')[-1][:-1])
		diam=(diam*.2)+diam
        if diamIn > 1: 
		diam=diamIn
	shell=subprocess.Popen("echo $SHELL", shell=True, stdout=subprocess.PIPE).stdout.read()
        shell=shell.split('/')[-1][:-1]
	gpuin=''
	if gpu > 0: 
		gpuin=' --gpu'	
        cmd='#!/bin/%s -x\n' %(shell)
	#cmd+='cd %s/\n' %(directory)
	cmd+='mpirun -np %i relion_refine_mpi --o %s/Class2D/run --i %s/%s --dont_combine_weights_via_disc --pool 3 --iter 25 --tau2_fudge 2 --particle_diameter %i --K %i --flatten_solvent  --zero_mask  --oversampling 1 --psi_step 10 --offset_range 5 --offset_step 2 --norm --scale  --j %i --angpix %f %s %s %s' %(mpi,directory,directory,inputfile,diam,numClasses,j_threads,apix,extraOption,gpuin,ctf)

	print 'mpirun -np %i relion_refine_mpi --o %s/Class2D/run --i %s/%s --dont_combine_weights_via_disc --pool 3 --iter 25 --tau2_fudge 2 --particle_diameter %i --K %i --flatten_solvent  --zero_mask  --oversampling 1 --psi_step 10 --offset_range 5 --offset_step 2 --norm --scale  --j %i --angpix %f %s %s %s' %(mpi,directory,directory,inputfile,diam,numClasses,j_threads,apix,extraOption,gpuin,ctf)

	removefile('%s/relionrun.com' %(directory))
        removefile('%s/relionLog.log' %(directory))

        relionFile = open('%s/relionrun.com' %(directory),'w')
        relionFile.write(cmd)
        relionFile.close()

        cmd = 'chmod +x %s/relionrun.com' %(directory)
        subprocess.Popen(cmd,shell=True).wait()

	if wait is True:
        	cmd = '%s/relionrun.com > %s/relionLog.log' %(directory,directory)
		subprocess.Popen(cmd,shell=True).wait()
        	removefile('%s/relionrun.com'%(directory))
        	removefile('%s/relionLog.log' %(directory))
	if wait is False: 
		cmd = '%s/relionrun.com > %s/relionLog.log' %(directory,directory)
                subprocess.Popen(cmd,shell=True)

#============================
def extractRelion(starfile,newapix,finalbox,mpi,boxdir,boxext,chdir): 
	directory=boxdir
	boxsize=finalbox
	if boxsize %2 == 1: 
		boxsize=boxsize+1
	os.makedirs('%s/Extract/' %(directory))
	#Write relion extraction command: 
	shell=subprocess.Popen("echo $SHELL", shell=True, stdout=subprocess.PIPE).stdout.read()
        shell=shell.split('/')[-1][:-1]
	include_mpi=''
	relionmpi=''
	if mpi > 1: 
		include_mpi='mpirun -np %i ' %(mpi)
		relionmpi='_mpi '
        cmd='#!/bin/%s -x\n' %(shell)
	cmd+='cd %s\n' %(chdir)
	cmd+='%srelion_preprocess%s --i %s --coord_dir %s --coord_suffix %s --part_star %s/Extract/particles.star --part_dir %s/Extract/ --extract --extract_size %i  --norm --bg_radius %i --white_dust 5 --black_dust 5' %(include_mpi,relionmpi,starfile,boxdir,boxext,boxdir,boxdir,int(boxsize),round((0.6*boxsize)/2))
	removefile('relionrun.com')
        removefile('relionLog.log')

        relionFile = open('relionrun.com','w')
        relionFile.write(cmd)
        relionFile.close()
        
	cmd = 'chmod +x relionrun.com'
        subprocess.Popen(cmd,shell=True).wait()

        cmd = './relionrun.com > relionLog.log'
        subprocess.Popen(cmd,shell=True).wait()
	
	removefile('relionrun.com')
        removefile('relionLog.log')

	return round(0.6*boxsize)
#============================
def Gauss_pick_micrographs(workingdir,micros,apix,binning,multiplier,starting_diameter,final_diameter,xdim,overlap_multiplier,boxsizemultiplier,startingthresh,finalthresh,in_mics,mpi,j_threads): 

	parallel_relion_image_handler(micros,workingdir,' --angpix %f --rescale_angpix %f' %(apix,apix*binning),(mpi*j_threads))

	newpix=apix*binning

	current_diameter=starting_diameter
	stacks=[]
	pickedmicros=[]
	outpickdirs=[]

	maxnum=((((xdim*apix)/current_diameter)*((xdim*apix)/current_diameter))/overlap_multiplier)
	boxsize=round((current_diameter*boxsizemultiplier)/newpix)
	if boxsize % 2 == 1: 
		boxsize= boxsize+1
	miccounter=0
	current_thresh=startingthresh
        outdir='%s' %(workingdir)
	outpickdirs.append(outdir)
        os.makedirs('%s/Import' %(outdir))
	while miccounter < len(micros): 
		threadnum=0
                workingpicks=[]
                while threadnum < mpi*j_threads: 
                	if miccounter+threadnum<len(micros): 
                       		mic=micros[miccounter+threadnum].split('/')[-1]
                        	diam=current_diameter/newpix
                        	cmd='ApDogPicker.py --image=%s/%s  --outfile=%s/%s_picks.txt --max-peaks=2000 --thresh=%f --diam=%f > tmpout.log' %(workingdir,mic,workingdir,mic[:-4],current_thresh,diam)
				subprocess.Popen(cmd,shell=True)
                        	workingpicks.append('%s/%s_picks.txt' %(workingdir,mic[:-4]))
                       		pickedmicros.append('%s/%s_picks.txt' %(workingdir,mic[:-4]))
			threadnum=threadnum+1
                for checkmic in workingpicks:
		 	isdone=0
			while isdone == 0:
				if os.path.exists(checkmic):
                                	isdone=1
		for mic in workingpicks: 		
			if os.path.exists(mic):
				#cmd='ln -s %s.mrc %s/%s.mrc' %(mic[:-10],outdir,mic.split('/')[-1][:-10])
				#subprocess.Popen(cmd,shell=True).wait()
				boxfile=picktxt_to_box(mic,boxsize,math.floor(xdim/binning),math.floor(ydim/binning),'%s/Import/%s.box' %(outdir,mic.split('/')[-1][:-10]))
						#if 'stack_thresh%.2f_diam%.0f' %(current_thresh, current_diameter) not in stacks: 
						#	stacks.append('stack_thresh%.2f_diam%.0f' %(current_thresh, current_diameter))
				removefile(mic)
				os.remove('%s-finalmap.jpg' %(mic[:-10]))
				removefile('%s-finalmap.jpg' %(mic[:-10]))
				removefile('%s-finalpicks.jpg' %(mic[:-10]))
				removefile('%s-map01.jpg' %(mic[:-10]))
				removefile('%s-map02.jpg' %(mic[:-10]))
				removefile('%s-picks01.jpg' %(mic[:-10]))
        		       	removefile('%s-picks02.jpg' %(mic[:-10]))
				removefile('tmpout.log')
		miccounter=miccounter+(mpi*j_threads)
	return stacks,micros,boxsize,outpickdirs,current_diameter-diameter_increment

#===============================
def convertToRelionCTF(workingdir,ctffile,phase,ampcontrast,cs,kev,reslim,mag,detectorpix):

	relionOut = writeRelionHeader(phase)

	out = open('%s/all_micrographs_ctf.star' %(workingdir),'w')

	ctf = open(ctffile,'r')

	if phase is False:
		FlagAdd=0
	if phase is True:
		FlagAdd=1
	
	for line in ctf:
		l = line.split()
		if float(l[-2]) > reslim:
			continue
		micro=l[0]	
		#Prepare micrograph name
		microname = '%s' %(micro.split()[0])

		#Get defocus information
		df1 = float(l[1])
		df2 = float(l[2])
		astig = float(l[3])
		crosscorr = float(l[4])
		if phase is True:
			phaseShift=float(l[6])
		if phase is False:
			phaseShift=float(l[6])
			relionOut+='%s  %.6f  %.6f  %.6f  %.6f  %.6f  %.6f  %.6g  %.6f  %.6f\n' %(microname,df1,df2,astig,kev,cs,ampcontrast,mag,detectorpix,crosscorr)
		if phase is True:
			relionOut+='%s  %.6f  %.6f  %.6f  %.6f  %.6f  %.6f  %.6g  %.6f  %.6f  %.6f\n' %(microname,df1,df2,astig,kev,cs,ampcontrast,mag,detectorpix,crosscorr,phaseShift)
		ctflog = '%s/%s_ctffind3.log' %(workingdir,microname[:-4])

                #Open new ctf log file
                ctf='\n'
                ctf+=' CTF DETERMINATION, V3.5 (9-Mar-2013)\n'
                ctf+=' Distributed under the GNU General Public License (GPL)\n'
                ctf+='\n'
                ctf+=' Parallel processing: NCPUS =         4\n'
                ctf+='\n'
                ctf+=' Input image file name\n'
                ctf+='%s\n' %(microname)
                ctf+='\n'
                ctf+='\n'
                ctf+=' Output diagnostic file name\n'
                ctf+='%s.ctf\n'%(microname[:-4])
                ctf+='\n'
                ctf+='\n'
                ctf+=' CS[mm], HT[kV], AmpCnst, XMAG, DStep[um]\n'
                ctf+='  %.1f    %.1f    %.2f   %.1f    %.3f\n' %(cs,kev,ampcontrast,mag,detectorpix)
                ctf+='\n'
                ctf+='\n'
                ctf+='      DFMID1      DFMID2      ANGAST          CC\n'
                ctf+='\n'
                ctf+='    %.2f\t%.2f\t%.2f\t%.5f\t%.4f\tFinal Values\n' %(df1,df2,astig,crosscorr,phaseShift)

                outctf = open(ctflog,'w')
                outctf.write(ctf)
                outctf.close()

 
	out.write(relionOut)

#================================
def writeRelionHeader(phase):

	relion='\n'
	relion+='data_\n'
	relion+='\n'
	relion+='loop_\n'
	relion+='_rlnMicrographName #1\n'
	relion+='_rlnDefocusU #2\n'
	relion+='_rlnDefocusV #3\n'
	relion+='_rlnDefocusAngle #4\n'
	relion+='_rlnVoltage #5\n'
	relion+='_rlnSphericalAberration #6\n' 
	relion+='_rlnAmplitudeContrast #7\n'
	relion+='_rlnMagnification #8\n'
	relion+='_rlnDetectorPixelSize #9\n'
	relion+='_rlnCtfFigureOfMerit #10\n' 
	if phase is True:
		relion+='_rlnPhaseShift #11\n'

	return relion

#==========
def removefile(infile): 
	if os.path.exists(infile): 
		os.remove(infile)

#==========
def estimate_CTF_gctf(ctfdir,apix,kev,cs,phase,wait,gpuid):

	removefile('micrographs_all_gctf.star') 
	ctfphase=''
	if phase is True: 
		ctfphase='--phase_shift_L 0 --phase_shift_H 180 --phase_shift_T 2'

        cmd = 'Gctf --apix %f --kV %i --Cs %f %s --ctfstar %s/micrographs_all_gctf.star --gid %i --logsuffix _ctffind3.log --do_unfinished 0 %s/*.mrc > %s/gctf.log' %(apix,kev,cs,ctfphase,ctfdir,gpuid, ctfdir,ctfdir)
	if wait is True: 
		subprocess.Popen(cmd,shell=True).wait()
		return '%s/micrographs_all_gctf.star' %(ctfdir)
	if wait is False: 
		subprocess.Popen(cmd,shell=True)
		return ''

#============================
def compare_avgs_to_avgs(avglist,relionDir,apix,resthresh,resthresh_incr,MinMatch):
	refavg=0
	if os.path.exists(relionDir): 
		shutil.rmtree(relionDir)
	os.makedirs(relionDir)
	print avglist
	goodAvgs=open('%s/template_candidates.txt' %(relionDir),'w')
	while refavg < len(avglist):
	        inavg=0
		while inavg < len(avglist): 
			if inavg == refavg: 
				inavg=inavg+1
				continue
        		TemplateAvgs=avglist[refavg]
			QueryAvgs=avglist[inavg]
			num1,num2,newapsh,boxsize=align_avgs_to_avgsSPI(QueryAvgs,TemplateAvgs,'%s' %(relionDir),'match_iter%iRef_vs_iter%iQ.spi' %(refavg,inavg),0,0,0)	
			fsclist=getFSC('%s/match_iter%iRef_vs_iter%iQ' %(relionDir,refavg,inavg),num1,newapsh,apix)
			recordGoodAvgs(newapsh,fsclist,resthresh,thresh_incr,goodAvgs,refavg,inavg)
			inavg=inavg+1
		refavg=refavg+1
	goodAvgs.close()
	uniqList=getUniqueGoodAvgs('%s/template_candidates.txt' %(relionDir),MinMatch)
	groupList=groupMatchedAvgs(uniqList,'%s/template_candidates.txt' %(relionDir),MinMatch)
	newoutdir='%s/' %(relionDir.split('DoG')[0]) 
	showAvgs(groupList,'%s/template_stack' %(relionDir),newoutdir,boxsize,avglist)

#===============================
def showAvgs(groupList,outfile,aligndir,box,avglist):
	counter=1
	partcounter=0
	for item in groupList: 
		if item == '//':
			align_avgs_to_avgsSPI('%s_%i_unaligned.mrcs' %(outfile,counter),'%s_%i_unaligned.mrcs' %(outfile,counter),'','%s_%i_aligned.spi' %(outfile,counter),partcounter,1,box)
			#Create stack with only the aligned averages
			avgcounter=1
			while avgcounter <= partcounter*2: 
				if avgcounter == 1: 
					cmd='e2proc2d.py %s_%i_aligned.spi %s.mrcs --first=0 --last=0' %(outfile,counter,outfile)
					subprocess.Popen(cmd,shell=True).wait()
				cmd='e2proc2d.py %s_%i_aligned.spi %s_%i.mrcs --first=%i --last=%i' %(outfile,counter,outfile,counter,avgcounter-1,avgcounter-1)
				subprocess.Popen(cmd,shell=True).wait()
				avgcounter=avgcounter+2
			counter=counter+1
			partcounter=0
			continue
		iteration=int(item.split('-')[0])+1
		avgnum=int(item.split('-')[1])-1
		cmd='e2proc2d.py %s %s_%i_unaligned.mrcs --first=%i --last=%i' %(avglist[iteration-1],outfile,counter,avgnum,avgnum)
		subprocess.Popen(cmd,shell=True).wait()	
		partcounter=partcounter+1

#===============================
def groupMatchedAvgs(uniqList,inlist,MinMatch):
	grouped=[]
	for uniq in uniqList: 
		for line in open(inlist,'r'): 
			if line.split()[0] == uniq:
				if uniq not in grouped: 
					grouped.append(uniq)
				if line.split()[1] not in grouped: 
					grouped.append(line.split()[1])
			if line.split()[1] == uniq:
                                if uniq not in grouped:
                                        grouped.append(uniq)
                                if line.split()[0] not in grouped:
                                        grouped.append(line.split()[0])
		if grouped[-1] != '//': 
			grouped.append('//')
	return grouped
		
#================================
def searchArray(array,qq): 
	x=0
	y=0
	xP='None'
	yP='None'
	xdim,ydim=array.shape
	while x < xdim: 
		while y < ydim: 
			if array[x,y] == qq: 
				xP=x
				yP=y
			y=y+1
		x=x+1
	return xP,yP

#===============================
def getUniqueGoodAvgs(inlist,MinMatch):
	in1=open(inlist,'r')
	in1list=[]
	uniqList=[]
	for line in in1: 
		entry=line.split()[0]
		entry2=line.split()[1]
		in1list.append(entry)
		in1list.append(entry)	
		if entry not in uniqList: 
			uniqList.append(entry)
		if entry2 not in uniqList: 
			uniqList.append(entry2)
	in1.close()
	return uniqList

#===============================
def recordGoodAvgs(apsh,fsc,resthresh,resthresh_incr,outfile,templatenum,querynum):

	fopen=open(fsc,'r')	
	counter=1
	for line in fopen: 
		if float(line.split()[0]) < resthresh:
#			if float(line.split()[0]) > resthresh-resthresh_incr: 
			apshline=linecache.getline(apsh,counter+2)
			refnum=int(float(apshline.split()[5]))
			outfile.write('%i-%i\t%i-%i\n' %(templatenum,refnum,querynum,counter))
		counter=counter+1		
	fopen.close()

#===============================
def getFSC(matchedAvgs,numrefs,newapsh,apix):

	counter=1
	o1=open('%s_frc.spi' %(matchedAvgs),'w')	
	while counter <= numrefs*2: 

		cmd='RF\n'
		cmd+='%s@%i\n' %(matchedAvgs,counter)
		cmd+='%s@%i\n' %(matchedAvgs,counter+1)
		cmd+='1\n'
		cmd+='0.8,1.2\n'
		cmd+='%s_fsc\n' %(matchedAvgs)
		runSpider(cmd)

		resolution=getRes('%s_fsc.spi' %(matchedAvgs))
		resolution=apix/resolution
		o1.write('%f\n' %(resolution))
		removefile('%s_fsc.spi' %(matchedAvgs))
		counter=counter+2
	o1.close()
	return '%s_frc.spi' %(matchedAvgs)
#===================
def getRes(infile): 

	in1=open(infile,'r')
	resolution=0
	flag=0
	for line in in1: 
		if line[1] == ';': 
			continue
		freq=float(line.split()[2])
		frc=float(line.split()[4])
		if frc < 0.5: 
			if flag == 0:	
				resolution=freq
				flag=1
	in1.close()
	return resolution 

#====================
def runSpider(lines):
       spifile = "currentSpiderScript.spi"
       if os.path.isfile(spifile):
               os.remove(spifile)
       spi=open(spifile,'w')
       spi.write("MD\n")
       spi.write("TR OFF\n")
       spi.write("MD\n")
       spi.write("VB OFF\n")
       spi.write("MD\n")
       spi.write("SET MP\n")
       spi.write("(4)\n")
       spi.write("\n")
       spi.write(lines)

       spi.write("\nEN D\n")
       spi.close()
       spicmd = "spider spi @currentSpiderScript"
       spiout = subprocess.Popen(spicmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).stderr.read()
       output = spiout.strip().split()
       if "ERROR" in output:
               print "Spider Error, check 'currentSpiderScript.spi'\n"
               sys.exit()
       # clean up
       os.remove(spifile)
       if os.path.isfile("LOG.spi"):
               os.remove("LOG.spi")
       resultf = glob.glob("results.spi.*")
       if resultf:
               for f in resultf:
                       os.remove(f)
#===============================
def align_avgs_to_avgsSPI(stack1,stack2,workingfolder,outfile,num1,num2,box):

	if not os.path.exists('%s.spi' %(stack1[:-5])): 
		cmd='e2proc2d.py %s %s.spi --outtype=spi --writejunk' %(stack1,stack1[:-5])
		subprocess.Popen(cmd,shell=True).wait()
	if not os.path.exists('%s.spi' %(stack2[:-5])):
		cmd='e2proc2d.py %s %s.spi --outtype=spi --writejunk'%(stack2,stack2[:-5])
		subprocess.Popen(cmd,shell=True).wait()	

	if num1 == 0: 
		num1=int(subprocess.Popen('cat %s_model.star | grep _rlnNrClasses' %(stack1[:-13]), shell=True, stdout=subprocess.PIPE).stdout.read().split()[-1])	
	if num2 == 0: 
		num2=int(subprocess.Popen('cat %s_model.star | grep _rlnNrClasses' %(stack2[:-13]), shell=True, stdout=subprocess.PIPE).stdout.read().split()[-1])
	if box == 0: 
		box=int(subprocess.Popen('cat %s_model.star | grep _rlnOriginalImageSize' %(stack2[:-13]), shell=True, stdout=subprocess.PIPE).stdout.read().split()[-1])

	removefile('inputfilespitmp.txt')
	if os.path.exists('tmpdirspi'):
		shutil.rmtree('tmpdirspi')		
	#Write spider inputs into temp file
	o1=open('inputfilespitmp.txt','w')
	print 'INput file inputfilespitmp.txt: '
	o1.write('%s\n' %(stack1[:-5]))
	print '%s' %(stack1[:-5])
	o1.write('%i\n' %(num1))
	print '%i' %(num1)
	o1.write('%s\n' %(stack2[:-5]))
	print '%s' %(stack2[:-5])
	o1.write('%i\n' %(num2))
	print '%i' %(num2)
	o1.write('%i\n' %(box))
	print '%i' %(box)
	o1.write('tmpdirspi\n')
	o1.close()

	
	cmd='spider spi @compare_avgs_to_avgs < inputfilespitmp.txt'
	subprocess.Popen(cmd,shell=True).wait()

	outpath='%s/%s' %(workingfolder,outfile)
	if outpath[0] == '/': 
		outpath=outpath[1:]
	if not os.path.exists('tmpdirspi/aligned_avgs_model.spi'): 
		sys.exit()
	cmd='mv tmpdirspi/aligned_avgs_model.spi /%s' %(outpath)
	subprocess.Popen(cmd,shell=True).wait()

	cmd='mv tmpdirspi/apsh.spi /%s_apsh.spi' %(outpath[:-4])
        subprocess.Popen(cmd,shell=True).wait()

	shutil.rmtree('tmpdirspi')
	removefile('inputfilespitmp.txt')

	resultf = glob.glob("results.spi.*")
	if resultf:
               for f in resultf:
                       os.remove(f)

	return num1,num2,'/%s_apsh.spi'%(outpath[:-4]),box
#===========
def picktxt_to_box(txtpick,boxsize,xlim,ylim,outbox): 

	newfile=open(outbox,'w')
	for line in open(txtpick,'r'): 
		x=float(line.split()[0])-(boxsize/2)
		y=float(line.split()[1])-(boxsize/2)
		if x < 10: 
			continue
		if y < 10: 
			continue
		if x>xlim-10:
			continue
		if y>ylim-10: 
			continue
		newfile.write('%.f\t%.f\t%i\t%i\t-3\n' %(x,y,boxsize,boxsize))

	newfile.close()
	return outbox

#======================
def writeSpiderScript():

	spi='FR\n'
	spi+='?Input normalized stack of class avgs? <stack>	\n'
	spi+='RR [numParts]\n'
	spi+='?Number of class averages to be match?\n'
	spi+='FR \n'
	spi+='?Input OTHER class averages that you would like aligned? <other>\n'
	spi+='RR [numParts2]\n'
	spi+='?Number of class average in second stack?\n'
	spi+='RR [boxSize]\n'
	spi+='?Box size?\n'
	spi+='FR\n'
	spi+='?Output folder? <output>\n'
	spi+='MD\n'
	spi+='SET MP\n'
	spi+='(0)\n'
	spi+='[radius]=[boxSize]/2\n'
	spi+='[imageCenter]=(INT([boxSize]/2)+1)\n'
	spi+='[allowedShift]=0.15*[boxSize]\n'
	spi+='[thetaMin]=0\n'
	spi+='[thetaMax]=90\n'
	spi+='[psiMin]=0\n'
	spi+='[psiMax]=359.9\n'
	spi+='VM\n'
	spi+='mkdir [output]\n'
	spi+='DOC CREATE\n'
	spi+='<output>/angles\n'
	spi+='(3)\n'
	spi+='1-[numParts2]\n'
	spi+='AP SH				; HUGE function.  Is the heart of proj matching.  Compares parts to ref\n'
	spi+='<other>@*****\n'
	spi+='(1-[numParts2])\n'
	spi+='([boxSize]*0.0625,1)		; allow shift up to ~1/20 box size at one pixel intervals (can edit this!)\n'
	spi+='(1,[boxSize]*0.3,1,1)		; allow rings from center to 90% of box to be used for rotational alignment\n'
	spi+='<output>/angles			; angles of the reprojections\n'
	spi+='<stack>@******			; particles to be matched\n'
	spi+='(1-[numParts])			; number of particles\n'
	spi+='*				; can include previous alignment params here to limit search\n'
	spi+='(0.0,0.0)			; amount to restrict\n'
	spi+='(1)				; check mirrors\n'
	spi+='<output>/apsh		; AP SH output text file with all alignment info\n'
	spi+='; [avgCC] - average correlation coefficient. \n'
	spi+='[avgCC]=0.0                        	;This is set to zero, but will be added to at the end of the master loop\n'
	spi+='DO LB5 [particle]=1,[numParts]           \n'
	spi+='UD IC [particle],[psi],[theta],[phi],[matchedProj],[expNum],[cumPsi],[cumX],[cumY],[numSearched],[angChange],[CC],[curPsi],[curX],[curY],[Mir]\n'   
	spi+='<output>/apsh\n'		
	spi+='[curPsi]=-[curPsi]                    		;Since the rotangle is calculated for the particle, its sign\n'
	spi+='IF([Mir].GT.0) GOTO LB61                	;see previous comment\n'
	spi+='MR                            			;mirror image command\n'
	spi+='<other>@{*****[matchedProj]}        ;selects input reference projection\n'
	spi+='_3                            			;temporary output image\n'
	spi+='Y                            			;mirrors across Y axis\n'
	spi+='RT SQ                            		;rotate and shift image       \n'
	spi+='_3                            			;temporary output from the mirroring command previous\n'
	spi+='_2                            			;output rotated reference projection\n'
	spi+='[curPsi]                        		;rotation angle\n'
	spi+='(0,0)                            		;no shifting!\n'
	spi+='goto LB62\n'
	spi+='LB61\n'
	spi+=';for particles that did not align to the mirrored reference, the reference is simply rotated (like above)\n'
	spi+='RT SQ                            		;rotate and shift image\n'
	spi+='<other>@{*****[matchedProj]}        	;input reference projection\n'
	spi+='_2                            			;output  rotated reference projection\n'
	spi+='[curPsi]                        		;rotation angle\n'
	spi+='(0,0)                            		;no shifting!\n'
	spi+='LB62\n'
	spi+=';Now we have the reference projection correctly rotated and mirrored (if necessary) and we will prepare it for\n'
	spi+=';alignment to the particle\n'
	spi+='MA                                		;mask image\n'
	spi+='_2                                		;reference projection\n'
	spi+='_3                                		;output - masked reference projection\n'
	spi+='[radius]                            		;radius of mask\n'
	spi+='D                                		;shape of mask:  D=disc\n'
	spi+='E                                		;choice of mask pixels: E=externally determined\n'
	spi+='(0)                                		;value for mask pixels, the externally determined value\n'
	spi+='[imageCenter],[imageCenter]                	;center coordinates for mask\n'
	spi+=';The masked reference projection will now be padded with extra pixels, making it twice the box size\n'
	spi+='PD                                		;pad image\n'
	spi+='_3                                		;input  masked reference projection\n'
	spi+='_2                                		;output padded, masked, reference projection\n'
	spi+='([boxSize]*2),([boxSize]*2)                    	;dimensions of padding\n'
	spi+='N                                		;Average?  Ans: N  no.  Value will be inputted below\n'
	spi+='(0.000E+00)                            		;value of padded pixels\n'
	spi+='(1,1)                                		;top left coordinates:  not sure why the image is placed in the top\n'
	spi+=';left corner instead of the middle of the padded image\n'
	spi+=';The unfiltered particle will also be padded before alignment to reference\n'
	spi+='CP                                		;copy image\n'
	spi+='<stack>@{******[particle]}                    	;retrieve particle from stack\n'
	spi+='_3                                		;output  unfiltered particle from stack\n'
	spi+='PD                                		;pad image\n'
	spi+='_3                                		;input  unfiltered particle from stack\n'
	spi+='_1                                		;output  padded unfiltered particle from stack\n'
	spi+='([boxSize]*2),([boxSize]*2)                  	;dimensions of padding\n'
	spi+='B                                		;B(order)  background value will be equal to avg. pixel value\n'
	spi+='(1,1)                                		;top left corner coordinates:  still not sure about it.\n'
	spi+='CC N                                		;calculate a normalized cross correlation image\n'
	spi+='_1                               		;input  padded unfiltered particle\n'
	spi+='_2                                		;ref input  padded, masked, reference projection\n'
	spi+='_1                                		;output  replaces the input with the cross correlation image\n'
	spi+=';The cross correlation image contains peak(s) corresponding to the X & Y shifts necessary to align the particle to\n'
	spi+=';the reference.  We need to read off this peak with the X&Y coordinates, and the coordinates will be shifts.\n'
	spi+='WI                                		;Window command\n'
	spi+='_1                                		;input  cross correlation output from CC N\n'
	spi+='_2                                		;output  a larger box size\n'
	spi+='(([allowedShift]*2)+1),(([allowedShift]*2)+1)      ;new window size\n'
	spi+='([boxSize]-[allowedShift]+1),([boxSize]-[allowedShift]+1)  ;new top left coordinates\n'
	spi+=';Peak search (PK) is now performed on the windowed cross correlation image so that the X & Y coordinates correspond\n'
	spi+=';to X & Y shifts.\n'
	spi+='PK [PK_X],[PK_Y],[PK_value],[PK_ratio],[PK_coords_X],[PK_coords_Y],[PK_coords_max]\n'
	spi+='_2                                		;input  windowed cross correlation image\n'
	spi+='(0)                                		;number of peaks\n'
	spi+='; In case there was no peak found for translation, the peak value [PK_coords_max] is 0.\n'
	spi+='; In this case get a value at the origin of the CCF and store in in shift document file.\n'
	spi+='; This value is needed for sorting and calculation of the average correlation coeff.\n'
	spi+='IF([PK_coords_max].EQ.0.0) THEN            	;If the value of the max peak equal zero, get the origin pixel\n'
	spi+='GP [PK_coords_max]                	;get pixel value command\n'
	spi+='_2                            		;input  windowed cross correlation image\n'
	spi+='([allowedShift]+1),([allowedShift]+1)   ;center coordinates\n'
	spi+=';Since there were no peaks found, the particle does not need to be shifted.  Therefore, the unshifted\n'
	spi+=';unfiltered particle is copied into the shifted particle stack\n'
	spi+='CP                            		;copy\n'
	spi+='_3                            		;input  unfiltered particle from stack\n'
	spi+='<output>/parts_shifted_<model>@{******[particle]}            ;output  new shifted stack\n'
	spi+='ELSE                                		;for all other particles that need to be shifted\n'
	spi+='SH F                            		;shift image using fourier interpolation (accurate!)\n'
	spi+='_3                            			;input  unfiltered particle from stack\n'
	spi+='<output>/parts_shifted_<model>@{******[particle]}        	;output placement into new shifted particle stack\n'
	spi+='-[PK_coords_X],-[PK_coords_Y]            	;coordinates from PK used to shift the particles.  Note that the\n'
	spi+=';sign is negative (spider convention)\n'
	spi+='ENDIF\n'
	spi+=';The final step is to start a running count of the correlation coefficients.  Initially, the [avgCC] is equal to\n'
	spi+=';zero, but all of the new CC values are added together to get the total sum of cross correlation values.\n'
	spi+='[avgCC]=[avgCC]+[PK_coords_max]\n'
	spi+=';This will output the shift parameters\n'
	spi+='SD [particle],[PK_coords_X],[PK_coords_Y],[PK_coords_max]\n'
	spi+='<output>/shifts_<model>\n'
	spi+='UD ICE \n'
	spi+='<output>/apsh\n'
	spi+='LB5                                    				;end of the loop for each particle\n'
	spi+=';Calculate the average cross correlation coefficient by dividing the total cross correlation value by the number of\n'
	spi+=';particles\n'
	spi+='[avgCC]=[avgCC]/[numParts]\n'
	spi+='SD -1,[avgCC]                            			;this saves the [avgCC] value as a comment in the last line\n'
	spi+='<output>/shifts_<model>			;of the shifts file using 1 register (column)\n'
	spi+='SD E                                    			;closes the shifts file\n'
	spi+='<output>/shifts_<model>\n'
	spi+='; remove inline files\n'
	spi+='DE\n'
	spi+='_1\n'
	spi+='DE\n'
	spi+='_2\n'
	spi+='DE\n'
	spi+='_3\n'
	spi+='lb2								;end of loop for each reference\n'
	spi+='SD IC NEW\n'
	spi+='incore_model_1\n'
	spi+='4,[numParts]\n'
	spi+='MS\n'
	spi+='_3@\n'
	spi+='[boxSize] [boxSize] 1 \n'
	spi+='([numParts]*2)\n'
	spi+='[count]=1\n'
	spi+='[num]=1\n'
	spi+='do lb1 [particle]=1,[numParts]			;loop over particles\n'
	spi+='UD IC [particle],[shiftX],[shiftY],[CC_Norm1]\n'
	spi+='<output>/shifts_<model>\n'
	spi+='UD IC [particle] [psi1] [theta1] [phi1] [matchedProj1] [expNum] [cumPsi] [cumX] [cumY] [numSearched] [angChange] [CC] [curPsi] [curX] [curY] [Mir]\n'
	spi+='<output>/apsh\n'
	spi+='SD IC [particle] [matchedProj1] [curPsi] [theta1] [phi1]\n'
	spi+='incore_model_1\n'
	spi+='RT SQ\n'
	spi+='<stack>@{*****[particle]}\n'
	spi+='_1\n'
	spi+='[curPsi],1\n'
	spi+='[curX],[curY]\n'
	spi+='IF([Mir].LT.0)THEN\n'
	spi+='MR\n'
	spi+='_1\n'
	spi+='_2\n'
	spi+='Y\n'
	spi+='CP \n'
	spi+='_2\n'
	spi+='_1\n'
	spi+='ENDIF\n'
	spi+='[next]=[count]+1\n'
	spi+='CP\n'
	spi+='_1\n'
	spi+='_3@{*****[count]}\n'
	spi+='CP\n'
	spi+='<other>@{*****[matchedProj1]}\n'
	spi+='_3@{*****[next]}\n'
	spi+='[count]=[count]+2\n'
	spi+='DE\n'
	spi+='_1\n'
	spi+='DE\n'
	spi+='_2\n'
	spi+='lb1\n'
	spi+='UD ICE\n'
	spi+='<output>/shifts_<model>\n'
	spi+='UD ICE\n'
	spi+='<output>/apsh\n'
	spi+='SD IC COPY\n'
	spi+='incore_model_1\n'
	spi+='<output>/model_assignment\n'
	spi+='CP \n'
	spi+='_3@\n'
	spi+='<output>/aligned_avgs_model@\n'
	spi+='END \n'

	if os.path.exists('compare_avgs_to_avgs.spi'):
		os.remove('compare_avgs_to_avgs.spi')
	o1=open('compare_avgs_to_avgs.spi','w')
	o1.write(spi)
	o1.close()
#==============================
def Gctf_estimation(params,outdir,firstmicnum,lastmicnum,gpus,reslim): 

	#Create ctf calculation directory
        os.makedirs('%s/gctf' %(outdir))

	#Check master list to see if micro has been analyzed.  

        #Select random subset of micrographs = large enough to be divided in two groups for validation of template
        allmics=glob.glob('%s/*.mrc' %(params['dir']))
        listMicsToGet=[]
	starmiccounter=1
	for inmic in allmics: 
		if starmiccounter >=firstmicnum: 
			if starmiccounter <= lastmicnum: 
				listMicsToGet.append(inmic)
		starmiccounter=starmiccounter+1
	workingMicList=[]
        counter=0
        maxnum=len(listMicsToGet)
        cs=params['cs']
        if maxnum > len(in_mics):
                maxnum=len(in_mics)
        while counter< maxnum:
                gpucounter=0
                while gpucounter < gpus:
                        if counter+gpucounter <= maxnum: 
				if not os.path.exists('%s/gctf/process%i' %(outdir,gpucounter)):
                                	os.makedirs('%s/gctf/process%i' %(outdir,gpucounter))
				os.symlink(listMicsToGet[counter+gpucounter],'%s/gctf/process%i/%s' %(outdir,gpucounter,listMicsToGet[counter+gpucounter].split('/')[-1]))
                        	workingMicList.append('%s/gctf/process%i/%s' %(outdir,gpucounter,listMicsToGet[counter+gpucounter].split('/')[-1]))
                        gpucounter=gpucounter+1
                counter=counter+gpus
        gpucounter=0
        wait=False
        while gpucounter < gpus:
                microstarfile=estimate_CTF_gctf('%s/gctf/process%i' %(outdir,gpucounter),params['apix'],params['kev'],params['cs'],params['phasePlate'],wait,gpucounter-1)
                gpucounter=gpucounter+1

        #Wait for Gctf to finish the combine star files into folder [outdir]/gctf/
	for mic in workingMicList:
                isdone=0
                while isdone == 0:
                        if os.path.exists('%s.ctf' %(mic[:-4])):
                                isdone=1
	starlist=[]
        gpucounter=0
        while gpucounter <  gpus:
                starlist.append('%s/gctf/process%i/micrographs_all_gctf.star' %(outdir,gpucounter))
                gpucounter=gpucounter+1

        combineStarFiles(starlist,'%s/gctf/micrographs_all_gctf.star' %(outdir))
	goodMicList=CreateGoodMicList('%s/gctf/micrographs_all_gctf.star' %(outdir),reslim,listMicsToGet)
	return goodMicList,'%s/gctf/micrographs_all_gctf.star' %(outdir)

#==============================
def CreateGoodMicList(microstar,reslim,inlist): 
	for entry in inlist: 
		inmic=entry.split('/')[-1]	
		o1=open(microstar,'r')
		for line in o1: 
			if len(line) < 40: 
				continue
			if float(line.split()[11]) > reslim: 
				inlist.remove(entry)
		o1.close()
	return inlist

#==============================
def runGaussianPick2(params,workingdir,initial_micros,negativeStain,partlim,in_mics,mpi,j_threads,gpus,starting_diameter,final_diameter):
	
	#Loop over all diameters that are going to be checked.
	cur_diam=starting_diameter
        while cur_diam <= final_diameter:
		#First, Confirm that the number of starting micrographs will have enough particles for initial classification
		twodcheck=0
		while twodcheck == 0: 
			if os.path.exists('%s/testpicks' %(workingdir)): 
				shutil.rmtree ('%s/testpicks' %(workingdir))
			os.makedirs('%s/testpicks' %(workingdir))
			goodmiclist,outmicstar=Gctf_estimation(params,'%s/testpicks' %(workingdir),1,initial_micros,gpus,params['microresmin'])
			Gauss_pick_micrographs('%s/testpicks' %(workingdir),goodmiclist,params['apix'],binning,multiplier,cur_diam,cur_diam,xdim,overlap_multiplier,boxsizemultiplier,startingthresh,finalthresh,in_mics,mpi,j_threads)
			numParticles=float(subprocess.Popen('cat %s/testpicks/Import/*.box | wc -l' %(workingdir),shell=True, stdout=subprocess.PIPE).stdout.read().strip())	
			if numParticles >= partlim: 
				twodcheck = 1

		print 'made it'
		sys.exit()
			

		cur_diam=cur_diam+diameter_increment

	'''	
	numPartsPicked=0
        initial_mic_adjustcounter=0
	while initial_mic_adjustcounter <= 20: 
		cur_diam=starting_diameter
		while cur_diam <= final_diameter: 
			os.makedirs('%s/testpicks' %(workingdir))
			Gauss_pick_micrographs('%s/testpicks' %(workingdir),initial_micros,1,params['apix'],binning,multiplier,cur_diam,cur_diam,xdim,overlap_multiplier,boxsizemultiplier,startingthresh,finalthresh,in_mics,mpi,j_threads)	
			sys.exit()
			#Gauss_pick_micrographs(workingdir,num_micros,startingMicNum,apix,binning,multiplier,starting_diameter,final_diameter,xdim,overlap_multiplier,boxsizemultiplier,startingthresh,finalthresh,in_mics,mpi,j_threads)
			cur_diam=cur_diam+diameter_increment
	'''
#==============================
def runGaussianPick(params,workingdir,initial_micros,negativeStain,partlim,in_mics,mpi,j_threads,gpus):

	#First, Confirm that the number of starting micrographs will have enough particles for initial classification
	numPartsPicked=0
	initial_mic_adjustcounter=0
	while initial_mic_adjustcounter <= 20: 
		numPartsPicked,in_mics_new=testPicking('%s/testpicks' %(workingdir),initial_micros+initial_mic_adjustcounter,1,negativeStain,params['apix'],mpi,j_threads)
		initial_mic_adjustcounter=initial_mic_adjustcounter+1
		if numPartsPicked>partlim: 
			initial_mic_adjust=initial_mic_adjustcounter
			initial_mic_adjustcounter = 400
	cmd='rm -rf %s/testpicks' %(workingdir)
	subprocess.Popen(cmd,shell=True).wait()

	#First attempt, select beginning, first quarter, halfway, third quarter, and last micrographs for picking, extraction, and 2D classification. 
	totmics=len(glob.glob(in_mics))
	firstQ=round((totmics/4))
	halfway=round(totmics/2)
	thirdQ=round(totmics/4)*3
	lastBit=totmics-initial_micros-1

	#Create empty list: finished_avgs. This list will have the list of all class average output stacks from these runs.
	finished_avgs=[]
	if totmics < 2*initial_micros: 
		avg,diameter,finalbox=pick_run_2D('%s/iteration1' %(workingdir),initial_micros+initial_mic_adjust,1,negativeStain,params['apix'],in_mics,mpi,j_threads,gpus,reslim)
		finished_avgs.append(avg)
		avg,diameter,finalbox=ppick_run_2D('%s/iteration2' %(workingdir),initial_micros+initial_mic_adjust,lastBit,negativeStain,params['apix'],in_mics,mpi,j_threads,gpus,reslim)
		finished_avgs.append(avg)
		maxiteration=2

	if totmics>=2*initial_micros and totmics <4*initial_micros: 
		avg,diameter,finalbox=pick_run_2D('%s/iteration1' %(workingdir),initial_micros+initial_mic_adjust,1,negativeStain,params['apix'],in_mics,mpi,j_threads,gpus,reslim)
		finished_avgs.append(avg)
		avg,diameter,finalbox=pick_run_2D('%s/iteration2' %(workingdir),initial_micros+initial_mic_adjust,halfway,negativeStain,params['apix'],in_mics,mpi,j_threads,gpus,reslim)
		finished_avgs.append(avg)
		avg,diameter,finalbox=pick_run_2D('%s/iteration3' %(workingdir),initial_micros+initial_mic_adjust,lastBit,negativeStain,params['apix'],in_mics,mpi,j_threads,gpus,reslim)	
		finished_avgs.append(avg)
		maxiteration=3

	if totmics >= 4*initial_micros:
		avg,diameter,finalbox=pick_run_2D('%s/iteration1' %(workingdir),initial_micros+initial_mic_adjust,1,negativeStain,params['apix'],in_mics,mpi,j_threads,gpus,reslim)
		finished_avgs.append(avg)
		avg,diameter,finalbox=pick_run_2D('%s/iteration2' %(workingdir),initial_micros+initial_mic_adjust,halfway,negativeStain,params['apix'],in_mics,mpi,j_threads,gpus,reslim)
	     	finished_avgs.append(avg)
		avg,diameter,finalbox=pick_run_2D('%s/iteration3' %(workingdir),initial_micros+initial_mic_adjust,lastBit,negativeStain,params['apix'],in_mics,mpi,j_threads,gpus,reslim)
		finished_avgs.append(avg)
		avg,diameter,finalbox=pick_run_2D('%s/iteration4' %(workingdir),initial_micros+initial_mic_adjust,firstQ,negativeStain,params['apix'],in_mics,mpi,j_threads,gpus,reslim)
		finished_avgs.append(avg)
		#avg,diameter,finalbox=pick_run_2D('%s/iteration5' %(workingdir),initial_micros+initial_mic_adjust,thirdQ,negativeStain,params['apix'],in_mics,mpi,j_threads,gpus)
		#finished_avgs.append(avg)
		maxiteration=4

	#finished_avgs - list of each directory with dog picker picks
	#
	#Get diameters to use for classification - these are in ANgstroms 
	diamlist=[]
	for diam in finished_avgs[0]: 
		diamlist.append(diam.split('diam')[-1])
	
	for diam in diamlist: 
		os.makedirs('%s/diameter%s' %(workingdir,diam))
		
	itercounter=1
	while itercounter <= maxiteration: 
		
		dogpicks

	sys.exit()

	#Now that relion has finished for all of these different runs, all averages will be compared against eachother. Then, for the top match, an FSC (FRC) will be calculated using SPIDER. The resolution of this will help to determine if the average is 'good'.

	#Averages with 1) a match and 2) have an FSC < 30Angstroms will be considered 'good'/
	resthreshmax=resthresh
	resthresh_incr=20
	resmax=resthresh_incr*3
	#Write out spider compare avgs script
	writeSpiderScript()
	while resthreshmax <= resmax+1: 
		compare_avgs_to_avgs(finished_avgs,'%s/DoGPicked_Avgs_Res%i' %(workingdir,resthreshmax),params['apix']*binning,resthreshmax,resthresh_incr,MinMatch)
		spilist=glob.glob('%s/iter*/Class2D/*.spi' %(workingdir))
		for spi in spilist: 
			removefile(spi)
		resthreshmax=resthreshmax+resthresh_incr

	removefile('compare_avgs_to_avgs.spi')

	return '%s/DoGPicked_Avgs_Res30/template_stack.mrcs' %(workingdir),initial_micros+initial_mic_adjust,diameter,finalbox

#==============================
def parallel_relion_image_handler(miclist,outdirectory,options,numthreads):
        mic=0
	maxnum=len(miclist)
        while mic<len(miclist):
                thread=0
		workinglist=[]
                while thread<numthreads:
                        if mic+thread < maxnum: 
				cmd='relion_image_handler --i %s --o %s/%s %s' %(miclist[mic+thread],outdirectory,miclist[mic+thread].split('/')[-1],options)
				subprocess.Popen(cmd,shell=True)
				workinglist.append('%s/%s' %(outdirectory,miclist[mic+thread].split('/')[-1]))
                        thread=thread+1
		numcomplete=0
		while numcomplete <= numthreads: 
			for checkmic in workinglist: 
				if os.path.exists(checkmic): 
					numcomplete=numcomplete+1
			time.sleep(0.5)			
                mic=mic+numthreads

#==============================
def transferStar(instar,outstar,apix,newapix,diradd,reslim):

	ostar=open(outstar,'w')
	instaropen=open(instar,'r')
	for line in instaropen: 
		if len(line) < 40: 
			ostar.write(line)
			continue
		l=line.split()
		detector=float(line.split()[9])*(newapix/apix)
		#mic='%s/%s' %(diradd,line.split()[0].split('/')[-1])
		mic=line.split()[0].split('/')[-1]
		res=float(line.split()[11])
		if res > reslim: 
			continue 
		l[9]=str(detector)
		l[0]=mic
		newline='\t'.join(l)
		ostar.write('%s\n' %(newline))	

	instaropen.close()
	ostar.close()

	return outstar

#==============================
def runGautoMatch(params,templates,in_mics,outdir,negativeStain,nummics,mpi,j_threads,gpus,apix,newapix,diam,finalbox,gautothreshstart,gautothreshfinal,compareAvgs):
	
	#Create working directory
	os.makedirs(outdir)
	#Create ctf calculation directory
	os.makedirs('%s/gctf' %(outdir))
	
	#Select random subset of micrographs = large enough to be divided in two groups for validation of template
	allmics=glob.glob('%s/*.mrc' %(params['dir']))
	listMicsToGet=[random.randint(0,nummics-1) for r in xrange(nummics)] 
	workingMicList=[]
	counter=1
	maxnum=len(listMicsToGet)
	cs=params['cs']
	if maxnum > len(in_mics): 
		maxnum=len(in_mics)
	while counter<= maxnum:
                gpucounter=1
		while gpucounter <= gpus: 
			if not os.path.exists('%s/gctf/process%i' %(outdir,gpucounter)): 
				os.makedirs('%s/gctf/process%i' %(outdir,gpucounter))
			os.symlink(in_mics[counter+gpucounter],'%s/gctf/process%i/%s' %(outdir,gpucounter,in_mics[counter+gpucounter].split('/')[-1]))
			workingMicList.append('%s/gctf/process%i/%s' %(outdir,gpucounter,in_mics[counter+gpucounter].split('/')[-1]))
			gpucounter=gpucounter+1
		counter=counter+gpus
	gpucounter=1
	wait=False
	while gpucounter <= gpus: 
		microstarfile=estimate_CTF_gctf('%s/gctf/process%i' %(outdir,gpucounter),apix,kev,cs,phase,wait,gpucounter-1)
		gpucounter=gpucounter+1
	
	#Wait for Gctf to finish the combine star files into folder [outdir]/gctf/
	for mic in workingMicList: 
		isdone=0
		while isdone == 0: 
			if os.path.exists('%s.ctf' %(mic[:-4])): 
				isdone=1
			
	starlist=[]
	gpucounter=1
	while gpucounter <=  gpus: 
		starlist.append('%s/gctf/process%i/micrographs_all_gctf.star' %(outdir,gpucounter))
		gpucounter=gpucounter+1

	combineStarFiles(starlist,'%s/gctf/micrographs_all_gctf.star' %(outdir)) 
	
	#Bin micrographs
	parallel_relion_image_handler(workingMicList,outdir,' --angpix %f --rescale_angpix %f' %(apix,newapix),(mpi*j_threads)-1)

	#Center template stack if more than average
	removefile('%s_centered.mrcs' %(templates[:-5]))
	cmd='relion_image_handler --i %s --o centered --shift_com' %(templates)
	subprocess.Popen(cmd,shell=True).wait()
 
	#Run gautomatch on all micrographs
	thresh=gautothreshstart
	finalthresh=gautothreshfinal
	maxgputhread=gpus
	allthreads=[]
	if ((finalthresh-thresh)/0.1) < maxgputhread: 
		maxgputhread=((finalthresh-thresh)/0.1)+1
		print maxgputhread
	while thresh < finalthresh: 
		thread=1
		workingthreads=[]
		while thread <= maxgputhread: 
			print thresh 
			os.makedirs('%s/thresh%.2f/' %(outdir,thresh))
			for mic in glob.glob('%s/*.mrc' %(outdir)): 
				os.symlink(mic,'%s/thresh%.2f/%s' %(outdir,thresh,mic.split('/')[-1]))
			cmd='Gautomatch --cc_cutoff %f --apixT %f --apixM %f --diameter %i -T %s %s/thresh%.2f/*.mrc > %s/thresh%.2f/Gautomatch.log' %(thresh,newapix,newapix,diam,templates,outdir,thresh,outdir,thresh)
			if negativeStain == 1: 	
				cmd=cmd+' --dont_invertT'
			subprocess.Popen(cmd,shell=True)
			workingthreads.append('%s/thresh%.2f/' %(outdir,thresh))
			allthreads.append('%s/thresh%.2f/' %(outdir,thresh))
			thread=thread+1
			thresh=thresh+0.1

		tot=0
		for workingfiles in workingthreads: 
			tot=tot+len(glob.glob('%s/*.mrc' %(workingfiles)))
		numcomplete=0	
			
		miclist=glob.glob('%s/*.mrc' %(workingfiles))
		for mic in miclist: 
			isdone=0
			while isdone == 0: 
				if os.path.exists('%s_automatch.box' %(mic[:-4])): 
					isdone=1
	#Extract particles
	for threshdir in allthreads:
		newstar=transferStar('%s/gctf/micrographs_all_gctf.star' %(outdir),'%s/micrographs_all_gctf.star' %(threshdir),apix,apix*binning,'',params['microresmin'])
		diameter=extractRelion(newstar,apix*binning,finalbox,mpi*j_threads,threshdir,'_automatch.box',outdir)

	#Run 2D alignment 
	groupnum=0
        wait=True
        avglist=[]
	if negativeStain == 0:
                ctfcorrection='--ctf'
        if negativeStain == 1: #YES
                ctfcorrection='--ctf --only_flip_phases'

        while groupnum < len(allthreads):
		numParticles=float(subprocess.Popen('cat %s/Extract/particles.star | grep Extract | wc -l' %(allthreads[groupnum]),shell=True, stdout=subprocess.PIPE).stdout.read().strip())
                numClasses=round(numParticles/150) 
		if numClasses < 10: 
			numClasses=10
		if numParticles > 1000: 
			runRelion2D(allthreads[groupnum],'Extract/particles.star',apix*binning,mpi,numClasses,diameter*apix*binning,'--dont_check_norm',j_threads,' --gpu ',ctfcorrection,wait)
			avglist.append('%s/Class2D/run_it025_classes.mrcs' %(allthreads[groupnum]))
		groupnum=groupnum+1

	if compareAvgs is True: 
		#Averages with 1) a match and 2) have an FSC < 30Angstroms will be considered 'good'
        	resthreshmax=resthresh
        	resthresh_incr=20
        	resmax=resthresh_incr*3
        	#Write out spider compare avgs script
       	 	writeSpiderScript()
        	while resthreshmax <= resmax+1:
                	compare_avgs_to_avgs(avglist,'%s/DoGPicked_Avgs_Res%i' %(outdir,resthreshmax),params['apix']*binning,resthreshmax,resthresh_incr,MinMatch)
                	spilist=glob.glob('%s/iter*/Class2D/*.spi' %(workingdir))
                	for spi in spilist:
                        	removefile(spi)
                	resthreshmax=resthreshmax+resthresh_incr
        	removefile('compare_avgs_to_avgs.spi')
		return '%s/DoGPicked_Avgs_Res30/template_stack.mrcs' %(outdir),0.5

	if compareAvgs is False: 
		return avglist[0],999

#==============================
def combineStarFiles(inlist,outfile): 

	staropen=open(outfile,'w')
	#Write header
	isdoneflag=0
	while isdoneflag == 0:
		firstopen=open(inlist[0],'r')
		for line in firstopen:  
			if len(line) <40: 
				staropen.write(line)
				isdoneflag=1
		firstopen.close()

	for infile in inlist:
		readlist=open(infile,'r')
		for line2 in readlist: 
			if len(line2) > 40: 
				staropen.write(line2)
		readlist.close()
	staropen.close()

#==============================
def checkStartingPoint(workingdir): 

	matchinfo=0
	dirlength=glob.glob('%s/*' %(workingdir))
	for dirs in dirlength: 
		if dirs == '%s/gaussPick' %(workingdir): 
			matchinfo=1	
		if dirs == '%s/gautoPick_initial' %(workingdir): 
			matchinfo=2
#		if dirs == '%s/gautoPick' %(workingdir):
#			matchinfo=3

	print 'matchinfo %i' %(matchinfo)
	#Info to return: matchinfo,templateStack,nummics,diam,finalbox

	if matchinfo == 0: 
		shutil.rmtree(workingdir)
		return matchinfo,'','','',''			

	if matchinfo == 1: 
		if os.path.exists('%s/gaussPick/DoGPicked_Avgs_Res30/template_stack_centered.mrcs' %(workingdir)): 
			nummics=len(glob.glob('%s/gaussPick/iteration1/*.mrc' %(workingdir)))
			diam=0
			if os.path.exists('%s/gaussPick/iteration1/Class2D/run_it025_optimiser.star' %(workingdir)): 
				diam=int(subprocess.Popen('cat %s/gaussPick/iteration1/Class2D/run_it025_optimiser.star | grep _rlnParticleDiameter' %(workingdir), shell=True, stdout=subprocess.PIPE).stdout.read().split()[-1])
				finalbox=int(subprocess.Popen('cat %s/gaussPick/iteration1/Class2D/run_it025_optimiser.star | grep _rlnCoarseImageSize' %(workingdir), shell=True, stdout=subprocess.PIPE).stdout.read().split()[-1])
			if diam==0: 
				return 0,'','','','',
			return matchinfo,'%s/gaussPick/DoGPicked_Avgs_Res30/template_stack_centered.mrcs' %(workingdir),nummics,diam,finalbox
		else: 
			shutil.rmtree(workingdir)
			return 0,'','','',''  

	if matchinfo == 2: 
		if os.path.exists('%s/gautoPick_initial/DoGPicked_Avgs_Res30/template_stack.mrcs' %(workingdir)):
			nummics=len(glob.glob('%s/gautoPick_initial/*.mrc' %(workingdir)))
                        diam=0
                        if os.path.exists('%s/gautoPick_initial/Class2D/run_it025_optimiser.star' %(workingdir)):
                                diam=int(subprocess.Popen('cat %s/gautoPick_initial/Class2D/run_it025_optimiser.star | grep _rlnParticleDiameter' %(workingdir), shell=True, stdout=subprocess.PIPE).stdout.read().split()[-1])
                                finalbox=int(subprocess.Popen('cat %s/gautoPick_initial/Class2D/run_it025_optimiser.star | grep _rlnCoarseImageSize' %(workingdir), shell=True, stdout=subprocess.PIPE).stdout.read().split()[-1])
                        if diam==0:
                                return 0,'','','','', 
			return matchinfo,'%s/gautoPick_initial/DoGPicked_Avgs_Res30/template_stack.mrcs' %(workingdir),nummics,diams,finalbox
		else: 
			if os.path.exists('%s/gaussPick/DoGPicked_Avgs_Res30/template_stack_centered.mrcs' %(workingdir)): 
	                        nummics=len(glob.glob('%s/gaussPick/iteration1/*.mrc' %(workingdir)))
                        	diam=0
                        	print '%s/gaussPick/iteration1/Class2D/run_it025_optimiser.star' %(workingdir)
				if os.path.exists('%s/gaussPick/iteration1/Class2D/run_it025_optimiser.star' %(workingdir)):
                                	diam=int(float(subprocess.Popen('cat %s/gaussPick/iteration1/Class2D/run_it025_optimiser.star | grep _rlnParticleDiameter' %(workingdir), shell=True, stdout=subprocess.PIPE).stdout.read().split()[-1]))
                                	finalbox=int(subprocess.Popen('cat %s/gaussPick/iteration1/Class2D/run_it025_optimiser.star | grep _rlnCoarseImageSize' %(workingdir), shell=True, stdout=subprocess.PIPE).stdout.read().split()[-1])
                        	if diam==0: 
                                	return 0,'','','','',
                        	return 1,'%s/gaussPick/DoGPicked_Avgs_Res30/template_stack_centered.mrcs' %(workingdir),nummics,diam,finalbox
	shutil.rmtree(workingdir)
        return 0,'','','',''  

#==============================
if __name__ == "__main__":

        params=setupParserOptions()
        checkConflicts(params)

	#######################	
	#########INPUTS	#######
	#######################
	if params['dir'][-1] == '/':
                params['dir']=params['dir'][:-1]
	if params['outdir'] == 'AutoPick': 
		workingdir='%s/AutoPick' %(params['dir'])
	else: 
		workingdir='%s/'%(params['outdir'])
	in_mics='%s/*.mrc'%(params['dir'])
	gpu=len(subprocess.Popen('lspci | grep NVIDIA', shell=True, stdout=subprocess.PIPE).stdout.read().split('NVIDIA'))-1
	if gpu == 0: 
		print '\nError: No NVIDIA GPUs detected. Exiting\n' 
		sys.exit()
	mpi=gpu+1
	totcores=float(subprocess.Popen('grep -c ^processor /proc/cpuinfo', shell=True, stdout=subprocess.PIPE).stdout.read())
	j_threads=int(math.ceil((totcores-mpi)/mpi))
	if params['dogStartDiam'] == 0: 
		starting_diameter=params['MW']/3#Angstroms
	else: 
		starting_diameter=params['dogStartDiam']
	if params['dogFinDiam'] == 0: 
                final_diameter=params['MW']#Angstroms
        else:
                final_diameter=params['dogFinDiam']
	binningres=5
	binningcheck=binningres/params['apix']
	if binningcheck >2 and binningcheck <=4: 
		binning=4
	if binningcheck <= 2: 
		binning=2
	if binningcheck >4 and binningcheck <=8: 
		binning=8 
	if params['binning'] > 0: 
		binning=params['binning']
	diameter_increment=100
	#Get X & Y dimensions
	xdim=int(subprocess.Popen('relion_image_handler --i %s --stats' %(glob.glob(in_mics)[0]), shell=True, stdout=subprocess.PIPE).stdout.read().split('(x,y,z,n)=')[-1].split('x')[0])
	ydim=int(subprocess.Popen('relion_image_handler --i %s --stats' %(glob.glob(in_mics)[0]), shell=True, stdout=subprocess.PIPE).stdout.read().split('(x,y,z,n)=')[-1].split('x')[1])
	overlap_multiplier=0.8
	startingthresh=params['dogThresh']
	finalthresh=startingthresh #not used right now.
	thresh_incr=0.4
	minnum=20
	reslim=params['microresmin']
	boxsizemultiplier=1.8
	kev=params['kev']
	cs=params['cs']
	phase=params['phasePlate']
	initial_mic_adjust=0
	resthresh=params['resmin']
	MinMatch=3
	if params['debug'] is True: 
		print 'GPUs=%i' %(gpu)
		print 'MPI=%i' %(mpi)
		print 'J_THREADS=%i' %(j_threads)
		print 'XDIM=%i' %(xdim)
		print 'YDIM=%i' %(ydim)
	#Set boundary limits for number of particles required to do first round of classification
	if params['negstain'] is False:  
		initial_micros=20
		if params['partLimit'] == 0: 
			partlim=20000
		else: 
			partlim=params['partLimit']
		if params['invert'] is False: 
                        flipcontrast=True
                if params['invert'] is True:
                        flipcontrast=False
		contrast=0.07
		negativeStain=0
		multiplier=-1
	if params['negstain'] is True:
		initial_micros=10
		if params['partLimit'] == 0: 
                        partlim=5000
                else:
                        partlim=params['partLimit']
		if params['invert'] is False: 
			flipcontrast=False
		if params['invert'] is True: 
			flipcontrast=True
		contrast=0.15
		negativeStain=1
		multiplier=1
	
	#########################################################
	##Check if directory exists and restart run if possible##
	#########################################################
	matchinfo=0
	if params['continue'] is True: 
		if os.path.exists(workingdir): 
			matchinfo,templateStack,nummics,diam,finalbox=checkStartingPoint(workingdir)
			if matchinfo > 0: 
				print '\nContinuing unfinished run in %s...\n' %(workingdir)

	######################################
	##Run initial picking with DogPicker##
	######################################
	if matchinfo == 0: #If matchinfo > 0 this step will be skipped, which assumes this has already run to completion
		if os.path.exists(workingdir): 
			shutil.rmtree(workingdir)
		templateStack,nummics,diam,finalbox=runGaussianPick2(params,'%s/gaussPick' %(workingdir),initial_micros,negativeStain,partlim,in_mics,mpi,j_threads,gpu,starting_diameter,final_diameter) 
		matchinfo=1

	################################################
        ##Run initial template picking with Gautomatch##
        ################################################
	if matchinfo == 1: 
		if os.path.exists('%s/gautoPick_initial' %(workingdir)): 
			shutil.rmtree('%s/gautoPick_initial' %(workingdir))
		templateStack,gautothresh=runGautoMatch(params,templateStack,glob.glob(in_mics),'%s/gautoPick_initial' %(workingdir),negativeStain,nummics*3,mpi,j_threads,gpu,params['apix'],params['apix']*binning,diam,finalbox,0.4,0.8,True)
		matchinfo = 2

	#######################################################
        ##Run template picking with Gautomatch on all micros ##
        #######################################################
	if matchinfo == 2: 
		finishedStack,blank=runGautoMatch(params,templateStack,glob.glob(in_mics),'%s/gautoPick' %(workingdir),negativeStain,len(glob.glob(in_mics)),mpi,j_threads,gpu,params['apix'],params['apix']*binning,diam,finalbox,gautothresh,gautothresh,False)

	##############################
	##Finished output summary ####
	##############################
	fin_numparts=int(subprocess.Popen("cat %sdata.star | grep gauto"%(finishedStack[:-12]), shell=True, stdout=subprocess.PIPE).stdout.read().strip())
	fin_nummics=int(subprocess.Popen("ls %s/*.mrc | wc -l" %(finishedStack.split('Class2D')[0])))
	print 'Finished alignment. Check final class averages here:'
	print '%s' %(finishedStack)
	print '\nTotal number of particles picked = %i' %(fin_numparts)
	print '\nTotal number of micrographs used = %i' %(fin_nummics)

	'''
	To do: 
		-Write restart number to jump into last picking phase with gautomatch
		-Estimate ratio of 'good' to 'bad' particles from guatomatch picking runs, Use to decide on thresh for full picking
		-Kill & restart any failed GPU jobs.
		-add gctf estimation before dogpicking so that bad micros are tossed.
		-check to see if there are at least two runs in order to assess 'good' averages. 
	
	Big to do: 
		-Incorporate workflow: 
			-Remove micrographs that have too few particles (comared to avg. num)
			-Re-run 2D classification without bad images
	
		-Remove micrographs with strong edges -> use feature detection
		-Remove micrographs if too many pixels are >4 sigma (ice)
		-Include waiting feature to return dog picking / initial gauto match if new micros appeared 
	'''	
