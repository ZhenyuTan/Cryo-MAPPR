#!/usr/bin/env python

import optparse
import glob
import time
import stat
import math
import linecache
import os
import sys
import subprocess
import shutil
import datetime

#=========================
def setupParserOptions():
        parser = optparse.OptionParser()
        parser.set_usage("This program will assess micrograph quality and return a list of good micrographs in .star format")
        parser.add_option("--dir",dest="dir",type="string",metavar="Directory",default='blank',
                    help="Provide directory containing micrographs")
        parser.add_option("--angpix",dest="apix",type="float",metavar='Angpix',default=0.9,
                    help='Pixel size of micrographs in .mrc format (Default=0.9)')
        parser.add_option("--cs",dest="cs",type="float",metavar='Cs',default=2.7,
                    help='Spherical aberration of microscope in mm (Default=2.7)')
	parser.add_option("--kev",dest="kev",type="int",metavar='Kev',default=300,
                    help='Accelerating voltage of microscope (Default=300)')
	parser.add_option("--wildcard",dest='wildcard',type='string',metavar='wildcard',default='',
		    help='Optional: Provide wildcard suffix for input micrographs. Default is none')
        parser.add_option("-v", action="store_true",dest="version",default=False,
                    help="Print version and exit.")
	parser.add_option("-d", action="store_true",dest="debug",default=False,
	            help="debug")
        options,args = parser.parse_args()

        if len(args) > 0:
                parser.error("Unknown commandline options: " +str(args))
        if len(sys.argv) <= 1:
                parser.print_help()
                sys.exit()
        params={}
        for i in parser.option_list:
                if isinstance(i.dest,str):
                        params[i.dest] = getattr(options,i.dest)
        return params

#=============================
def checkConflicts(params):

	#Get number of micrographs in input directory
	numMics=len(glob.glob('%s/*%s.mrc' %(params['dir'],params['wildcard'])))

	#If zero micrographs, throw error & exit
	if numMics == 0: 
		print '\nError: No micrographs found in input directory %s with extension *%s.mrc. Exiting.\n' %(params['dir'],params['wildcard'])
		sys.exit()

#============================
if __name__ == "__main__":
	
	#Version number
	version='0.0.1'
	
	#Get input parameters from the command line
        params=setupParserOptions()

	#Print version and exit
	if params['version'] is True: 
		print '\nVersion=%s\n' %(version)
		sys.exit()

	#Check inputs exist
	checkConflicts(params)

	if params['debug'] is True: 
		print 'finished'
