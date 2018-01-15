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
                    help='Pixel size of micrographs (Default=0.9)')
        parser.add_option("--cs",dest="cs",type="float",metavar='Cs',default=2.7,
                    help='Spherical aberration of microscope in mm (Default=2.7)')
	parser.add_option("--kev",dest="kev",type="int",metavar='Kev',default=300,
                    help='Accelerating voltage of microscope (Default=300)')
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
#def checkConflicts(params):


if __name__ == "__main__":

        params=setupParserOptions()

	print 'first test'