# Cryo-MAPPR
Cryo-EM micrograph assessment and particle picking routine

## Goal
The goal of this project is to fullly automate micrograph assessment and particle picking through empirical determination of 'good' parameters for particle picking. This will entail shape based picking (DoG picking) and then later reference-based picking. 

## Installation

Software dependencies: 
* mrcfile - python package to read MRC Files ($ pip install mrcfile --user)
* CTFFIND4 - CTF estimation program (http://grigoriefflab.janelia.org/ctffind4)
* scipy & numpy ($ pip install numpy scipy --user) 
* Python image library ($ pip install pillow --user)
* mahota ($ pip install mahota --user)

Cryo-MAPPR requires that the Github repo be added into the user PYTHONPATH. For bash shells, include this in your .bashrc: 

<pre>export PYTHONPATH=/path/to/Cryo-MAPPR/lib:$PYTHONPATH</pre>

## Usage

To available options, run run_mapper.py without any inputs: 
<pre>
$ Cryo-MAPPR/run_mappr.py 
Usage: This program will assess micrograph quality and return a list of good micrographs in .star format

Options:
  -h, --help           show this help message and exit
  --dir=Directory      Provide directory containing micrographs
  --diam=Diam          Particle diameter in Angstroms (Default=115A)
  --angpix=Angpix      Pixel size of micrographs in .mrc format (Default=0.9)
  --cs=Cs              Spherical aberration of microscope in mm (Default=2.7)
  --kev=Kev            Accelerating voltage of microscope (Default=300)
  --wildcard=wildcard  Optional: Provide wildcard suffix for input
                       micrographs. Default is none
  -v                   Print version and exit.
  -d                   debug</pre>

**Note**: Program only reads MRC image files and users can specify wildcard flags to select subsets of micrographs within a given directory.

## Outputs

Currently, the program will output two files into the micrograph directory: 

<pre>good_micrographs.txt
bad_micrographs.txt</pre>

## Development notes
 
The majority of the code must be written in the Cryo-MAPPR/lib python library and then imported, in order to run. 

For example, to run the defition 'testing' from lib/input_check.py: 

<pre>>>> import input_check
>>> input_check.testing()</pre>

For all debugging, print outputs using the following: 

<pre>if params['debug'] is True:
	print 'text here'</pre>


