# Cryo-MAPPR
Cryo-EM micrograph assessment and particle picking routine

## Goal
The goal of this project is to fullly automate micrograph assessment and particle picking through empirical determination of 'good' parameters for particle picking. This will entail shape based picking (DoG picking) and then later reference-based picking. 

## Installation

Software dependencies: 
* mrcfile - python package to read MRC Files ($ pip install mrcfile --user)
* CTFFIND4 - CTF estimation program (http://grigoriefflab.janelia.org/ctffind4)
* scipy & numpy 

Cryo-MAPPR requires that the Github repo be added into the user PYTHONPATH. For bash shells, include this in your .bashrc: 

<pre>export PYTHONPATH=/path/to/Cryo-MAPPR/lib:$PYTHONPATH</pre>

## Development notes
 
The majority of the code must be written in the Cryo-MAPPR/lib python library and then imported, in order to run. 

For example, to run the defition 'testing' from lib/input_check.py: 

<pre>>>> import input_check
>>> input_check.testing()</pre>

For all debugging, print outputs using the following: 

<pre>if params['debug'] is True:
	print 'text here'</pre>


