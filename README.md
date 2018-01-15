# Cryo-MAPPR
Cryo-EM micrograph assessment and particle picking routine

## Goal
The goal of this project is to fullly automate micrograph assessment and particle picking through empirical determination of 'good' parameters for particle picking. This will entail shape based picking (DoG picking) and then later reference-based picking. 

## Installation

Software dependencies: 
* CTFFIND4 

Cryo-MAPPR requires that the Github repo be added into the user PYTHONPATH. For bash shells, include this in your .bashrc: 

export PYTHONPATH=/path/to/Cryo-MAPPR/lib:$PYTHONPATH

## Development notes
 
The majority of the code must be written in the Cryo-MAPPR/lib python library and then imported, in order to run. 

For example, to run the defition 'testing' from lib/input_check.py: 

>> import input_check
>> input_check.testing()



