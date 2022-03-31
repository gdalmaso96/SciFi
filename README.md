## SciFi
SciFi  simulation

# Macros
 - signalsLite.C: macro to process MC output and mimic DAQ
 - merge.C: macro to merge files
 - parallel.py: parallelization script. To be improved: so far there is no queue, the next processes are started only after all the previous batch has been terminated. With this version is better to run geant4 macros that have similar processing time 

# Installation:

git clone this repo and create a build directory.

<code> cd path-to-build

<code> cmake path-to-repo

<code> make -jN

With N <= number of machine cores

