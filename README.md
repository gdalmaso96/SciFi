## SciFi
SciFi  simulation

# Macros
 - signalsLite.C: macro to process MC output and mimic DAQ
 - merge.C: macro to merge files
 - parallel.py: parallelization script. To be improved: so far there is no queue, the next processes are started only after all the previous batch has been terminated. With this version is better to run geant4 macros that have similar processing time 

# Installation:

git clone this repo and create a build directory.

<code> cd path-to-build </code>

<code> cmake path-to-repo </code>

<code> make -jN </code>

With N <= number of machine cores

# Output signals
To see the signals you must use signals.C macro:

<code> bash>> root </code>

<code> root-bash>> .L signals.C </code>

<code> root-bash>> signals() </code> 

or 

<code> signals(path-to-file) </code>

or 

<code> signals(path-to-file, path-to-output-file<no file extension>) </code> 

or 

<code> signals(path-to-file, threashold, path-to-output-file<no file extension>)</code> (threashold is in number of pixel activated per channel)

These commands will produce a out.txt file. Here you have the number of signals per channel, the proton current (always at max 0.22) and the run duration. To extract beam information from datas you can use readgio.cpp macro:

<code> bash>> root </code>

<code> root-bash>> .L readgio.cpp </code>

<code> root-bash>> readMatrix("out.txt") </code> (it requires a "fig" directory to exist in the build folder)

