# Invariance_Entropy_upper_bounds
SCOTS from: https://gitlab.lrz.de/matthias/SCOTSv0.2 

dtControl from: https://dtcontrol.model.in.tum.de/

# Setup:

Download Boost library: https://sourceforge.net/projects/boost/files/boost/1.72.0/

Download LEMON Graph Library: http://lemon.cs.elte.hu/trac/lemon/wiki/Downloads

Compile the mex function 'comp_hgtm_dtControl': In MATLAB open the folder SCOTSv0.2/mfiles/mexfiles/. Then execute:
	
	mex -R2018a COPTIMFLAGS="-O3 -Oy- -DNDEBUG" comp_hgtm_dtControl.c

In the folder 'dtControl': run the makefile.

# Examples

Go to a folder in: SCOTSv0.2/examples/
	
	Execute the makefile in the selected example folder. 

