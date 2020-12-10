# Invariance_Entropy_upper_bounds
SCOTS: https://gitlab.lrz.de/matthias/SCOTSv0.2 

dtControl: https://dtcontrol.model.in.tum.de/

# Setup:

Need to adjust paths of files and folders.

Compile the mex function 'comp_hgtm_dtControl': In MATLAB open the folder SCOTSv0.2/mfiles/mexfiles/. Then execute:
	
	mex -R2018a COPTIMFLAGS="-O3 -Oy- -DNDEBUG" comp_hgtm_dtControl.c

In the folder 'dtControl': run the makefile.

Linear 2 dimensional example: SCOTSv0.2/examples/linearTwoD_v2/
