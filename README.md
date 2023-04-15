# Invariance_Entropy_upper_bounds
Requires [SCOTSv0.2](https://github.com/mahendrasinghtomar/Invariance_Entropy_upper_bounds/tree/main/SCOTSv0.2) and [dtControl](https://dtcontrol.model.in.tum.de/); copy of both of them is provided here.




# Setup:

Download Boost library: https://sourceforge.net/projects/boost/files/boost/1.72.0/

Download LEMON Graph Library: http://lemon.cs.elte.hu/trac/lemon/wiki/Downloads

Compile the mex function 'comp_hgtm_dtControl': In MATLAB open the folder SCOTSv0.2/mfiles/mexfiles/. Then execute:
	
	mex -R2018a COPTIMFLAGS="-O3 -Oy- -DNDEBUG" comp_hgtm_dtControl.c

In the folder 'dtControl': run the makefile.

# Examples

Go to a folder in: SCOTSv0.2/examples/
	
	Execute the makefile in the selected example folder. 

# Publication

https://www.sciencedirect.com/science/article/abs/pii/S0167691122001724

Mahendra Singh Tomar, C. Kawan, and M. Zamani, ``Numerical over-approximation of invariance entropy via finite abstractions" Systems \& Control Letters, vol.170,pp.105395,2022
