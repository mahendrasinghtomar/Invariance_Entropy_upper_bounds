
ifeq ($(OS),Windows_NT)
	COS = win10
else
	UNAME_S := $(shell uname -s)
    ifeq ($(UNAME_S),Linux)
		OS = linux
	endif
	ifeq ($(UNAME_S),Darwin)
		OS = Mac
	endif
endif
# https://stackoverflow.com/questions/714100/os-detecting-makefile

TARGETFOLDER = hyperbolic_v2_minPosts
TARGET = hyperbolic
MINPOST = 0
#DTCONTROL = 1 or 0
DTCONTROL = 1
#cart, logreg, linsvm, oc1
C_METHOD = cart
DETERMINIZER = maxfreq
#venvDtcontrolMST decimal, venvDtcontrolMST_05April2020 float
dtcontrolVersionFolder = venvDtcontrolMST

#
# compiler
#
CC        = g++
MCC       	  = clang++
# CXXFLAGS 		= -g -Wall -Wextra -std=c++11
# CXXFLAGS 		= -Wall -Wextra -std=c++11 -O3 -DNDEBUG -fsanitize=address
CXXFLAGS 		= -Wall -Wextra -std=c++14 -O3 -DNDEBUG

#
# scots 
#
SCOTSROOT		= ../..
SCOTSINC		= -I$(SCOTSROOT)/src2  -I$(SCOTSROOT)/utils/

# lenovo pc
BOOSTINC        = -I/mnt/c/Users/mahen/StorageLocal/boost_1_72_0/
FILESINC		= -I/mnt/c/Users/mahen/GoogleDrive/Currently_working_on_these/lyf_MST_documents/lyf_pc_mac/cpp_packages/SCOTSv0.2/examples/$(TARGETFOLDER)
DTCONTROLREP	= /mnt/c/Users/mahen/StorageLocal/dtcontrol_rep
# XCODEPROJECTFOLDER = /mnt/c/Users/mahen/GoogleDrive/Currently_working_on_these/lyf_MST_documents/lyf_pc_mac/Xcode/SCOTSv02MST/SCOTSv02MST/
MATLAB_linux = "/mnt/c/Program Files/MATLAB/R2019b/bin/matlab.exe"
MATLABFOLDER_linux = $(SCOTSROOT)/examples/$(TARGETFOLDER)


# Mac
MBOOSTINC        = -I/usr/local/Cellar/boost/1.72.0/include
MFILESINC		= -I/Users/mst/GoogleDrive/Currently_working_on_these/lyf_MST_documents/lyf_pc_mac/cpp_packages/SCOTSv0.2/examples/$(TARGETFOLDER)
MDTCONTROLREP		= /Users/mst/dtcontrol_rep
# MXCODEPROJECTFOLDER = /Users/mst/GoogleDrive/Currently_working_on_these/lyf_MST_documents/lyf_pc_mac/Xcode/SCOTSv02MST/SCOTSv02MST
MATLAB_Mac = /Applications/MATLAB_R2018b.app/bin/matlab
MATLABFOLDER_Mac = $(SCOTSROOT)/examples/$(TARGETFOLDER)

# pfaces -CG -d 1 -k gb_fp@../../kernel-pack -cfg hyperbolic_intersection.cfg -p
PFACES_CLOUD_KEY = /mnt/c/Users/mahen/Downloads/AWS_KEY_PFACESDEMO.pem
PFACES_ADDRESS = ubuntu@ec2-34-221-234-199.us-west-2.compute.amazonaws.com
PFACES_LOCAL_FOLDER = /mnt/c/Users/mahen/GoogleDrive/Currently_working_on_these/lyf_MST_documents/lyf_pc_mac/cpp_packages/SCOTSv0.2/examples
PFACES_CLOUD_FOLDER = /home/ubuntu/pFaces-SymbolicControl/ex_gb_fp


# .PHONY: dcdc simulate
#.PHONY: $(TARGET)

#nodtcontrol2: $(TARGET)_Run1_splitA_$(OS) $(TARGET)_Run1_splitB_$(OS) sall_Run2_$(OS) rmatlab

#nodtcontrol3: $(TARGET)_Run1_splitB_$(OS) sall_Run2_$(OS) rmatlab

# For win10-linux shell
#1.run SCOTS, 2.run dtc
#all: sall_Run1_$(OS) dall_$(OS) sall_Run2_$(OS) rmatlab

ifeq ($(DTCONTROL),0)
	USE_DTCONTROL = 
else
	USE_DTCONTROL = dall_$(OS)
endif

#longer length control input
all: $(TARGET)_Run1_splitA_$(OS) $(TARGET)_Run1_splitB_$(OS) $(USE_DTCONTROL) sall_Run2_$(OS) rmatlab

splitB_onwards: $(TARGET)_Run1_splitB_$(OS) $(USE_DTCONTROL) sall_Run2_$(OS) rmatlab

splitA_and_splitB: $(TARGET)_Run1_splitA_$(OS) $(TARGET)_Run1_splitB_$(OS)

splitB: $(TARGET)_Run1_splitB_$(OS)

splitA: $(TARGET)_Run1_splitA_$(OS)

TESTFILE = hyperbolic_fwd
test:
	$(CC) $(CXXFLAGS) $(SCOTSINC) $(FILESINC) $(TESTFILE).cc -o $(TESTFILE)
	./$(TESTFILE)

PFACESfilename = hyperbolic_fwd

pfpull:
	sudo scp -i $(PFACES_CLOUD_KEY) $(PFACES_ADDRESS):$(PFACES_CLOUD_FOLDER)/hyperbolic/$(PFACESfilename).raw /mnt/c/Users/mahen/StorageLocal/pFaces_rawFiles
	
pfpull_log:
	sudo scp -i $(PFACES_CLOUD_KEY) $(PFACES_ADDRESS):/home/ubuntu/pFaces-SymbolicControl/kernel-pack/gb_fp.gpu.cl.build.log $(PFACES_LOCAL_FOLDER)/$(TARGETFOLDER)/
	
pfpush:
	sudo scp -i $(PFACES_CLOUD_KEY) $(PFACES_LOCAL_FOLDER)/$(TARGETFOLDER)/$(PFACESfilename).cfg $(PFACES_ADDRESS):$(PFACES_CLOUD_FOLDER)/hyperbolic/ 

scotspush:
	sudo scp -i $(PFACES_CLOUD_KEY) $(PFACES_LOCAL_FOLDER)/$(TARGETFOLDER)/$(PFACESfilename).cc $(PFACES_ADDRESS):/home/ubuntu/SCOTSv0.2/examples/hyperbolic/ 

# pfpull:
	# scp -i /mnt/c/Users/mahen/Downloads/AWS_KEY_PFACESDEMO.pem ubuntu@ec2-34-220-133-106.us-west-2.compute.amazonaws.com:/home/ubuntu/pFaces-SymbolicControl/README.md $(MATLABFOLDER_linux)
	
#pfpush:
	# scp -i /mnt/c/Users/mahen/Downloads/AWS_KEY_PFACESDEMO.pem /mnt/c/Users/mahen/GoogleDrive/Currently_working_on_these/lyf_MST_documents/lyf_pc_mac/cpp_packages/SCOTSv0.2/examples/aircraftThrusted/aircraftThrusted.cfg ubuntu@ec2-34-220-133-106.us-west-2.compute.amazonaws.com:/home/ubuntu/pFaces-SymbolicControl/ex_gb_fp/aircraftThrusted/ 

#intersection
#allint: sall_Run1_$(OS) 

simu:
	$(CC) $(CXXFLAGS) $(SCOTSINC) $(BOOSTINC) $(FILESINC) $(TARGET)_simulate.cpp -o $(TARGET)_simulate
	./$(TARGET)_simulate

intersec: sall_Run1_$(OS) 

#dtcontrol onwards
dtcon: dall_$(OS) sall_Run2_$(OS) rmatlab

sall_Run1_linux: $(TARGET)_Run1_linux runSEg 

sall_Run2_linux: $(TARGET)_Run2_linux

sall_Run1_Mac: $(TARGET)_Run1_Mac runSEg 

sall_Run2_Mac: $(TARGET)_Run2_Mac

dall_linux: dtcClean_linux copy_linux dtc_linux

dall_Mac: dtcClean_Mac copy_Mac dtc_Mac

runSEg:
	./$(TARGET)

$(TARGET)_Run1_linux: $(TARGET).h
	-rm dtControlClassify.h
	-rm DataMatlab/*.txt
	-rm $(TARGET).scs
	-rm -r decision_trees
	-rm $(TARGET)
	$(CC) $(CXXFLAGS) -DRUN_INDEX=1 -DMINPOST=$(MINPOST) $(SCOTSINC) $(BOOSTINC) $(FILESINC) $(TARGET)_$(OS).cc -o $(TARGET)

$(TARGET)_Run1_splitA_linux: $(TARGET).h
	-rm dtControlClassify.h
	-rm DataMatlab/*.txt
	-rm $(TARGET).scs
	-rm -r decision_trees
	-rm $(TARGET)
	$(CC) $(CXXFLAGS) -DRUN_INDEX=1 -DHYPERBOLIC_STEP=1 $(SCOTSINC) $(BOOSTINC) $(FILESINC) $(TARGET)_$(OS).cc -o $(TARGET)
	./$(TARGET)
	
$(TARGET)_Run1_splitB_linux: $(TARGET).h
	$(CC) $(CXXFLAGS) -DRUN_INDEX=1 -DMINPOST=$(MINPOST) -DHYPERBOLIC_STEP=2 $(SCOTSINC) $(BOOSTINC) $(FILESINC) $(TARGET)_$(OS).cc -o $(TARGET)
	./$(TARGET)

$(TARGET)_Run2_linux: $(TARGET).h
	$(CC) $(CXXFLAGS) -DRUN_INDEX=2 -DHYPERBOLIC_STEP=2 -DDTCONTROL=$(DTCONTROL) $(SCOTSINC) $(BOOSTINC) $(FILESINC) $(TARGET)_$(OS).cc -o $(TARGET)
	./$(TARGET)

$(TARGET)_Run1_Mac: $(TARGET).h
	-rm dtControlClassify.h
	-rm DataMatlab/*.txt
	-rm $(TARGET).scs
	-rm -r decision_trees
	-rm $(TARGET)
	$(MCC) $(CXXFLAGS) -DRUN_INDEX=1 -DMINPOST=$(MINPOST) $(SCOTSINC) $(MBOOSTINC) $(MFILESINC) $(TARGET)_$(OS).cc -o $(TARGET)
	
$(TARGET)_Run1_splitA_Mac: $(TARGET).h
	-rm dtControlClassify.h
	-rm DataMatlab/*.txt
	-rm $(TARGET).scs
	-rm -r decision_trees
	-rm $(TARGET)
	$(MCC) $(CXXFLAGS) -DRUN_INDEX=1 -DHYPERBOLIC_STEP=1 $(SCOTSINC) $(MBOOSTINC) $(MFILESINC) $(TARGET)_$(OS).cc -o $(TARGET)
	./$(TARGET)
	
$(TARGET)_Run1_splitB_Mac: $(TARGET).h
	$(MCC) $(CXXFLAGS) -DRUN_INDEX=1 -DMINPOST=$(MINPOST) -DHYPERBOLIC_STEP=2 $(SCOTSINC) $(MBOOSTINC) $(MFILESINC) $(TARGET)_$(OS).cc -o $(TARGET)
	./$(TARGET)

$(TARGET)_Run2_Mac: $(TARGET).h
	$(MCC) $(CXXFLAGS) -DRUN_INDEX=2 -DHYPERBOLIC_STEP=2 -DDTCONTROL=$(DTCONTROL) $(SCOTSINC) $(MBOOSTINC) $(MFILESINC) $(TARGET)_$(OS).cc -o $(TARGET)
	./$(TARGET)


# simulate: 
#	 $(CC) $(CXXFLAGS) $(SCOTSINC) simulate.cc -o simulate

clean:
	rm $(TARGET)
	
#1. setup virtual env if not done anytime before.
#2. activate virtual env for the terminal window once for the session.
#3. install dtcontrol if not done anytime before.
dtcSetupVirtualEnv_linux:
	virtualenv -p python3 /mnt/c/Users/mahen/StorageLocal/dtcontrol-venv

# https://stackoverflow.com/questions/24736146/how-to-use-virtualenv-in-makefile
VEnvActivateRunDtc_linux:
	(\
	  . /mnt/c/Users/mahen/StorageLocal/dtcontrol-venv/bin/activate; \
	  dtcontrol -i $(DTCONTROLREP)/examples/dcdc/dcdc.scs -m cart -d $(DETERMINIZER) -t 60m --artifact # MaxFreq; \
	)

# 'make' don't activate the virtualenv in the current shell, but opens a new shell, activate and then closes the shell.
# dtcVEnvActivate:
#  . /mnt/c/Users/mahen/StorageLocal/dtcontrol-venv/bin/activate
	
dtcClean_linux:
	-rm benchmark.html benchmark.json $(DTCONTROLREP)/examples/dcdc/dcdc.scs plot.png

dtcClean_Mac:
	-rm benchmark.html benchmark.json $(MDTCONTROLREP)/examples/dcdc.scs plot.png

copy_linux: 
	cp -f $(TARGET).scs $(DTCONTROLREP)/examples/dcdc/dcdc.scs

copy_Mac: 
	cp -f $(TARGET).scs $(MDTCONTROLREP)/examples/dcdc.scs

C_METHOD_DETERMINIZER = -m $(C_METHOD) -d $(DETERMINIZER)

ifeq ($(DETERMINIZER), nil)
	C_METHOD_DETERMINIZER = -m $(C_METHOD)
endif
ifeq ($(C_METHOD), cart)
	ifeq ($(DETERMINIZER), nil)
		CLASSIFY_FOLDERNAME = CART
	endif
	ifeq ($(DETERMINIZER), maxfreq)
		CLASSIFY_FOLDERNAME = MaxFreqDT
	endif
	ifeq ($(DETERMINIZER), minnorm)
		CLASSIFY_FOLDERNAME = MinNormDT
	endif
endif
ifeq ($(C_METHOD), logreg)
	ifeq ($(DETERMINIZER), nil)
		CLASSIFY_FOLDERNAME = LinearClassifierDT-LogisticRegression
	endif
	ifeq ($(DETERMINIZER), maxfreq)
		CLASSIFY_FOLDERNAME = MaxFreq-LinearClassifierDT-LogisticRegression
	endif
	ifeq ($(DETERMINIZER), minnorm)
		CLASSIFY_FOLDERNAME = MinNorm-LinearClassifierDT
	endif
endif
ifeq ($(C_METHOD), linsvm)
	ifeq ($(DETERMINIZER), nil)
		CLASSIFY_FOLDERNAME = LinearClassifierDT-LinearSVC
	endif
	ifeq ($(DETERMINIZER), maxfreq)
		CLASSIFY_FOLDERNAME = MaxFreq-LinearClassifierDT-LinearSVC
	endif
	ifeq ($(DETERMINIZER), minnorm)
		CLASSIFY_FOLDERNAME = MinNorm-LinearClassifierDT
	endif
endif
ifeq ($(C_METHOD), oc1)
	CLASSIFY_FOLDERNAME = OC1
endif

dtc_linux:
	(\
	  . /mnt/c/Users/mahen/StorageLocal/dtcontrol_MSTedited/$(dtcontrolVersionFolder)/bin/activate; \
	  mprof clean; \
	  mprof run dtcontrol -i $(DTCONTROLREP)/examples/dcdc/dcdc.scs $(C_METHOD_DETERMINIZER) -t 300m --artifact; \
	  mprof plot --output=plot.png; \
	  cp -f decision_trees/$(CLASSIFY_FOLDERNAME)/dcdc/$(CLASSIFY_FOLDERNAME).h dtControlClassify.h; \
	#   "/mnt/c/Program Files (x86)/Microsoft/Edge/Application/msedge.exe" C:/Users/mahen/GoogleDrive/Currently_working_on_these/lyf_MST_documents/lyf_pc_mac/cpp_packages/SCOTSv0.2/examples/$(TARGETFOLDER)/benchmark.html; \
	)
	
dtc_Mac:
	(\
	source /Users/mst/Downloads/Local_Storage/dtcontrolMSTedited/$(dtcontrolVersionFolder)/bin/activate; \
	mprof clean; \
	mprof run dtcontrol -i $(MDTCONTROLREP)/examples/dcdc.scs $(C_METHOD_DETERMINIZER) -t 300m --artifact;	\
	mprof plot --output=plot.png; \
	cp -f decision_trees/$(CLASSIFY_FOLDERNAME)/dcdc/$(CLASSIFY_FOLDERNAME).h dtControlClassify.h; \
	)
	
rmatlab:
	$(MATLAB_$(OS)) -nodisplay -nosplash -nodesktop -r "run('$(MATLABFOLDER_$(OS))/SCOTS_dtCon_hyperbolic.m');"

# MATLABFILE = /Users/mst/GoogleDrive/Currently_working_on_these/Matlab/SCOTS_dtCon.m
# . /Applications/MATLAB_R2018b.app/bin/matlab -nodisplay -nosplash -nodesktop -r "run('/Users/mst/GoogleDrive/Currently_working_on_these/Matlab/SCOTS_dtCon.m');exit;"



