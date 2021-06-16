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

TARGETFOLDER = hyperbolic_v2_BGL
TARGET = hyperbolic
MINPOST = 0
#cart, logreg, linsvm, oc1
C_METHOD = linsvm
DETERMINIZER = maxfreq

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
BOOSTINC        = -I/mnt/c/Users/..../StorageLocal/boost_1_72_0/
FILESINC		= -I/mnt/c/Users/..../SCOTSv0.2/examples/$(TARGETFOLDER)
DTCONTROLREP	= /mnt/c/Users/..../dtcontrol_rep

MATLAB_linux = "/mnt/c/Program Files/MATLAB/R2019b/bin/matlab.exe"
MATLABFOLDER_linux = $(SCOTSROOT)/examples/$(TARGETFOLDER)
LEMONROOT = /mnt/c/Users/..../lemon-1.3.1_wsl/lemon-1.3.1
LEMONINC = -I $(LEMONROOT)/build -I $(LEMONROOT)  -L $(LEMONROOT)/build/lemon -lemon


# Mac
MBOOSTINC        = -I/usr/local/Cellar/boost/1.72.0/include
MFILESINC		= -I/Users/..../SCOTSv0.2/examples/$(TARGETFOLDER)
MDTCONTROLREP		= /Users/..../dtcontrol_rep

MATLAB_Mac = /Applications/MATLAB_R2018b.app/bin/matlab
MATLABFOLDER_Mac = $(SCOTSROOT)/examples/$(TARGETFOLDER)


mew_mcm=0
# To compute max edge-wt and cycle mean without dtControl.
all: $(TARGET)_splitB_$(OS)

all_o: $(TARGET)_Run1_splitA_$(OS) $(TARGET)_Run1_splitB_$(OS)

# Use dtControl with max edge-wt and cycle mean.
a1d2: $(TARGET)_Run1_splitB_$(OS) dall_$(OS) sall_Run2_$(OS) 

a1: $(TARGET)_Run1_splitB_$(OS)

a1d: $(TARGET)_Run1_splitB_$(OS) dall_$(OS)

a2: sall_Run2_$(OS)

#all3: sall_Run1_$(OS)

dot: 
	dot -Tpng gvz_file.txt -o gvz_file.png
	cmd.exe /C start gvz_file.png

#dtcontrol onwards
dtcon: dall_$(OS) sall_Run2_$(OS)

sall_Run1_linux: $(TARGET)_Run1_linux runSEg 

sall_Run2_linux: $(TARGET)_Run2_linux

sall_Run1_Mac: $(TARGET)_Run1_Mac runSEg 

sall_Run2_Mac: $(TARGET)_Run2_Mac

dall_linux: dtcClean_linux copy_linux dtc_linux

dall_Mac: dtcClean_Mac copy_Mac dtc_Mac

run:
	./$(TARGET)

$(TARGET)_Run1_linux: $(TARGET).h
	-rm dtControlClassify.h
	-rm DataMatlab/*.txt
	-rm $(TARGET).scs
	-rm -r decision_trees
	-rm $(TARGET)
	$(CC) $(CXXFLAGS) -DRUN_INDEX=1 $(SCOTSINC) $(BOOSTINC) $(FILESINC) $(TARGET)_$(OS).cc -o $(TARGET)
	
$(TARGET)_Run1_splitA_linux: $(TARGET).h
	-rm dtControlClassify.h
	-rm DataMatlab/*.txt
	-rm $(TARGET).scs
	-rm -r decision_trees
	-rm $(TARGET)
	$(CC) $(CXXFLAGS) -DHYPERBOLIC_STEP=1 $(SCOTSINC) $(BOOSTINC) $(FILESINC) $(TARGET)_$(OS).cc -o $(TARGET)
	./$(TARGET)

$(TARGET)_splitB_linux: $(TARGET).h
	$(CC) $(CXXFLAGS) -DHYPERBOLIC_STEP=2 -DDTCONTROL=0 -DRUN_INDEX=1 -DMAXEDGEWT_MCM=$(mew_mcm) $(SCOTSINC) $(BOOSTINC) $(FILESINC) $(TARGET)_$(OS).cc -o $(TARGET) $(LEMONINC)
	./$(TARGET)
	
$(TARGET)_Run1_splitB_linux: $(TARGET).h
	-rm dtControlClassify.h
	-rm -r decision_trees
	$(CC) $(CXXFLAGS) -DHYPERBOLIC_STEP=2 -DDTCONTROL=1 -DRUN_INDEX=1 -DMINPOST=$(MINPOST) $(SCOTSINC) $(BOOSTINC) $(FILESINC) $(TARGET)_$(OS).cc -o $(TARGET)
	./$(TARGET)


$(TARGET)_Run2_linux: $(TARGET).h
	$(CC) $(CXXFLAGS) -DDTCONTROL=1 -DRUN_INDEX=2 -DMAXEDGEWT_MCM=$(mew_mcm) $(SCOTSINC) $(BOOSTINC) $(FILESINC) $(TARGET)_$(OS).cc -o $(TARGET) $(LEMONINC)
	./$(TARGET)

$(TARGET)_Run1_Mac: $(TARGET).h
	-rm dtControlClassify.h
	-rm DataMatlab/*.txt
	-rm $(TARGET).scs
	-rm -r decision_trees
	-rm $(TARGET)
	$(MCC) $(CXXFLAGS) -DRUN_INDEX=1 $(SCOTSINC) $(MBOOSTINC) $(MFILESINC) $(TARGET)_$(OS).cc -o $(TARGET)
	
$(TARGET)_Run1_splitA_Mac: $(TARGET).h
	-rm dtControlClassify.h
	-rm DataMatlab/*.txt
	-rm $(TARGET).scs
	-rm -r decision_trees
	-rm $(TARGET)
	$(MCC) $(CXXFLAGS) -DHYPERBOLIC_STEP=1 $(SCOTSINC) $(MBOOSTINC) $(MFILESINC) $(TARGET)_$(OS).cc -o $(TARGET)
	./$(TARGET)
	
$(TARGET)_Run1_splitB_Mac: $(TARGET).h
	$(MCC) $(CXXFLAGS) -DHYPERBOLIC_STEP=2 $(SCOTSINC) $(MBOOSTINC) $(MFILESINC) $(TARGET)_$(OS).cc -o $(TARGET)
	./$(TARGET)

$(TARGET)_Run2_Mac: $(TARGET).h
	$(MCC) $(CXXFLAGS) -DRUN_INDEX=2 $(SCOTSINC) $(MBOOSTINC) $(MFILESINC) $(TARGET)_$(OS).cc -o $(TARGET)
	./$(TARGET)


clean:
	rm $(TARGET)
	
	
dtcClean_linux:
	-rm benchmark.html benchmark.json $(DTCONTROLREP)/examples/dcdc/dcdc.scs

dtcClean_Mac:
	-rm benchmark.html benchmark.json $(MDTCONTROLREP)/examples/dcdc.scs

copy_linux: 
	cp -f $(TARGETFOLDER).scs $(DTCONTROLREP)/examples/dcdc/dcdc.scs

copy_Mac: 
	cp -f $(TARGETFOLDER).scs $(MDTCONTROLREP)/examples/dcdc.scs

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
	  . /mnt/c/Users/..../venvDtcontrolMST/bin/activate; \
	  dtcontrol -i $(DTCONTROLREP)/examples/dcdc/dcdc.scs $(C_METHOD_DETERMINIZER) -t 60m --artifact; \
	  cp -f decision_trees/$(CLASSIFY_FOLDERNAME)/dcdc/$(CLASSIFY_FOLDERNAME).h dtControlClassify.h; \
	)
	
dtc_Mac:
	source /Users/..../venvDtcontrolMST/bin/activate;
	dtcontrol -i $(MDTCONTROLREP)/examples/dcdc.scs $(C_METHOD_DETERMINIZER) -t 60m --artifact;	
	cp -f decision_trees/$(CLASSIFY_FOLDERNAME)/dcdc/$(CLASSIFY_FOLDERNAME).h dtControlClassify.h
	# open -a "Google Chrome" benchmark.html
	
rmatlab:
	$(MATLAB_$(OS)) -nodisplay -nosplash -nodesktop -r "run('$(MATLABFOLDER_$(OS))/SCOTS_dtCon_hyperbolic.m');"




