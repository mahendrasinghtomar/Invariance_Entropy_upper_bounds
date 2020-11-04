#--------------------------------------------------------------------#
#$$$$ Adjust     $$$$$#
#cart, logreg, linsvm, oc1
C_METHOD = cart
DETERMINIZER = maxfreq
MINPOST = 0

BOOSTINC        = -I/AdjustPath/boost_1_72_0/
MATLAB_linux = "/mnt/c/Program Files/MATLAB/R2019b/bin/matlab.exe"
#$$$$ end Adjust $$$$#
#--------------------------------------------------------------------#


#
# compiler
#
CC        = g++
#CC       	  = clang++
CXXFLAGS 		= -Wall -Wextra -std=c++14 -O3 -DNDEBUG


TARGETFOLDER = linearTwoD_v2
TARGET = linearTwoD_dtCon
TARGETM = linearTwoD

#
# scots 
#
SCOTSROOT		= ../..
SCOTSINC		= -I$(SCOTSROOT)/src2  -I$(SCOTSROOT)/utils/

FILESINC		= -I$(SCOTSROOT)/examples/$(TARGETFOLDER)
DTCONTROLREP	= /AdjustPath/dtcontrol_rep
MATLABFOLDER_linux = $(SCOTSROOT)/examples/$(TARGETFOLDER)


all: sall_Run1 dall sall_Run2 rmatlab

dtcon: dall sall_Run2 rmatlab

sall_Run1: $(TARGET)_Run1 runSEg 

sall_Run2: $(TARGET)_Run2

dall_linux: dtcClean copy dtc

runSEg:
	./$(TARGET)

$(TARGET)_Run1: $(TARGET).h
	-rm dtControlClassify.h
	-rm DataMatlab/*.txt
	-rm $(TARGET).scs
	-rm -r decision_trees
	-rm $(TARGET)
	$(CC) $(CXXFLAGS) -DRUN_INDEX=1 -DMINPOST=$(MINPOST) $(SCOTSINC) $(BOOSTINC) $(FILESINC) $(TARGET).cc -o $(TARGET)

$(TARGET)_Run2: $(TARGET).h
	$(CC) $(CXXFLAGS) -DRUN_INDEX=2 $(SCOTSINC) $(BOOSTINC) $(FILESINC) $(TARGET).cc -o $(TARGET)
	./$(TARGET)

clean:
	rm $(TARGET)
	
dtcSetupVirtualEnv_linux:
	virtualenv -p python3 /mnt/c/Users/mahen/StorageLocal/dtcontrol-venv

VEnvActivateRunDtc_linux:
	(\
	  . /AdjustPath/dtcontrol-venv/bin/activate; \
	  dtcontrol -i $(DTCONTROLREP)/examples/dcdc/dcdc.scs -m cart -d $(DETERMINIZER) -t 60m --artifact # MaxFreq; \
	)

dtcClean:
	-rm benchmark.html benchmark.json $(DTCONTROLREP)/examples/dcdc/dcdc.scs plot.png

copy: 
	cp -f $(TARGET).scs $(DTCONTROLREP)/examples/dcdc/dcdc.scs

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
	  . /AdjustPath/venvDtcontrolMST/bin/activate; \
	  mprof clean; \
	  mprof run dtcontrol -i $(DTCONTROLREP)/examples/dcdc/dcdc.scs $(C_METHOD_DETERMINIZER) -t 60m --artifact # MaxFreq; \
	  mprof plot --output=plot.png; \
	  cp -f decision_trees/$(CLASSIFY_FOLDERNAME)/dcdc/$(CLASSIFY_FOLDERNAME).h dtControlClassify.h; \
	)
	
rmatlab:
	$(MATLAB_$(OS)) -nodisplay -nosplash -nodesktop -r "run('$(MATLABFOLDER_$(OS))/SCOTS_dtCon_$(TARGETM).m');"


