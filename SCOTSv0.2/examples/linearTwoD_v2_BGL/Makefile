#--------------------------------------------------------------------#
#$$$$ AdjustPath     $$$$$#
#cart, logreg, linsvm, oc1
C_METHOD = cart
DETERMINIZER = maxfreq

BOOSTINC        = -I/AdjustPath/boost_1_72_0/
MATLAB = "/mnt/c/Program Files/MATLAB/R2019b/bin/matlab.exe"
#$$$$ end Adjust $$$$#
#--------------------------------------------------------------------#


#
# compiler
#
CC        = g++
MCC       	  = clang++
# CXXFLAGS 		= -g -Wall -Wextra -std=c++11
CXXFLAGS 		= -Wall -Wextra -std=c++14 -O3 -DNDEBUG


TARGETFOLDER = linearTwoD_v2_BGL
TARGET = linearTwoD_dtCon

#
# scots 
#
SCOTSROOT		= ../..
SCOTSINC		= -I$(SCOTSROOT)/src2  -I$(SCOTSROOT)/utils/

FILESINC		= -I$(SCOTSROOT)/examples/$(TARGETFOLDER)
DTCONTROLREP	= ../../../dtControl


all: sall_Run1

dot: 
	dot -Tpng gvz_file.txt -o gvz_file.png
	cmd.exe /C start gvz_file.png

all2: sall_Run1_$(OS) dall_$(OS) sall_Run2_$(OS) 

dtcon: dall_$(OS) sall_Run2_$(OS)

sall_Run1: $(TARGET)_Run1 runSEg 

sall_Run2_linux: $(TARGET)_Run2_linux

sall_Run1_Mac: $(TARGET)_Run1_Mac runSEg 

sall_Run2_Mac: $(TARGET)_Run2_Mac 

dall_linux: dtcClean_linux copy_linux dtc_linux

dall_Mac: dtcClean_Mac copy_Mac dtc_Mac

runSEg:
	./$(TARGET)

$(TARGET)_Run1: $(TARGET).cc
	-rm dtControlClassify.h
	-rm DataMatlab/*.txt
	-rm $(TARGET).scs
	-rm -r decision_trees
	-rm $(TARGET)
	$(CC) $(CXXFLAGS) -DRUN_INDEX=1 $(SCOTSINC) $(BOOSTINC) $(FILESINC) $(TARGET).cc -o $(TARGET)

$(TARGET)_Run2_linux: $(TARGET).h
	$(CC) $(CXXFLAGS) -DRUN_INDEX=2 $(SCOTSINC) $(BOOSTINC) $(FILESINC) $(TARGET)_$(OS).cc -o $(TARGET)
	./$(TARGET)

$(TARGET)_Run1_Mac: $(TARGET).h
	-rm dtControlClassify.h
	-rm DataMatlab/*.txt
	-rm $(TARGET).scs
	-rm -r decision_trees
	-rm $(TARGET)
	$(MCC) $(CXXFLAGS) -DRUN_INDEX=1 $(SCOTSINC) $(MBOOSTINC) $(MFILESINC) $(TARGET)_$(OS).cc -o $(TARGET)

$(TARGET)_Run2_Mac: $(TARGET).h
	$(MCC) $(CXXFLAGS) -DRUN_INDEX=2 $(SCOTSINC) $(MBOOSTINC) $(MFILESINC) $(TARGET)_$(OS).cc -o $(TARGET)
	./$(TARGET)


clean:
	rm $(TARGET)
	




