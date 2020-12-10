
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


dtcontrolVersionFolder = venvDtcontrol

all: instEvn envActiInsDtc

clean:
	rm -r $(dtcontrolVersionFolder)

instEvn:
	virtualenv -p python3 $(dtcontrolVersionFolder)

envActiInsDtc: envActiInsDtc_$(OS)

envActiInsDtc_linux:
	(\
	. $(dtcontrolVersionFolder)/bin/activate; \
	pip install dtcontrol_tum-1.0.0rc1-py3-none-any.whl; \
	)

envActiInsDtc_Mac:
	source $(dtcontrolVersionFolder)/bin/activate;
	pip install dtcontrol_tum-1.0.0rc1-py3-none-any.whl

installDtc:
	pip install dtcontrol_tum-1.0.0rc1-py3-none-any.whl
