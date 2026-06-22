REQUIRED_VAR = APPTAINER_HPC

#Ensure path to HPC apptainer library is set
ifndef $(REQUIRED_VAR)
$(error WARNING: $(REQUIRED_VAR) is not set. Please export it before running make.)
endif

all: localintegrationtest

localintegrationtest:
	pixi run -e test-local test-local
