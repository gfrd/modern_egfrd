.PHONY: all samples visualize seedgenerator tests egfrd libs support clean

TARGETPATH := ./bin

ifdef MATRIXSIZE
        export DEFMATSIZE := -DMATRIXSIZE=$(MATRIXSIZE)
endif

DEBUG ?= 0
ifeq ($(DEBUG), 1)
	export CXXOPTIONS := -g -D_DEBUG $(DEFMATSIZE) -Wall
endif

all: tests visualize samples seedgenerator

samples: egfrd
	@$(MAKE) --no-print-directory -C Support/RunGfrd
	@cp ./Support/RunGfrd/RunGfrd $(TARGETPATH)

visualize: egfrd
	@$(MAKE) --no-print-directory -C Support/gfrdVisualizer
	@cp ./Support/gfrdVisualizer/gfrdVisualizer $(TARGETPATH)

seedgenerator:   
	@$(MAKE) --no-print-directory -C Support/SeedGenerator
	@cp ./Support/SeedGenerator/SeedGenerator $(TARGETPATH)
   
tests: egfrd
	@$(MAKE) --no-print-directory -C UnitTest/TestGreensFunctions
	@cp ./UnitTest/TestGreensFunctions/TestGreensFunctions $(TARGETPATH)
	@$(MAKE) --no-print-directory -C UnitTest/TestGFRD
	@cp ./UnitTest/TestGFRD/TestGFRD $(TARGETPATH)

egfrd: libs
	@mkdir -p $(TARGETPATH)
	@cp ./Support/genBesselTables/CylindricalBesselTable.bin $(TARGETPATH)
	@cp ./Support/genBesselTables/SphericalBesselTable.bin $(TARGETPATH)
	@cp ./Support/Logger/libLogger.so $(TARGETPATH)
	@cp ./GreensFunctions/libGreensFunctions.so $(TARGETPATH)
	@cp ./eGFRD/libeGFRD.so $(TARGETPATH)
	@echo export LD_LIBRARY_PATH=${PWD}/$(TARGETPATH)

libs: support
	@$(MAKE) --no-print-directory -C GreensFunctions
	@$(MAKE) --no-print-directory -C eGFRD

support:
	@$(MAKE) --no-print-directory -C Support/genBesselTables
	@$(MAKE) --no-print-directory -C Support/Logger

clean:
	@$(MAKE) --no-print-directory -C Support/genBesselTables $@
	@$(MAKE) --no-print-directory -C Support/Logger $@
	@$(MAKE) --no-print-directory -C GreensFunctions $@
	@$(MAKE) --no-print-directory -C eGFRD $@
	@$(MAKE) --no-print-directory -C UnitTest/TestGFRD $@
	@$(MAKE) --no-print-directory -C UnitTest/TestGreensFunctions $@
	@$(MAKE) --no-print-directory -C Support/RunGfrd $@
	@$(MAKE) --no-print-directory -C Support/gfrdVisualizer $@
	@$(MAKE) --no-print-directory -C Support/SeedGenerator $@
	@$(RM) $(TARGETPATH)/*BesselTable.bin $(TARGETPATH)/lib*.so $(TARGETPATH)/Test* $(TARGETPATH)/gfrdVisualizer $(TARGETPATH)/RunGfrd
	@$(RM) -rf $(TARGETPATH)

