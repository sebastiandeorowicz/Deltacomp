all: dc

DC_ROOT_DIR = .
DC_MAIN_DIR = ./DeltaComp
DC_LIBS_DIR = ./libs

CC 	= g++
CFLAGS	= -Wall -O3 -m64 -std=c++14 -mavx -I ~/cfitsio -I ~/CCfits -I DeltaComp/libs/ppmd
CLINK	= -lm -O3 -std=c++11 -lCCfits -lcfitsio -L~/cfitsio -L/usr/local/lib 

.cpp.o:
	$(CC) $(CFLAGS) -c $< -o $@

dc: $(DC_MAIN_DIR)/delta_comp.o \
	$(DC_MAIN_DIR)/DeltaComp.o \
	$(DC_LIBS_DIR)/ppmd/Model.o \
	$(DC_LIBS_DIR)/ppmd/PPMd.o 
	$(CC) $(CLINK) -o $(DC_ROOT_DIR)/$@  $(DC_MAIN_DIR)/delta_comp.o \
	$(DC_MAIN_DIR)/DeltaComp.o \
	$(DC_LIBS_DIR)/ppmd/Model.o \
	$(DC_LIBS_DIR)/ppmd/PPMd.o 

clean:
	-rm $(DC_MAIN_DIR)/*.o
	-rm $(DC_LIBS_DIR)/ppmd/*.o
	-rm DeltaComp
	
