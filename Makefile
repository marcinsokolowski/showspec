# CC_FLAGS=-D_SUPERFAST_OPTIMIZED_
# CC_FLAGS=

# HDF5 file library :
HDF5_INCLUDE_DIR=/usr/include/hdf5/serial/
HDF5_LIB_DIR=/usr/lib/x86_64-linux-gnu/hdf5/serial/

# INCLUDES=-D_UNIX -I/opt/bighorns/ext/dload/HealPix/src/cxx/generic_gcc/include/ -I/opt/bighorns/ext/dload/HealPix/src/cxx/alice/ $(CC_FLAGS) -I./
INCLUDES=-D_UNIX -I/opt/bighorns/ext/dload/HealPix/src/cxx/generic_gcc/include/ $(CC_FLAGS) -I./ -I $(HDF5_INCLUDE_DIR)
# LIBS=-L/opt/bighorns/ext/dload/HealPix/src/cxx/generic_gcc/lib/ -lhealpix_cxx  -lcxxsupport  -lnova -lpsht /usr/lib/libgomp.so.1 # /usr/lib/libgomp.so.1 -lpsht
LIBS=-L/opt/bighorns/ext/dload/HealPix/src/cxx/generic_gcc/lib/ -lhealpix_cxx  -lcxxsupport  -lnova 
COMMON_LIBS=-lcfitsio
PROGNAME=showspec_standalone
PROGNAME_MWA=showspec_standalone_mwa
PREPROC=-D_NO_ROOT_

OBJECTS=antenna_pattern.o user_defined_antenna_pattern.o skymap_cache.o ce300e.o showspec_math.o MollweideSkyMap.o \
        bg_components.o bg_array.o bg_fits.o bg_geo.o libnova_interface.o bg_total_power.o bg_units.o bg_globals.o  bg_date.o weather_station.o bg_bedlam.o cvalue_vector.o bg_stat.o \
	myparser.o mystring.o basestring.o myfile.o mysafekeytab.o cexcp.o mycmnglobals.o mydate.o mystrtable.o mykeytab.o myfract.o myenv.o basestructs.o random.o mylock.o mydate2.o \
	myfits.o mathfunc.o mwa_beam_interface.o # common libraries 

C_OBJECTS=gal2hor.o

CPP_COMP=g++
C_COMP=gcc
# OPT=-O3

all : $(PROGNAME) $(PROGNAME_MWA)
#	cp showspec_standalone $(BIGHORNS)/bin/ 
#	cp showspec_standalone_test $(BIGHORNS)/bin/	
#	cp showspec_standalone_intnew $(BIGHORNS)/bin/	
#	cp libnovatest $(BIGHORNS)/bin/
#	cp get_ant_gain $(BIGHORNS)/bin/
#	cp get_ant_gain_azh $(BIGHORNS)/bin/
#	cp get_ant_gain_bulk $(BIGHORNS)/bin/
#	cp test_ant_pattern $(BIGHORNS)/bin/
#	cp datVStxt $(BIGHORNS)/bin/
#	cp datVStxt! $(BIGHORNS)/bin/

$(PROGNAME) : showmap.cpp ce300e.cpp bg_components.cpp antenna_pattern.cpp skymap_cache.cpp showmap_test.cpp showmap_intnew.cpp $(OBJECTS) $(C_OBJECTS)
	@echo
	@echo "Compiling $(PROGNAME)"
	g++  $(OPT) mwa_beam_interface.cpp -c      
	g++ $(OPT) showmap.cpp $(C_OBJECTS) $(OBJECTS) -o $(PROGNAME) $(INCLUDES) $(LIBS) $(COMMON_LIBS)

# separate compilation to include MWA beam model as it requires more external libraries which might not be available everywhere 
$(PROGNAME_MWA) : showmap.cpp ce300e.cpp bg_components.cpp antenna_pattern.cpp skymap_cache.cpp showmap_test.cpp showmap_intnew.cpp mwa_beam_interface.cpp beam2016implementation.cpp $(OBJECTS) $(C_OBJECTS) 
	@echo
	@echo "Compiling $(PROGNAME_MWA)"	        
	g++  $(OPT) beam2016implementation.cpp -c -I $(HDF5_INCLUDE_DIR)
	g++  $(OPT) system.cpp -c -std=gnu++11
	g++  $(OPT) mwa_beam_interface.cpp -D_MWA_2016_BEAM_MODEL_ -c -I $(HDF5_INCLUDE_DIR)
	g++ $(OPT) showmap.cpp $(C_OBJECTS) $(OBJECTS) beam2016implementation.o system.o -o $(PROGNAME_MWA) $(INCLUDES) $(LIBS) $(COMMON_LIBS) -L $(HDF5_LIB_DIR) -lhdf5_cpp -lhdf5 -lboost_system -lboost_filesystem -D_MWA_2016_BEAM_MODEL_

beam2016test : beam2016test.cpp beam2016implementation.cpp 
	@echo
	@echo "Compiling beam2016test"
	g++  $(OPT) beam2016implementation.cpp -c -I $(HDF5_INCLUDE_DIR)
	g++  $(OPT) system.cpp -c -std=gnu++11
	g++ $(OPT) -I/usr/include/hdf5/serial beam2016test.cpp beam2016implementation.o system.o -o beam2016test $(INCLUDES) $(LIBS) $(COMMON_LIBS) -L $(HDF5_LIB_DIR) -lhdf5_cpp -lhdf5 -lboost_system -lboost_filesystem -D_MWA_2016_BEAM_MODEL_ -lz -lpthread -lsz -ldl 
   

$(OBJECTS) : %.o :  %.cpp
	time -p $(CPP_COMP) $(INCLUDES) $(CMN_INCLUDES) $(PREPROC) $(CCFLAGS) $(OPT) -c $<

$(C_OBJECTS) : %.o :  %.c
	time -p $(C_COMP) $(CCFLAGS) -c $<
        
        
clean : 
	rm -f *.o $(PROGNAME) $(PROGNAME_MWA)
		