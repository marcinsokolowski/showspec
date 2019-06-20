# showspec
Program integrating FEKO beam or MWA telescope 2016 beam model with Global Sky Model
Description:

Showspec software enables calculation of antenna temperature using Global Sky Model (http://space.mit.edu/home/angelica/gsm/)
and FEKO generated antenna pattern. The antenna pattern file should be in FEKO output format .fee.
However, a three column text file (with phi, theta and power) can also be handled.
It is also possible to use the new 2016 Murchison Wide-field Array (MWA) beam model. 

Referencing :

If you use GSM model please reference this paper : http://adsabs.harvard.edu/abs/2008MNRAS.388..247D 
and if you use showspec software please reference this paper : http://adsabs.harvard.edu/abs/2015PASA...32....4S


1/ required packages :
   sudo apt-get install libhdf5-dev


2/  Compilation :
   cd showspec_standalone
   mkdir build
   cd build/
   cmake ../
   make VERBOSE=1
   On Ubuntu 16 LTS still required to manually add : -lz -lpthread -lsz -ldl to linking line (shown after make VERBOSE=1)
   make all
   make install
   

3/ Possible ERRORS :
   make VERBOSE=1 

   -  Might require adding in the last command (c++/linking) : -lz -lpthread -lsz -ldl 
       Exacutable files are in bin/ subdirectory, they are : showspec_standalone and showspec_standalone_mwa (enabiling usage of H5 file with 2016 MWA beam model)

   - or on my system I needed to include /usr/lib/libgomp.so.1 to target_link_libraries  in CMakeLists.txt :

     target_link_libraries(showspec_standalone_mwa ${NOVA_LIB} ${FITSIO_LIB} ${HEALPIX_LIBRARIES} ${HDF5_C_LIBRARY} ${HDF5_LIB} ${Boost_SYSTEM_LIBRARY} ${Boost_FILESYSTEM_LIBRARY} /usr/lib/libgomp.so.1)




4/ For more information see : doc/showspec_standalone.pdf (or odt)

5/ The HD5 file with coefficients of spherical harmonics for MWA 2016 beam model can be downloaded by :
    wget http://cerberus.mwa128t.org/mwa_full_embedded_element_pattern.h5

6/ Global Sky Model (GSM) data files can be downloaded from http://xte.mit.edu/angelica/gsm/gsm.tar.gz 
   and other details about the GSM model can be found on this webpage : http://space.mit.edu/home/angelica/gsm/ .
   
   I have also generated gsm files at every 10 MHz step and put them here : 
   The directory gsm/ should contain these .txt files and placed where showspec program is executed (or symbolic link created : ln -s $GSM_FULL_PATH gsm)

   The source code for gsm was also included into this distribution of showspec software ( see gsm/README for details and proper referencing )
