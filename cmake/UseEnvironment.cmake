if(EXISTS $ENV{has_pythia6})
  add_definitions(-DPYTHIA6)
endif()
if(EXISTS $ENV{has_grv})
  add_definitions(-DGRVPDF)
endif()

if($ENV{HOSTNAME} MATCHES "^lxplus[0-9]+.cern.ch")
  set(BASE_DIR "/cvmfs/sft.cern.ch/lcg/external")
  set(GSL_DIR "${BASE_DIR}/GSL/1.14/x86_64-slc5-gcc44-opt")
  set(HEPMC_DIR "${BASE_DIR}/HepMC/2.06.08/x86_64-slc6-gcc48-opt")

  message(STATUS "Compiling on LXPLUS. Do not forget to source the environment variables!")
  #--- searching for GSL
  find_library(GSL_LIB gsl HINTS "${GSL_DIR}/lib")
  find_library(GSL_CBLAS_LIB gslcblas HINTS "${GSL_DIR}/lib")
  #--- searching for HepMC
  find_library(HEPMC_LIB HepMC HINTS "${HEPMC_DIR}/lib")
  find_library(HEPMC_FIO_LIB HepMCfio HINTS "${HEPMC_DIR}/lib")
  find_path(HEPMC_INCLUDE HepMC HINTS "${HEPMC_DIR}/include")
  find_library(LIBCONFIG config++)
else()
  find_library(GSL_LIB gsl)
  find_library(GSL_CBLAS_LIB gslcblas)
  find_library(HEPMC_LIB HepMC)
  find_path(HEPMC_INCLUDE HepMC)
  find_library(LIBCONFIG config++)
endif()

