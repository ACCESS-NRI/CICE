cmake_minimum_required(VERSION 3.18)

#[==============================================================================[
#                         Basic project definition                              #
#]==============================================================================]

project(CICE 
        DESCRIPTION "CICE Sea Ice Model"
        HOMEPAGE_URL https://github.com/ciCE-Consortium/cice
        LANGUAGES C Fortran)

#[==============================================================================[
#                                   Options                                     #
#]==============================================================================]

# Build options
option(ACCESS3_LIB_INSTALL   "Instal ACCESS3 libraries" OFF)
option(CICE_IO               "CICE IO Method" "Binary")
option(CESMCOUPLED           "CESMCOUPLED Build CPP" OFF)
# openmp ?

message(STATUS "Build options")
message(STATUS "  - ACCESS3_LIB_INSTALL   ${ACCESS3_LIB_INSTALL}")
message(STATUS "  - CICE_IO   ${CICE_IO}")
message(STATUS "  - CESMCOUPLED   ${CESMCOUPLED}")

#[==============================================================================[
#                           Project configuration                               #
#]==============================================================================]

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR})
include(CMakePackageConfigHelpers)
include(GNUInstallDirs)
include(FortranLib)

# Common compiler flags and definitions
if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fbacktrace -fconvert=big-endian -ffree-line-length-none -ffixed-line-length-none")
  if(${CMAKE_Fortran_COMPILER_VERSION} VERSION_GREATER_EQUAL 10)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fallow-argument-mismatch")
  endif()
  set(CMAKE_Fortran_FLAGS_RELEASE "-O")
  set(CMAKE_Fortran_FLAGS_DEBUG "-g -Wall -Og -ffpe-trap=zero,overflow -fcheck=bounds")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
  set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -qno-opt-dynamic-align  -convert big_endian -assume byterecl -ftz -traceback -assume realloc_lhs -fp-model source")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O2 -debug minimal")
  set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -g -check uninit -check bounds -check pointers -fpe0 -check noarg_temp_created")
else()
  message(WARNING "Fortran compiler with ID ${CMAKE_Fortran_COMPILER_ID} will be used with CMake default options")
endif()

if(CMAKE_C_COMPILER_ID MATCHES "GNU")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=gnu99")
  set(CMAKE_C_FLAGS_RELEASE "-O")
  set(CMAKE_C_FLAGS_DEBUG "-g -Wall -Og -fbacktrace -ffpe-trap=invalid,zero,overflow -fcheck=bounds")
elseif(CMAKE_C_COMPILER_ID MATCHES "Intel")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -traceback -qno-opt-dynamic-align -fp-model precise -std=gnu99")
  set(CMAKE_C_FLAGS_RELEASE "-O2 -debug minimal")
  set(CMAKE_C_FLAGS_DEBUG "-O0 -g")
else()
  message(WARNING "C compiler with ID ${CMAKE_C_COMPILER_ID} will be used with CMake default options")
endif()

if (CESMCOUPLED)
  add_compile_definitions(
    CESMCOUPLED
  )
endif()

## Fortran modules path; currently this is simply set to be the include dir
set(CMAKE_INSTALL_MODULEDIR ${CMAKE_INSTALL_INCLUDEDIR}
  CACHE STRING
  "Fortran module installation path (Not a cmake native variable)"
)

#[==============================================================================[
#                              External packages                                #
#]==============================================================================]

if(ACCESS3_LIB_INSTALL)
  find_package(Access3Share) #should this have a REQUIRED?
  find_package(ESMF 8.3.0 MODULE REQUIRED)
else()
  find_package(MPI REQUIRED)
endif()

if(CICE_IO MATCHES "^(NetCDF|PIO)$")
  find_package(NetCDF 4.7.3 REQUIRED Fortran)
endif()
if(CICE_IO MATCHES "PIO")
  find_package(PIO 2.5.3 REQUIRED COMPONENTS C Fortran)
endif()

#[==============================================================================[
#                             Main definitions                                  #
#]==============================================================================]

### Targets

## CICE library
add_fortran_library(cicelib mod STATIC)

if(CICE_IO MATCHES "^(NetCDF|PIO)$")
  # which is this not just an "add_compile_definitions?"
  target_compile_definitions(cicelib PRIVATE FORTRANUNDERSCORE ncdf)
  target_compile_definitions(cicelib PRIVATE USE_NETCDF)
endif()
if(ACCESS3_LIB_INSTALL)
  target_link_libraries(cicelib
    PUBLIC esmf
    PRIVATE Access3Share::nuopc_cap_share Access3Share::share Access3Share::timing Access3Share::cdeps-common
  )
endif()
if(CICE_IO MATCHES "^(NetCDF|PIO)$")
  target_link_libraries(cicelib PRIVATE NetCDF::NetCDF_Fortran)
  if(CICE_IO MATCHES "PIO")
    target_link_libraries(cicelib PRIVATE PIO::PIO_Fortran)
  endif()
endif()

if(OpenMP_Fortran_FOUND)
  target_link_libraries(cicelib PRIVATE OpenMP::OpenMP_Fortran)
endif()

set(CICE_CORE "${CMAKE_SOURCE_DIR}/../cicecore")
set(ICEPACK "${CMAKE_SOURCE_DIR}/../icepack")

target_sources(cicelib PRIVATE
  # Shared List:
  ${CICE_CORE}/shared/ice_arrays_column.F90
  ${CICE_CORE}/shared/ice_calendar.F90
  ${CICE_CORE}/shared/ice_constants.F90
  ${CICE_CORE}/shared/ice_domain_size.F90
  ${CICE_CORE}/shared/ice_fileunits.F90
  ${CICE_CORE}/shared/ice_init_column.F90
  ${CICE_CORE}/shared/ice_kinds_mod.F90
  ${CICE_CORE}/shared/ice_restart_column.F90
  ${CICE_CORE}/shared/ice_restart_shared.F90
  ${CICE_CORE}/shared/ice_spacecurve.F90
  ${CICE_CORE}/shared/ice_distribution.F90

  # Analysis
  ${CICE_CORE}/cicedyn/analysis/ice_diagnostics.F90
  ${CICE_CORE}/cicedyn/analysis/ice_diagnostics_bgc.F90
  ${CICE_CORE}/cicedyn/analysis/ice_history.F90
  ${CICE_CORE}/cicedyn/analysis/ice_history_bgc.F90
  ${CICE_CORE}/cicedyn/analysis/ice_history_drag.F90
  ${CICE_CORE}/cicedyn/analysis/ice_history_fsd.F90
  ${CICE_CORE}/cicedyn/analysis/ice_history_mechred.F90
  ${CICE_CORE}/cicedyn/analysis/ice_history_pond.F90
  ${CICE_CORE}/cicedyn/analysis/ice_history_shared.F90
  ${CICE_CORE}/cicedyn/analysis/ice_history_snow.F90

  # Dynamics
  ${CICE_CORE}/cicedyn/dynamics/ice_dyn_core1d.F90
  ${CICE_CORE}/cicedyn/dynamics/ice_dyn_eap.F90
  ${CICE_CORE}/cicedyn/dynamics/ice_dyn_evp.F90
  ${CICE_CORE}/cicedyn/dynamics/ice_dyn_evp1d.F90
  ${CICE_CORE}/cicedyn/dynamics/ice_dyn_shared.F90
  ${CICE_CORE}/cicedyn/dynamics/ice_dyn_vp.F90
  ${CICE_CORE}/cicedyn/dynamics/ice_transport_driver.F90
  ${CICE_CORE}/cicedyn/dynamics/ice_transport_remap.F90

  # General
  ${CICE_CORE}/cicedyn/general/ice_init.F90
  ${CICE_CORE}/cicedyn/general/ice_flux.F90
  ${CICE_CORE}/cicedyn/general/ice_flux_bgc.F90
  ${CICE_CORE}/cicedyn/general/ice_forcing.F90
  ${CICE_CORE}/cicedyn/general/ice_forcing_bgc.F90
  ${CICE_CORE}/cicedyn/general/ice_state.F90
  ${CICE_CORE}/cicedyn/general/ice_step_mod.F90

  # Infrastructure
  ${CICE_CORE}/cicedyn/infrastructure/ice_blocks.F90
  ${CICE_CORE}/cicedyn/infrastructure/ice_grid.F90
  ${CICE_CORE}/cicedyn/infrastructure/ice_memusage.F90
  ${CICE_CORE}/cicedyn/infrastructure/ice_memusage_gptl.c
  ${CICE_CORE}/cicedyn/infrastructure/ice_read_write.F90
  ${CICE_CORE}/cicedyn/infrastructure/ice_restart_driver.F90
  ${CICE_CORE}/cicedyn/infrastructure/ice_restoring.F90
  ${CICE_CORE}/cicedyn/infrastructure/ice_domain.F90

  # Icepack
  ${ICEPACK}/columnphysics/icepack_aerosol.F90
  ${ICEPACK}/columnphysics/icepack_age.F90
  ${ICEPACK}/columnphysics/icepack_algae.F90
  ${ICEPACK}/columnphysics/icepack_atmo.F90
  ${ICEPACK}/columnphysics/icepack_brine.F90
  ${ICEPACK}/columnphysics/icepack_firstyear.F90
  ${ICEPACK}/columnphysics/icepack_flux.F90
  ${ICEPACK}/columnphysics/icepack_fsd.F90
  ${ICEPACK}/columnphysics/icepack_intfc.F90
  ${ICEPACK}/columnphysics/icepack_isotope.F90
  ${ICEPACK}/columnphysics/icepack_itd.F90
  ${ICEPACK}/columnphysics/icepack_kinds.F90
  ${ICEPACK}/columnphysics/icepack_mechred.F90
  ${ICEPACK}/columnphysics/icepack_meltpond_lvl.F90
  ${ICEPACK}/columnphysics/icepack_meltpond_topo.F90
  ${ICEPACK}/columnphysics/icepack_mushy_physics.F90
  ${ICEPACK}/columnphysics/icepack_ocean.F90
  ${ICEPACK}/columnphysics/icepack_orbital.F90
  ${ICEPACK}/columnphysics/icepack_parameters.F90
  ${ICEPACK}/columnphysics/icepack_shortwave_data.F90
  ${ICEPACK}/columnphysics/icepack_shortwave.F90
  ${ICEPACK}/columnphysics/icepack_snow.F90
  ${ICEPACK}/columnphysics/icepack_therm_bl99.F90
  ${ICEPACK}/columnphysics/icepack_therm_itd.F90
  ${ICEPACK}/columnphysics/icepack_therm_mushy.F90
  ${ICEPACK}/columnphysics/icepack_therm_shared.F90
  ${ICEPACK}/columnphysics/icepack_therm_vertical.F90
  ${ICEPACK}/columnphysics/icepack_tracers.F90
  ${ICEPACK}/columnphysics/icepack_warnings.F90
  ${ICEPACK}/columnphysics/icepack_wavefracspec.F90
  ${ICEPACK}/columnphysics/icepack_zbgc.F90
  ${ICEPACK}/columnphysics/icepack_zbgc_shared.F90

  # Shared C
  ${CICE_CORE}/cicedyn/infrastructure/ice_shr_reprosum86.c

  # MPI
  ${CICE_CORE}/cicedyn/infrastructure/comm/mpi/ice_boundary.F90
  ${CICE_CORE}/cicedyn/infrastructure/comm/mpi/ice_broadcast.F90
  ${CICE_CORE}/cicedyn/infrastructure/comm/mpi/ice_communicate.F90
  ${CICE_CORE}/cicedyn/infrastructure/comm/mpi/ice_exit.F90
  ${CICE_CORE}/cicedyn/infrastructure/comm/mpi/ice_gather_scatter.F90
  ${CICE_CORE}/cicedyn/infrastructure/comm/mpi/ice_global_reductions.F90
  ${CICE_CORE}/cicedyn/infrastructure/comm/mpi/ice_reprosum.F90
  ${CICE_CORE}/cicedyn/infrastructure/comm/mpi/ice_timers.F90
)

if(ACCESS3_LIB_INSTALL)
  target_sources(cicelib PRIVATE
    # NUOPC CMEPS driver
    ${CICE_CORE}/drivers/nuopc/cmeps/CICE_FinalMod.F90
    ${CICE_CORE}/drivers/nuopc/cmeps/CICE_InitMod.F90
    ${CICE_CORE}/drivers/nuopc/cmeps/CICE_RunMod.F90
    ${CICE_CORE}/drivers/nuopc/cmeps/cice_wrapper_mod.F90
    ${CICE_CORE}/drivers/nuopc/cmeps/ice_comp_nuopc.F90
    ${CICE_CORE}/drivers/nuopc/cmeps/ice_import_export.F90
    ${CICE_CORE}/drivers/nuopc/cmeps/ice_mesh_mod.F90
    ${CICE_CORE}/drivers/nuopc/cmeps/ice_prescribed_mod.F90
    ${CICE_CORE}/drivers/nuopc/cmeps/ice_scam.F90
    ${CICE_CORE}/drivers/nuopc/cmeps/ice_shr_methods.F90
  )
else()
  target_sources(cicelib PRIVATE
    # CICE standalone
    ${CICE_CORE}/drivers/standalone/cice/CICE_FinalMod.F90
    ${CICE_CORE}/drivers/standalone/cice/CICE_InitMod.F90
    ${CICE_CORE}/drivers/standalone/cice/CICE_RunMod.F90
  )
endif()

# Select IO source files based on CICE_IO
if(CICE_IO MATCHES "NetCDF")
  target_sources(cicelib PRIVATE
    ${CICE_CORE}/cicedyn/infrastructure/io/io_netcdf/ice_history_write.F90
    ${CICE_CORE}/cicedyn/infrastructure/io/io_netcdf/ice_restart.F90
  )
elseif(CICE_IO MATCHES "PIO")
  target_sources(cicelib PRIVATE
    ${CICE_CORE}/cicedyn/infrastructure/io/io_pio2/ice_history_write.F90
    ${CICE_CORE}/cicedyn/infrastructure/io/io_pio2/ice_pio.F90
    ${CICE_CORE}/cicedyn/infrastructure/io/io_pio2/ice_restart.F90
  )
elseif(CICE_IO MATCHES "Binary")
  target_sources(cicelib PRIVATE
    ${CICE_CORE}/cicedyn/infrastructure/io/io_binary/ice_history_write.F90
    ${CICE_CORE}/cicedyn/infrastructure/io/io_binary/ice_restart.F90
  )
endif()

#[==============================================================================[
#                             Install or Export                                 #
#]==============================================================================]

if(ACCESS3_LIB_INSTALL)
  ## Library
  set_target_properties(cicelib PROPERTIES
    OUTPUT_NAME access-cice-cmeps
    EXPORT_NAME cice-cmeps
  )
  install(TARGETS cicelib
    EXPORT AccessCiceCmepsTargets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} 
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR} 
  )
  # Fortran module files are a special case, as currently there is no standard
  # way of handling them in CMake
  target_include_directories(cicelib PUBLIC "$<INSTALL_INTERFACE:${CMAKE_INSTALL_MODULEDIR}/access-cice>")
  get_target_property(cice_moddir cicelib Fortran_MODULE_DIRECTORY)
  install(FILES ${cice_moddir}/ice_comp_nuopc.mod
    DESTINATION ${CMAKE_INSTALL_MODULEDIR}/access-cice
  )

  # Make sure the dependencies get exported too
  configure_package_config_file(
    AccessCiceCmeps.cmake.in
    "${CMAKE_CURRENT_BINARY_DIR}/AccessCiceCmepsConfig.cmake"
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/AccessCiceCmeps
  )
  install(FILES ${CMAKE_SOURCE_DIR}/FindNetCDF.cmake ${CMAKE_SOURCE_DIR}/FindPIO.cmake ${CMAKE_CURRENT_BINARY_DIR}/AccessCiceCmepsConfig.cmake
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/AccessCiceCmeps
    COMPONENT AccessCiceCmeps
  )
  # install(EXPORT AccessCiceCmepsTargets
  #   FILE AccessCiceCmepsTargets.cmake
  #   NAMESPACE AccessCiceCmeps::
  #   DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/AccessCiceCmeps
  # )
else()
  set_target_properties(cicelib PROPERTIES
      OUTPUT_NAME cicelib-standalone
      EXPORT_NAME cicelib
  )
  add_executable(CICE ${CICE_CORE}/drivers/standalone/cice/CICE.F90)
  target_link_libraries(CICE PRIVATE cicelib)

  set_target_properties(CICE PROPERTIES
    LINKER_LANGUAGE Fortran
    OUTPUT_NAME cice
  )

  install(TARGETS CICE
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  )
endif()
