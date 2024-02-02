set(YEWPAR_BUILD_BP_APPS_MDJEEP "ON" CACHE BOOL "Build MDjeep")

if(YEWPAR_BUILD_BP_APPS_MDJEEP)

add_hpx_executable(mdjeep
  SOURCES main.cpp bp.hpp distance.c matrices.c objfun.c printfile.c pruningtest.c readfile.c spg.c splitime.c utils.c vertex.c bp.c mdjeep.cpp utils.cpp
  DEPENDENCIES YewPar
)
endif(YEWPAR_BUILD_BP_APPS_MDJEEP)