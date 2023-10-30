set(YEWPAR_BUILD_BP_APPS_DMDGP "ON" CACHE BOOL "Build DMDGP")

if(YEWPAR_BUILD_BP_APPS_DMDGP)

add_hpx_executable(dmdgp
  SOURCES main.cpp readfile.c bp.c distance.c matrices.c objfun.c printfile.c pruningtest.c splitime.c utils.c vertex.c bp.hpp
  DEPENDENCIES YewPar
)
endif(YEWPAR_BUILD_ENUMERATION_APPS_UTS)