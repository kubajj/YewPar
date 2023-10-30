set(YEWPAR_BUILD_BP_APPS_DMDGP "ON" CACHE BOOL "Build DMDGP")

if(YEWPAR_BUILD_BP_APPS_DMDGP)

add_hpx_executable(dmdgp
  SOURCES main.cpp readfile.c distance.c utils.c bp.hpp
  DEPENDENCIES YewPar
)
endif(YEWPAR_BUILD_BP_APPS_DMDGP)