set(YEWPAR_BUILD_BP_APPS_DMDGP "ON" CACHE BOOL "Build DMDGP")

if(YEWPAR_BUILD_BP_APPS_DMDGP)

add_hpx_executable(dmdgp
  SOURCES main.cpp readfile.cpp vertex.cpp bp.hpp
  DEPENDENCIES YewPar
)
endif(YEWPAR_BUILD_BP_APPS_DMDGP)