set(YEWPAR_BUILD_BP_APPS_DDGP "ON" CACHE BOOL "Build DDGP")

if(YEWPAR_BUILD_BP_APPS_DDGP)

add_hpx_executable(ddgp
  SOURCES main.cpp bp.hpp readfile.cpp
  DEPENDENCIES YewPar
)
endif(YEWPAR_BUILD_BP_APPS_DDGP)