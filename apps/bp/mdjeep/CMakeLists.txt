set(YEWPAR_BUILD_BP_APPS_MDJEEP "ON" CACHE BOOL "Build MDjeep")

if(YEWPAR_BUILD_BP_APPS_DMDGP)

add_hpx_executable(mdjeep
  SOURCES main.cpp bp.cpp
  DEPENDENCIES YewPar
)
endif(YEWPAR_BUILD_BP_APPS_MDJEEP)