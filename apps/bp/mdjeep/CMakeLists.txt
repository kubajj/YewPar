set(YEWPAR_BUILD_BP_APPS_MDJEEP "ON" CACHE BOOL "Build MDjeep")

if(YEWPAR_BUILD_BP_APPS_MDJEEP)

add_hpx_executable(mdjeep
  SOURCES main.cpp bp.hpp
  DEPENDENCIES YewPar
)
endif(YEWPAR_BUILD_BP_APPS_MDJEEP)