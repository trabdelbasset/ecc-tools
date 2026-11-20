if(DEBUG_IMJ_ANTENNA_REPAIRER)
    message(STATUS "MJ: DEBUG_IMJ_ANTENNA_REPAIRER")
    set(CMAKE_BUILD_TYPE "Debug")
else()
    message(STATUS "MJ: RELEASE_IMJ_ANTENNA_REPAIRER")
    set(CMAKE_BUILD_TYPE "Release")
endif()

add_library(imj_antenna_repairer
    ${IMJ_MODULE}/antenna_repairer/AntennaRepairer.cpp
)

target_link_libraries(imj_antenna_repairer
    PUBLIC
    geometry_db
    imj_data_manager
    imj_module
    imj_toolkit
    idm
)

target_include_directories(imj_antenna_repairer
    PUBLIC
    ${IMJ_MODULE}/antenna_repairer
)
