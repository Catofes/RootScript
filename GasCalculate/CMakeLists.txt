#SET(Garfield_INC_DIR /usr/include/garfieldpp/)
#SET(Garfield_LIB_DIR /usr/lib/)
SET(GARFIELD_HOME /home/herbertqiao/Work/garfield/build)
find_package(Garfield)
if (Garfield_FOUND)
    include_directories(${Garfield_INCLUDE_DIRS})
    add_executable(GasCalculate main.cpp)
    target_link_libraries(GasCalculate ${Garfield_LIBRARIES})

    add_executable(FMT FTM.C)
    target_link_libraries(GasCalculate ${Garfield_LIBRARIES})

    add_executable(SimpleFit fit.cpp)
endif ()