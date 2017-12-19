cmake_minimum_required(VERSION 2.8)
project(RootScript)

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})

find_package(ROOT)
include(cmake/MacroRootDict.cmake)

execute_process(
        COMMAND ${ROOT_CONFIG_EXECUTABLE} --libs --cflags
        OUTPUT_VARIABLE ROOT_FLAGS
        OUTPUT_STRIP_TRAILING_WHITESPACE)

set(CMAKE_CXX_FLAGS "-Wl,--no-as-needed ${CMAKE_CXX_FLAGS} ${ROOT_FLAGS} -std=c++11 -march=native")

include_directories(ArgParser)

add_subdirectory(LatexAnalyse)
add_subdirectory(CountingStationAnalyse)
add_subdirectory(MLData)
#add_subdirectory(MulitFuncFit)
add_subdirectory(LengthDistribution)
add_subdirectory(Bamboo2REST)
add_subdirectory(TimeSplit)
add_subdirectory(RootExam)
add_subdirectory(PrototypeCut)
add_subdirectory(BambooKeras)
add_subdirectory(PrototypeElectron)
add_subdirectory(GasCalculate)
add_subdirectory(ConvolutionFit)

