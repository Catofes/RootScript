INCLUDE_DIRECTORIES(Readout_Struct/include)
INCLUDE_DIRECTORIES(.)
set(ROOT_DICT_INPUT_HEADERS ReadoutWave.hh)
GEN_ROOT_DICT_SOURCE(ReadoutWave.cpp)
add_subdirectory(Readout_Struct)
add_library(ReadoutWave SHARED ${ROOT_DICT_OUTPUT_SOURCES})
add_executable(KerasDrift ReadoutWave.cpp KerasDrift.cpp)
target_link_libraries(KerasDrift readout ReadoutWave)

add_executable(AddUUID AddUUID.cpp)

add_executable(KerasReconstruct KerasReconstruct.cpp)
target_link_libraries(KerasReconstruct readout ReadoutWave)

add_executable(KerasPicture KerasPicture.cpp)
target_link_libraries(KerasPicture readout ReadoutWave)

add_executable(MergePicture MergePicture.cpp)
add_executable(DrawWave DrawWave.cpp)
target_link_libraries(DrawWave readout ReadoutWave)

add_executable(3DReconstruct 3DReconstruct.cpp)
target_link_libraries(3DReconstruct readout ReadoutWave)