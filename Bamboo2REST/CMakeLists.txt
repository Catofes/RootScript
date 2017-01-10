add_library(G4DataConvert SHARED G4DataConvert.cpp)
add_executable(main main.cpp)
target_link_libraries(main G4DataConvert)