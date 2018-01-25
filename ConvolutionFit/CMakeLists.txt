INCLUDE_DIRECTORIES(.)
INCLUDE_DIRECTORIES(/usr/include/thrust)

FIND_PACKAGE(CUDA REQUIRED)

set(ROOT_DICT_INPUT_HEADERS RooMyPdf.h LinkDef.h)
GEN_ROOT_DICT_SOURCE(RooMyPdfDict.cxx)
cuda_add_executable(OriginMethod OriginMethod.cpp RooMyPdf.cxx fastgl.cpp cuda_calculate.cu ${ROOT_DICT_OUTPUT_SOURCES})
target_link_libraries(OriginMethod RooFit RooFitCore)

