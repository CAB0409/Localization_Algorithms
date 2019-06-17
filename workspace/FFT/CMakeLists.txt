
# ------ Includes
include_directories("${PROJECT_SOURCE_DIR}/include")
#FFTW
find_package(FFTW REQUIRED)

add_definitions(-std=c++1y)
add_subdirectory(examples)

include_directories(${FFTW_INCLUDES})
include_directories("${ARMADILLO_INCLUDE_DIRS}")

add_executable(fft fft_test.cpp)
target_link_libraries(fft
                      DARC_waveformlib
                      ${ARMADILLO_LIBRARIES}
                      ${FFTW_LIBRARIES})
