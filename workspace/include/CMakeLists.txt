set(CMAKE_CXX_FILES "-Wall -03 ${CMAKE_CXX_FLAGS}")

add_library(spdlog INTERFACE)

install(FILES
        spdlog
        zhelpers.hpp
        zhelpers.h
        zmqThreadHelper.hpp
        DESTINATION include
        )
