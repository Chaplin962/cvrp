<<<<<<< HEAD
cmake_minimum_required(VERSION 3.15)
project(HGS_CVRP)
set(CMAKE_CXX_STANDARD 17)

set(
        src_files
        Program/Genetic.cu
        Program/Individual.cpp
        Program/LocalSearch.cu
        Program/Params.cpp
        Program/Population.cpp
        Program/Split.cpp
        Program/InstanceCVRPLIB.cpp
        Program/AlgorithmParameters.cpp
        Program/C_Interface.cpp)

if (MSVC)
	set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
endif (MSVC)

include_directories(Program)

# Build Executable
add_executable(bin
        Program/main.cu
        ${src_files})

set_target_properties(bin PROPERTIES OUTPUT_NAME hgs)

# Test Executable
include(CTest)
add_test(NAME    bin_test_X-n101-k25
         COMMAND ${CMAKE_COMMAND}   -DINSTANCE=X-n101-k25
                                    -DCOST=27591
                                    -DROUND=1
                                    -P ${CMAKE_CURRENT_SOURCE_DIR}/Test/TestExecutable.cmake)
add_test(NAME    bin_test_X-n106-k14
         COMMAND ${CMAKE_COMMAND}   -DINSTANCE=X-n110-k13
                                    -DCOST=14971
                                    -DROUND=1
                                    -P ${CMAKE_CURRENT_SOURCE_DIR}/Test/TestExecutable.cmake)

# Test Executable: Instances with Duration, without Rounding
add_test(NAME    bin_test_CMT6
        COMMAND ${CMAKE_COMMAND}    -DINSTANCE=CMT6
                                    -DCOST=555.43
                                    -DROUND=0
                                    -P ${CMAKE_CURRENT_SOURCE_DIR}/Test/TestExecutable.cmake)
add_test(NAME    bin_test_CMT7
        COMMAND ${CMAKE_COMMAND}    -DINSTANCE=CMT7
                                    -DCOST=909.675
                                    -DROUND=0
                                    -P ${CMAKE_CURRENT_SOURCE_DIR}/Test/TestExecutable.cmake)

# Build Library
add_library(lib SHARED ${src_files})
set_target_properties(lib PROPERTIES OUTPUT_NAME hgscvrp)

# Test Library
add_subdirectory(Test/Test-c/)
add_test(NAME       lib_test_c
         COMMAND    lib_test_c)

# Install
install(TARGETS lib
        DESTINATION lib)
install(TARGETS bin
        DESTINATION bin)
install(FILES Program/AlgorithmParameters.h Program/C_Interface.h
=======
cmake_minimum_required(VERSION 3.8 FATAL_ERROR)
project(HGS_CVRP LANGUAGES CXX CUDA)
set(CMAKE_CXX_STANDARD 17)

#[edit] for cuda
cmake_minimum_required(VERSION 3.9 FATAL_ERROR)

set(
        src_files
        Program/Genetic.cu
        Program/Individual.cpp
        Program/LocalSearch.cu
        Program/Params.cpp
        Program/Population.cpp
        Program/Split.cpp
        Program/InstanceCVRPLIB.cpp
        Program/AlgorithmParameters.cpp
        Program/C_Interface.cpp)

if (MSVC)
	set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
endif (MSVC)

include_directories(Program)

# Build Executable
add_executable(bin
        Program/main.cu
        ${src_files})

set_target_properties(bin PROPERTIES OUTPUT_NAME hgs)

# Test Executable
include(CTest)
add_test(NAME    bin_test_X-n101-k25
         COMMAND ${CMAKE_COMMAND}   -DINSTANCE=X-n101-k25
                                    -DCOST=27591
                                    -DROUND=1
                                    -P ${CMAKE_CURRENT_SOURCE_DIR}/Test/TestExecutable.cmake)
add_test(NAME    bin_test_X-n106-k14
         COMMAND ${CMAKE_COMMAND}   -DINSTANCE=X-n110-k13
                                    -DCOST=14971
                                    -DROUND=1
                                    -P ${CMAKE_CURRENT_SOURCE_DIR}/Test/TestExecutable.cmake)

# Test Executable: Instances with Duration, without Rounding
add_test(NAME    bin_test_CMT6
        COMMAND ${CMAKE_COMMAND}    -DINSTANCE=CMT6
                                    -DCOST=555.43
                                    -DROUND=0
                                    -P ${CMAKE_CURRENT_SOURCE_DIR}/Test/TestExecutable.cmake)
add_test(NAME    bin_test_CMT7
        COMMAND ${CMAKE_COMMAND}    -DINSTANCE=CMT7
                                    -DCOST=909.675
                                    -DROUND=0
                                    -P ${CMAKE_CURRENT_SOURCE_DIR}/Test/TestExecutable.cmake)

# Build Library
add_library(lib SHARED ${src_files})
set_target_properties(lib PROPERTIES OUTPUT_NAME hgscvrp)

# Test Library
add_subdirectory(Test/Test-c/)
add_test(NAME       lib_test_c
         COMMAND    lib_test_c)

# Install
install(TARGETS lib
        DESTINATION lib)
install(TARGETS bin
        DESTINATION bin)
install(FILES Program/AlgorithmParameters.h Program/C_Interface.h
>>>>>>> 62dfb14936b294b99fa4505e58b7f299baacbf3e
        DESTINATION include)