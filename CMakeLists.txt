# cmake version to be used
cmake_minimum_required( VERSION 3.0 )
set(CMAKE_FIND_DEBUG_MODE 1)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# project name
project(watsonfit)

# target
add_library(watsonfit SHARED src/bonndit/utilc/watsonfit.cpp)

# external libs
set(Python_ADDITIONAL_VERSIONS 3.9)
find_package(PythonLibs 3 REQUIRED)
find_package(Python3 COMPONENTS Interpreter NumPy REQUIRED)
include_directories(${PYTHON_INCLUDE_DIRS})

find_package(Eigen3 REQUIRED)

find_package(Ceres REQUIRED)
find_package(OpenMP REQUIRED)

target_compile_options(watsonfit PRIVATE -DCYTHON_PEP489_MULTI_PHASE_INIT=0)

include_directories(${PYTHON_INCLUDE_DIRS})
target_include_directories(watsonfit
  PRIVATE
    ${CERES_INCLUDE_DIRS}
)
target_link_libraries(watsonfit 
  PRIVATE 
    ${PYTHON_LIBRARIES}
)
target_link_libraries(watsonfit PRIVATE ${PYTHON_LIBRARIES} Python3::NumPy)
target_link_libraries(watsonfit
  PRIVATE
    ${CERES_LIBRARIES}
)
target_link_libraries(watsonfit PRIVATE OpenMP::OpenMP_CXX)
target_link_libraries(watsonfit PRIVATE fftw3)
target_link_libraries(watsonfit PRIVATE cerf)