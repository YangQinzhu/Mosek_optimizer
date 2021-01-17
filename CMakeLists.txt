cmake_minimum_required(VERSION 2.8)

project(linear_optimization)

set(CMAKE_CXX_FLAGS "-std=c++11")

#设置头文件路径
include_directories("/home/yqz/data/homework/mosek/9.2/tools/platform/linux64x86/h")

#设置库文件路径
link_directories("/home/yqz/data/homework/mosek/9.2/tools/platform/linux64x86/bin")

add_executable(mosek_lianxi /home/yqz/data/homework/mosek/9.2/tools/examples/fusion/cxx/test_new/facility_location_yqz.cc)
target_link_libraries(mosek_lianxi libmosek64.so libfusion64.so)
