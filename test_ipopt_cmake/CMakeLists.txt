cmake_minimum_required(VERSION 3.0)
project(MX_TEST_IPOPT)

include_directories(G:/mx_ipopt/Ipopt-3.13/include/coin-or)
link_directories(G:/mx_ipopt/Ipopt-3.13/lib)

include_directories(./include)

aux_source_directory(./src DIR_SRC)

add_executable(test ./test/test.cpp ${DIR_SRC})
target_link_libraries(test ipopt.dll.lib ipoptamplinterface.dll.lib)
