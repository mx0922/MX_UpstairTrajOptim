set(LIB_NAME MX)
set(${LIB_NAME}_DIR ${ROOT_PATH}/test_ipopt_cmake)
message("-- [${LIB_NAME}]: Hello! I'm in ${${LIB_NAME}_DIR}")

include_directories(${${LIB_NAME}_DIR}/include)

aux_source_directory(${${LIB_NAME}_DIR}/src ANKLE_TRAJ_NLP_FILES)