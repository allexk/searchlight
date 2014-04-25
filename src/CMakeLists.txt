#includes

# cpdb library
add_library(libsl SHARED ???)
target_link_libraries(libsl ${ORTOOLS_LIBRARIES})
target_link_libraries(libsl ${GFLAGS_LIBRARIES})
target_link_libraries(libsl ${PROTOBUF_LIBRARIES})
target_link_libraries(libsl ${COIN_LIBRARIES})
set_target_properties(libsl PROPERTIES OUTPUT_NAME sl)

# example problems
set(examples sum sim point opt)
foreach (ex ${examples})
    add_executable(sl_${ex} examples/${ex}.cpp)
    target_link_libraries(sl_${ex} libsl ${SYS_LIBRARIES})
endforeach()
