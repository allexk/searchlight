# generate protobuf files
PROTOBUF_GENERATE_CPP(MSG_SRC MSG_HDRS searchlight/searchlight_messages.proto)
include_directories(${CMAKE_CURRENT_BINARY_DIR})

# searchlight SciDb UDO
add_library(searchlight_udo SHARED searchlight/adapter.cpp
    searchlight/array_access.cpp searchlight/array_desc.cpp
    searchlight/ops/LogicalSearchlight.cpp searchlight/ops/PhysicalSearchlight.cpp
    searchlight/sampler.cpp searchlight/searchlight_task.cpp
    searchlight/searchlight.cpp searchlight/validator.cpp
    searchlight/searchlight_collector.cpp searchlight/ops/plugin.cpp
    searchlight/default_search.cpp searchlight/searchlight_messenger.cpp ${MSG_SRC}
    searchlight/ops/LogicalDepart.cpp searchlight/ops/PhysicalDepart.cpp
    searchlight/ops/DepartArray.cpp)

# We have to link the libraries since SciDb will not load them during the
# dynamic operator loading. We don't need to do anything special for SciDb
# libs, though. We also link UDFs and tasks with the UDO for proper
# symbol resolution.
target_link_libraries(searchlight_udo ${ORTOOLS_LIBRARIES})
target_link_libraries(searchlight_udo ${GFLAGS_LIBRARIES})
target_link_libraries(searchlight_udo ${COIN_LIBRARIES})

# searchlight UDFs
add_library(searchlight_udfs SHARED searchlight/searchlight_udfs.cpp)
target_link_libraries(searchlight_udfs searchlight_udo)

# example tasks
add_library(searchlight_sw SHARED tasks/sw.cpp)
target_link_libraries(searchlight_sw searchlight_udo)

#install rules for the Searchlight libs
install(TARGETS searchlight_udo searchlight_udfs searchlight_sw
    DESTINATION lib/scidb/plugins)
