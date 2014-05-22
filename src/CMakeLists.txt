#includes

# searchlight SciDb UDO
add_library(searchlight_udo SHARED searchlight/adapter.cpp
    searchlight/array_access.cpp searchlight/array_desc.cpp
    searchlight/LogicalSearchlight.cpp searchlight/PhysicalSearchlight.cpp
    searchlight/sampler.cpp searchlight/searchlight_task.cpp
    searchlight/searchlight.cpp searchlight/validator.cpp
    searchlight/plugin.cpp)

# We have to link the libraries since SciDb will not load them during the
# dynamic operator loading. We don't need to do anything special for SciDb
# libs, though. Also, UDFs and tasks are run via the UDO, so we do not need
# to link them against or-tools either.
target_link_libraries(searchlight_udo ${ORTOOLS_LIBRARIES})
target_link_libraries(searchlight_udo ${GFLAGS_LIBRARIES})
target_link_libraries(searchlight_udo ${PROTOBUF_LIBRARIES})
target_link_libraries(searchlight_udo ${COIN_LIBRARIES})
set_target_properties(searchlight_udo PROPERTIES OUTPUT_NAME searchlight_udo)

# searchlight UDFs
add_library(searchlight_udfs SHARED searchlight/searchlight_udfs.cpp)
set_target_properties(searchlight_udfs PROPERTIES OUTPUT_NAME searchlight_udfs)

# example tasks
add_library(searchlight_sw SHARED tasks/common.cpp tasks/sw.cpp)
set_target_properties(searchlight_sw PROPERTIES OUTPUT_NAME searchlight_sw)

#install rules
install(TARGETS searchlight_udo searchlight_udfs searchlight_sw
    DESTINATION lib/scidb/plugins)
