# generate protobuf files
PROTOBUF_GENERATE_CPP(MSG_SRC MSG_HDRS searchlight/searchlight_messages.proto)
include_directories(${CMAKE_CURRENT_BINARY_DIR})

# searchlight SciDb UDO
add_library(searchlight_udo SHARED searchlight/ortools_model.cpp searchlight/base.cpp
    searchlight/ortools_inc.cpp searchlight/adapter.cpp searchlight/synopsis.cpp
    searchlight/array_access.cpp searchlight/array_desc.cpp
    searchlight/ops/LogicalSearchlight.cpp searchlight/ops/PhysicalSearchlight.cpp
    searchlight/sampler.cpp searchlight/searchlight_task.cpp
    searchlight/searchlight.cpp searchlight/validator.cpp
    searchlight/ops/plugin.cpp
    searchlight/default_search.cpp searchlight/searchlight_messenger.cpp ${MSG_SRC}
    searchlight/ops/LogicalDepart.cpp searchlight/ops/PhysicalDepart.cpp
    searchlight/ops/DepartArray.cpp
    searchlight/relax.cpp)

# We have to link the libraries since SciDb will not load them during the
# dynamic operator loading. We don't need to do anything special for SciDb
# libs, though. We also link UDFs and tasks with the UDO for proper
# symbol resolution.
target_link_libraries(searchlight_udo ${ORTOOLS_LIBRARIES})
target_link_libraries(searchlight_udo ${GFLAGS_LIBRARIES})
target_link_libraries(searchlight_udo ${COIN_LIBRARIES})
target_link_libraries(searchlight_udo ${FFTW_LIBRARIES})

# NOTE: we need to link all auxilliary libraries with the UDO library,
# since it is loaded with RTLD_LOCAL by SciDB. That means its dependencies
# on constraint solver and other principal libraries will not work for
# task libraries. Thus, explicit dependencies.

# searchlight UDFs
add_library(searchlight_udfs SHARED searchlight/searchlight_udfs.cpp)
target_link_libraries(searchlight_udfs searchlight_udo)

# example tasks
add_library(searchlight_sw SHARED tasks/sw.cpp)
target_link_libraries(searchlight_sw searchlight_udo)
add_library(searchlight_sdss SHARED tasks/sdss.cpp)
target_link_libraries(searchlight_sdss searchlight_udo)
add_library(searchlight_mimic SHARED tasks/mimic.cpp)
target_link_libraries(searchlight_mimic searchlight_udo)
add_library(searchlight_dist SHARED tasks/dist.cpp)
target_link_libraries(searchlight_dist searchlight_udo)

#install rules for the Searchlight libs
install(TARGETS searchlight_udo searchlight_udfs searchlight_sw
    searchlight_sdss searchlight_mimic searchlight_dist
    DESTINATION ${SL_LIB_DIR})

# utils
add_executable(dft_sampler utils/dft_sampler.cpp)
target_link_libraries(dft_sampler ${FFTW_LIBRARIES})
set_target_properties(dft_sampler PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

add_executable(paa_sampler utils/paa_sampler.cpp)
target_link_libraries(paa_sampler ${SCIDB_LIBRARIES})
target_link_libraries(paa_sampler ${GFLAGS_LIBRARIES})
target_link_libraries(paa_sampler ${CMAKE_THREAD_LIBS_INIT})
set_target_properties(paa_sampler PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
