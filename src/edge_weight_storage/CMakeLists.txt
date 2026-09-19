add_library(edge_weight_storage
        wgt_storage_interface.c wgt_storage_interface.h
        cmap.c cmap.h
        qisq2_map.c qisq2_map.h
        fast_hash.h fast_hash.c
        flt.h
        MurmurHash3.h MurmurHash3.c
        util.h)

target_include_directories(edge_weight_storage
    PUBLIC
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>"
        "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/..>"
        "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>"
)

# qisq2_map.h includes gmp.h and is public, so this propagates to anything
# linking the storage library.
target_link_libraries(edge_weight_storage PUBLIC GMP::GMP)

# Match the main library: PIC when Q-Sylvan is shared or explicitly asks for it.
if(BUILD_SHARED_LIBS OR SYLVAN_ENABLE_PIC)
    set_target_properties(edge_weight_storage PROPERTIES POSITION_INDEPENDENT_CODE ON)
endif()

option(SYLVAN_QUADMATH "128 bit arithmetic" OFF)
if(SYLVAN_QUADMATH)
    target_link_libraries(edge_weight_storage PUBLIC "-lquadmath")
endif()

include(GNUInstallDirs)
install(TARGETS edge_weight_storage
    EXPORT  qsylvan-targets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

if(BUILD_TESTING)
    add_subdirectory(test)
endif()
