# ── Platform configuration ──────────────────────────────────────────────────

# Configure platform features before declaring the generated header file set.
option(SYLVAN_USE_MMAP "Let Sylvan use mmap or VirtualAlloc to allocate (virtual) memory" ON)
if(SYLVAN_USE_MMAP AND NOT WIN32)
    include(CheckSymbolExists)
    check_symbol_exists(mmap "sys/mman.h" HAVE_MMAP)
    if(NOT HAVE_MMAP)
        message(WARNING "mmap not found: disabling mmap support")
        set(SYLVAN_USE_MMAP OFF)
    endif()
endif()

option(SYLVAN_STATS "Let Sylvan collect statistics at runtime" OFF)

include(DetectCacheLineSize)
detect_cache_line_size(SYLVAN_CACHE_LINE_SIZE DEFAULT 64)

configure_file(
    sylvan_config.h.in
    "${CMAKE_CURRENT_BINARY_DIR}/sylvan_config.h"
    @ONLY
)

# ── Library target ───────────────────────────────────────────────────────────

add_library(qsylvan)
add_library(qsylvan::sylvan ALIAS qsylvan)

# ── Sources ──────────────────────────────────────────────────────────────────

target_sources(qsylvan
  PRIVATE
    sha2.c
    sylvan_bdd.c
    sylvan_cache.c
    sylvan_common.c
    sylvan_hash.c
    sylvan_ldd.c
    sylvan_mt.c
    sylvan_mtbdd.c
    sylvan_platform.c
    sylvan_refs.c
    sylvan_sl.c
    sylvan_stats.c
    sylvan_table.c
    sylvan_zdd.c
    # Q-Sylvan additions
    qsylvan_gates.c
    qsylvan_gates_mtbdd_mpc.c
    qsylvan_simulator.c
    qsylvan_simulator_mtbdd.c
    sylvan_evbdd.c
    sylvan_edge_weights.c
    sylvan_edge_weights_complex.c
    sylvan_edge_weights_qisq2.c
    qsylvan_limdd_pauli.c
    qsylvan_limdd_lim.c
    qsylvan_limdd_node.c
    qsylvan_limdd_stab.c
    qsylvan_limdd_canon.c
    qsylvan_limdd_gc.c
    qsylvan_limdd_ops.c
    qsylvan_limdd_inspect.c
    qsylvan_bqd.c
    qsylvan_bqd_ops.c
  PUBLIC
    FILE_SET HEADERS
    BASE_DIRS
      "${CMAKE_CURRENT_SOURCE_DIR}"
      "${CMAKE_CURRENT_BINARY_DIR}"
    FILES
      sylvan.h
      sylvan_align.h
      sylvan_bdd.h
      sylvan_cache.h
      sylvan_common.h
      sylvan_hash.h
      sylvan_int.h
      sylvan_ldd.h
      sylvan_ldd_int.h
      sylvan_mt.h
      sylvan_mtbdd.h
      sylvan_mtbdd_int.h
      sylvan_platform.h
      sylvan_stats.h
      sylvan_table.h
      sylvan_zdd.h
      sylvan_zdd_int.h
      # Q-Sylvan additions
      qsylvan.h
      qsylvan_gates.h
      qsylvan_gates_mtbdd_mpc.h
      qsylvan_simulator.h
      qsylvan_simulator_mtbdd.h
      sylvan_evbdd.h
      sylvan_evbdd_int.h
      sylvan_edge_weights.h
      sylvan_edge_weights_complex.h
      sylvan_edge_weights_qisq2.h
      qsylvan_limdd_pauli.h
      qsylvan_limdd_lim.h
      qsylvan_limdd_node.h
      qsylvan_limdd_inspect.h
      qsylvan_bqd.h
      qsylvan_bqd_ops.h
      qsylvan_limdd_stab.h
      qsylvan_limdd_canon.h
      qsylvan_limdd_gc.h
      qsylvan_limdd_ops.h
      "${CMAKE_CURRENT_BINARY_DIR}/sylvan_config.h"
)

# ── C++ bindings ─────────────────────────────────────────────────────────────

target_sources(qsylvan
  PRIVATE sylvan_obj.cpp
  PUBLIC FILE_SET HEADERS FILES sylvan_obj.hpp
)

# ── Edge weight storage (Q-Sylvan) ───────────────────────────────────────────

# GMP is found before the storage library is added, because qisq2 numbers are
# exact algebraic values held in GMP integers and qisq2_map.h -- a PUBLIC
# header of that library -- includes gmp.h.
find_package(GMP REQUIRED)

add_subdirectory(edge_weight_storage)
target_link_libraries(qsylvan PUBLIC edge_weight_storage)

# ── MPFR / MPC support (required by Q-Sylvan) ────────────────────────────────

target_sources(qsylvan
  PRIVATE sylvan_gmp.c
  PUBLIC FILE_SET HEADERS FILES sylvan_gmp.h
)
target_link_libraries(qsylvan PUBLIC GMP::GMP)

find_package(MPFR REQUIRED)
find_package(MPC REQUIRED)
target_sources(qsylvan
  PRIVATE sylvan_mpc.c
  PUBLIC FILE_SET HEADERS FILES sylvan_mpc.h
)
target_include_directories(qsylvan PUBLIC ${MPFR_INCLUDE_DIR} ${MPC_INCLUDE_DIR})
target_link_libraries(qsylvan PUBLIC ${MPFR_LIBRARIES} ${MPC_LIBRARIES})

# ── Compile features ────────────────────────────────────────────────────────

target_compile_features(qsylvan PUBLIC c_std_11 cxx_std_11)

# ── Target properties ───────────────────────────────────────────────────────

# PIC: required for shared libraries, optional for static (e.g. embedding
# Q-Sylvan in a downstream shared library).
if(BUILD_SHARED_LIBS OR SYLVAN_ENABLE_PIC)
    set_target_properties(qsylvan PROPERTIES POSITION_INDEPENDENT_CODE ON)
endif()

# Shared library versioning (ignored for static builds)
if(BUILD_SHARED_LIBS)
    set_target_properties(qsylvan PROPERTIES
        VERSION   ${PROJECT_VERSION}
        SOVERSION ${PROJECT_VERSION_MAJOR}
    )
endif()

# ── LIMDD Pauli width ───────────────────────────────────────────────────────
#
# 64-bit words per component of a Pauli word, so the qubit limit is 64 times
# this. One is the historical width and generates the same code as before it
# was a variable; raising it costs 16 bytes per Pauli table entry per extra
# word and moves those entries off the bucket onto the heap.
set(LIMDD_PAULI_WORDS 1 CACHE STRING "64-bit words per Pauli component (qubit limit is 64x this)")
if(NOT LIMDD_PAULI_WORDS MATCHES "^[1-9][0-9]*$")
    message(FATAL_ERROR "LIMDD_PAULI_WORDS must be a positive integer, got '${LIMDD_PAULI_WORDS}'")
endif()
target_compile_definitions(qsylvan PUBLIC LIMDD_PAULI_WORDS=${LIMDD_PAULI_WORDS})

# ── Compiler flags ──────────────────────────────────────────────────────────

if(CMAKE_C_COMPILER_ID MATCHES "Clang|GNU")
    target_compile_options(qsylvan PRIVATE
        -pipe
        -fno-strict-aliasing
        -Wall -Wextra -Wno-deprecated
        # Errors, not warnings. A limdd_pauli_t is about to become an array
        # pair, and every brace initialiser written as { x, z } then means
        # { {x, z}, {0, 0} } -- the right values in the wrong component, with
        # no diagnostic above a warning. The library has no -Werror (the tests
        # do), so these two would otherwise pass silently.
        -Werror=missing-braces -Werror=missing-field-initializers
        $<$<CONFIG:Debug>:-O0 -Wpedantic -Wvla>
        $<$<BOOL:${SYLVAN_NATIVE_OPT}>:-march=native>
    )
elseif(MSVC)
    target_compile_options(qsylvan
      PRIVATE
        /Zc:__STDC__
        /experimental:c11atomics
        $<$<CONFIG:Debug>:/Od /W4 /Zi>
    )
endif()

# ── Link dependencies ───────────────────────────────────────────────────────

find_package(Threads REQUIRED)
target_link_libraries(qsylvan PUBLIC Threads::Threads)
target_link_libraries(qsylvan PUBLIC lace::lace)
if(NOT APPLE)
    target_link_libraries(qsylvan PUBLIC m)
endif()

# ── Hide tunables when consumed as a subproject ──────────────────────────────

if(NOT PROJECT_IS_TOP_LEVEL)
    mark_as_advanced(
        SYLVAN_ENABLE_PIC
        SYLVAN_NATIVE_OPT
        SYLVAN_USE_MMAP
        SYLVAN_STATS
    )
endif()

# ── Installation ─────────────────────────────────────────────────────────────

include(GNUInstallDirs)

install(TARGETS qsylvan
    EXPORT qsylvan-targets
    ARCHIVE          DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY          DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME          DESTINATION ${CMAKE_INSTALL_BINDIR}
    FILE_SET HEADERS DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(EXPORT qsylvan-targets
    FILE        qsylvan-targets.cmake
    NAMESPACE   qsylvan::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/qsylvan
)

# ── pkg-config ───────────────────────────────────────────────────────────────

set(PKGC_REQUIRES "lace >= 1.7.0 gmp mpfr")
set(PKGC_LIBS_PRIVATE "")
if(NOT APPLE)
    set(PKGC_LIBS_PRIVATE "-lm")
endif()

configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/../qsylvan.pc.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/qsylvan.pc"
    @ONLY
)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/qsylvan.pc
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig
)
