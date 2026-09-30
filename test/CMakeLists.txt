function(sylvan_test_warnings target)
    target_compile_options(${target} PRIVATE
        $<$<COMPILE_LANG_AND_ID:C,GNU,Clang,AppleClang>:-Wall -Wextra -Werror -Wno-deprecated>
        $<$<COMPILE_LANG_AND_ID:CXX,GNU,Clang,AppleClang>:-Wall -Wextra -Werror -Wno-deprecated>
        $<$<COMPILE_LANG_AND_ID:C,MSVC>:/W4>
        $<$<COMPILE_LANG_AND_ID:CXX,MSVC>:/W4>
    )
endfunction()

# ── Upstream Sylvan tests ────────────────────────────────────────────────────

add_executable(test_basic)
target_sources(test_basic PRIVATE test_basic.c)
target_link_libraries(test_basic PRIVATE qsylvan::sylvan)
target_compile_features(test_basic PRIVATE c_std_11)
sylvan_test_warnings(test_basic)

add_executable(test_zdd)
target_sources(test_zdd PRIVATE test_zdd.c)
target_link_libraries(test_zdd PRIVATE qsylvan::sylvan)
target_compile_features(test_zdd PRIVATE c_std_11)
sylvan_test_warnings(test_zdd)

add_executable(test_table)
target_sources(test_table PRIVATE test_table.c)
target_link_libraries(test_table PRIVATE qsylvan::sylvan)
target_compile_features(test_table PRIVATE c_std_11)
sylvan_test_warnings(test_table)

add_executable(test_cxx)
target_sources(test_cxx PRIVATE test_cxx.cpp)
target_link_libraries(test_cxx PRIVATE qsylvan::sylvan)
target_compile_features(test_cxx PRIVATE cxx_std_11)
sylvan_test_warnings(test_cxx)

add_executable(test_table_multi)
target_sources(test_table_multi PRIVATE test_table_multi.c)
target_link_libraries(test_table_multi PRIVATE qsylvan::sylvan)
target_compile_features(test_table_multi PRIVATE c_std_11)
sylvan_test_warnings(test_table_multi)

add_test(NAME test_table_multi COMMAND test_table_multi)

add_test(NAME test_basic COMMAND test_basic)
add_test(NAME test_zdd   COMMAND test_zdd)
add_test(NAME test_table COMMAND test_table)
add_test(NAME test_cxx   COMMAND test_cxx)

# ── Q-Sylvan MTBDD tests ─────────────────────────────────────────────────────

add_executable(test_mtbdd)
target_sources(test_mtbdd PRIVATE test_mtbdd.c)
target_link_libraries(test_mtbdd PRIVATE qsylvan::sylvan)
target_compile_features(test_mtbdd PRIVATE c_std_11)
sylvan_test_warnings(test_mtbdd)

add_executable(test_simulator_mtbdd)
target_sources(test_simulator_mtbdd PRIVATE test_simulator_mtbdd.c)
target_link_libraries(test_simulator_mtbdd PRIVATE qsylvan::sylvan)
target_compile_features(test_simulator_mtbdd PRIVATE c_std_11)
sylvan_test_warnings(test_simulator_mtbdd)

add_test(NAME test_mtbdd           COMMAND test_mtbdd)
add_test(NAME test_simulator_mtbdd COMMAND test_simulator_mtbdd)

# ── Q-Sylvan LIMDD tests ─────────────────────────────────────────────────────

add_executable(test_limdd_pauli)
target_sources(test_limdd_pauli PRIVATE test_limdd_pauli.c)
target_link_libraries(test_limdd_pauli PRIVATE qsylvan::sylvan)
target_compile_features(test_limdd_pauli PRIVATE c_std_11)
sylvan_test_warnings(test_limdd_pauli)

add_test(NAME test_limdd_pauli COMMAND test_limdd_pauli)

add_executable(test_limdd_lim)
target_sources(test_limdd_lim PRIVATE test_limdd_lim.c)
target_link_libraries(test_limdd_lim PRIVATE qsylvan::sylvan)
target_compile_features(test_limdd_lim PRIVATE c_std_11)
sylvan_test_warnings(test_limdd_lim)

add_test(NAME test_limdd_lim COMMAND test_limdd_lim)

add_executable(test_limdd_node)
target_sources(test_limdd_node PRIVATE test_limdd_node.c)
target_link_libraries(test_limdd_node PRIVATE qsylvan::sylvan)
target_compile_features(test_limdd_node PRIVATE c_std_11)
sylvan_test_warnings(test_limdd_node)

add_test(NAME test_limdd_node COMMAND test_limdd_node)

add_executable(test_limdd_stab)
target_sources(test_limdd_stab PRIVATE test_limdd_stab.c)
target_link_libraries(test_limdd_stab PRIVATE qsylvan::sylvan)
target_compile_features(test_limdd_stab PRIVATE c_std_11)
sylvan_test_warnings(test_limdd_stab)

add_test(NAME test_limdd_stab COMMAND test_limdd_stab)

add_executable(test_limdd_stab_node)
target_sources(test_limdd_stab_node PRIVATE test_limdd_stab_node.c)
target_link_libraries(test_limdd_stab_node PRIVATE qsylvan::sylvan)
target_compile_features(test_limdd_stab_node PRIVATE c_std_11)
sylvan_test_warnings(test_limdd_stab_node)

add_test(NAME test_limdd_stab_node COMMAND test_limdd_stab_node)

add_executable(test_limdd_coset)
target_sources(test_limdd_coset PRIVATE test_limdd_coset.c)
target_link_libraries(test_limdd_coset PRIVATE qsylvan::sylvan)
target_compile_features(test_limdd_coset PRIVATE c_std_11)
sylvan_test_warnings(test_limdd_coset)

add_test(NAME test_limdd_coset COMMAND test_limdd_coset)

add_executable(test_limdd_canon)
target_sources(test_limdd_canon PRIVATE test_limdd_canon.c)
target_link_libraries(test_limdd_canon PRIVATE qsylvan::sylvan)
target_compile_features(test_limdd_canon PRIVATE c_std_11)
sylvan_test_warnings(test_limdd_canon)

add_test(NAME test_limdd_canon COMMAND test_limdd_canon)

add_executable(test_limdd_gc)
target_sources(test_limdd_gc PRIVATE test_limdd_gc.c)
target_link_libraries(test_limdd_gc PRIVATE qsylvan::sylvan)
target_compile_features(test_limdd_gc PRIVATE c_std_11)
sylvan_test_warnings(test_limdd_gc)

add_test(NAME test_limdd_gc COMMAND test_limdd_gc)

add_executable(test_limdd_structure)
target_sources(test_limdd_structure PRIVATE test_limdd_structure.c)
target_link_libraries(test_limdd_structure PRIVATE qsylvan::sylvan)
target_compile_features(test_limdd_structure PRIVATE c_std_11)
sylvan_test_warnings(test_limdd_structure)

add_test(NAME test_limdd_structure COMMAND test_limdd_structure)

add_executable(test_limdd_shape)
target_sources(test_limdd_shape PRIVATE test_limdd_shape.c)
target_link_libraries(test_limdd_shape PRIVATE qsylvan::sylvan)
target_compile_features(test_limdd_shape PRIVATE c_std_11)
sylvan_test_warnings(test_limdd_shape)

add_test(NAME test_limdd_shape COMMAND test_limdd_shape)

add_executable(test_bqd)
target_sources(test_bqd PRIVATE test_bqd.c)
target_link_libraries(test_bqd PRIVATE qsylvan::sylvan)
target_compile_features(test_bqd PRIVATE c_std_11)
sylvan_test_warnings(test_bqd)

add_test(NAME test_bqd COMMAND test_bqd)

# Two cache operations with one id read each other's entries; nothing else
# would notice. Checked on the header itself, so it needs no build.
add_test(NAME test_cache_opids
         COMMAND ${CMAKE_COMMAND} -DHEADER=${PROJECT_SOURCE_DIR}/src/sylvan_int.h
                 -P ${CMAKE_CURRENT_SOURCE_DIR}/check_cache_opids.cmake)

add_executable(test_bqd_ops)
target_sources(test_bqd_ops PRIVATE test_bqd_ops.c)
target_link_libraries(test_bqd_ops PRIVATE qsylvan::sylvan)
target_compile_features(test_bqd_ops PRIVATE c_std_11)
sylvan_test_warnings(test_bqd_ops)

add_test(NAME test_bqd_ops COMMAND test_bqd_ops)

add_executable(test_bqd_gates)
target_sources(test_bqd_gates PRIVATE test_bqd_gates.c)
target_link_libraries(test_bqd_gates PRIVATE qsylvan::sylvan)
target_compile_features(test_bqd_gates PRIVATE c_std_11)
sylvan_test_warnings(test_bqd_gates)

add_test(NAME test_bqd_gates COMMAND test_bqd_gates)

add_executable(test_bqd_skip)
target_sources(test_bqd_skip PRIVATE test_bqd_skip.c)
target_link_libraries(test_bqd_skip PRIVATE qsylvan::sylvan)
target_compile_features(test_bqd_skip PRIVATE c_std_11)
sylvan_test_warnings(test_bqd_skip)

add_test(NAME test_bqd_skip COMMAND test_bqd_skip)

add_executable(test_bqd_ratio)
target_sources(test_bqd_ratio PRIVATE test_bqd_ratio.c)
target_link_libraries(test_bqd_ratio PRIVATE qsylvan::sylvan)
target_compile_features(test_bqd_ratio PRIVATE c_std_11)
sylvan_test_warnings(test_bqd_ratio)

add_test(NAME test_bqd_ratio COMMAND test_bqd_ratio)

add_executable(test_lowlevel_ops)
target_sources(test_lowlevel_ops PRIVATE test_lowlevel_ops.c)
target_link_libraries(test_lowlevel_ops PRIVATE qsylvan::sylvan)
target_compile_features(test_lowlevel_ops PRIVATE c_std_11)
sylvan_test_warnings(test_lowlevel_ops)

add_test(NAME test_lowlevel_ops COMMAND test_lowlevel_ops)

add_executable(test_wgt_table_race)
target_sources(test_wgt_table_race PRIVATE test_wgt_table_race.c)
target_link_libraries(test_wgt_table_race PRIVATE qsylvan::sylvan)
target_compile_features(test_wgt_table_race PRIVATE c_std_11)
sylvan_test_warnings(test_wgt_table_race)

add_test(NAME test_wgt_table_race COMMAND test_wgt_table_race)

add_executable(test_limdd_ops)
target_sources(test_limdd_ops PRIVATE test_limdd_ops.c)
target_link_libraries(test_limdd_ops PRIVATE qsylvan::sylvan)
target_compile_features(test_limdd_ops PRIVATE c_std_11)
sylvan_test_warnings(test_limdd_ops)

add_test(NAME test_limdd_ops COMMAND test_limdd_ops)

# ── Q-Sylvan QMDD tests ──────────────────────────────────────────────────────

foreach(_qmdd_test
        test_qmdd_normalization
        test_qmdd_basics
        test_qmdd_gates
        test_qmdd_circuits
        test_qmdd_matrix)
    add_executable(${_qmdd_test} ${_qmdd_test}.c)
    target_link_libraries(${_qmdd_test} PRIVATE qsylvan::sylvan)
    add_test(NAME ${_qmdd_test} COMMAND ${_qmdd_test})
endforeach()

add_executable(test_qmdd_gc test_qmdd_gc.c ../examples/grover.c)
target_link_libraries(test_qmdd_gc PRIVATE qsylvan::sylvan)
add_test(NAME test_qmdd_gc COMMAND test_qmdd_gc)
