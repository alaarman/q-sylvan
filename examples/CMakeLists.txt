if(APPLE)
    # add argp library for macOS
    find_package(Argp REQUIRED)
    include_directories(${ARGP_INCLUDE_PATH})
endif()

macro(set_compilation_settings NAME)
    target_include_directories(${NAME} PRIVATE ${CMAKE_CURRENT_LIST_DIR})
    target_link_libraries(${NAME} PRIVATE qsylvan::sylvan ${ARGP_LIBRARIES})
    target_compile_features(${NAME} PRIVATE c_std_11 cxx_std_11)
    if(CMAKE_C_COMPILER_ID MATCHES "Clang|GNU")
        target_compile_options(${NAME} PRIVATE -Wall -Wextra -Wno-unused -Wno-deprecated)
    endif()
endmacro(set_compilation_settings)

macro(add_example NAME SOURCE)
    add_executable(${NAME} ${SOURCE})
    set_compilation_settings(${NAME})
endmacro(add_example)

add_example(bell_state bell_state.c)
add_example(vqc vqc.c)

set(ALGORITHM_EXAMPLES
    grover_cnf.c
    grover.c
    random_circuit.c
    shor.c
    supremacy.c
)

add_example(hybrid_demo hybrid_demo.c)
add_example(alg_run alg_run.c)
target_sources(alg_run PRIVATE ${ALGORITHM_EXAMPLES})

add_example(test_algs test_algs.c)
target_sources(test_algs PRIVATE ${ALGORITHM_EXAMPLES})

add_executable(circuit_equivalence circuit_equivalence.c)
set_compilation_settings(circuit_equivalence)
target_link_libraries(circuit_equivalence PRIVATE qsylvan_qasm_parser)
