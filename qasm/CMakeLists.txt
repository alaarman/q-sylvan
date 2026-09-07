# 1. Compile qasm parser as separate library
#   (to allow the parser to be included in other places as a library, without 
#   compiling eval_expr.cpp again)
add_library(qsylvan_qasm_parser)
target_sources(qsylvan_qasm_parser 
    PRIVATE
        qsylvan_qasm_parser.cpp
        parse_math/eval_expr.cpp
    PUBLIC
        qsylvan_qasm_parser.h)

# 1a. Third-party header workaround (this target only).
#   parse_math/exprtk.hpp writes `expression_node<T>::ndb_t::template collect(...)`
#   i.e. the `template` keyword without an explicit template argument list.
#   Clang >= 17 (incl. Apple clang 17) diagnoses this as an ERROR BY DEFAULT under
#   -Wmissing-template-arg-list-after-template-kw. Silencing the diagnostic does not
#   change the generated code; the construct is accepted by every other compiler.
#   We probe the POSITIVE spelling because GCC silently accepts unknown -Wno-* flags,
#   which would make a probe of the negative spelling a false positive.
include(CheckCXXCompilerFlag)
check_cxx_compiler_flag("-Wmissing-template-arg-list-after-template-kw"
                        QSYLVAN_CXX_HAS_W_MISSING_TEMPLATE_ARG_LIST)
if(QSYLVAN_CXX_HAS_W_MISSING_TEMPLATE_ARG_LIST)
    target_compile_options(qsylvan_qasm_parser
        PRIVATE -Wno-missing-template-arg-list-after-template-kw)
endif()

# 1b. argp (not part of libc on macOS); the examples/ targets do the same.
if(APPLE)
    find_package(Argp REQUIRED)
    include_directories(${ARGP_INCLUDE_PATH})
endif()

# 2. QMDD QASM simulator
add_executable(run_qasm_on_qmdd run_qasm_on_qmdd.c)
target_link_libraries(run_qasm_on_qmdd qsylvan qsylvan_qasm_parser ${ARGP_LIBRARIES})

# 3. MTBDD QASM simulator
add_executable(run_qasm_on_mtbdd run_qasm_on_mtbdd.c)
target_link_libraries(run_qasm_on_mtbdd qsylvan qsylvan_qasm_parser ${ARGP_LIBRARIES})
