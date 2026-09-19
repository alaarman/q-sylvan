cmake_minimum_required(VERSION 3.23...3.31)

project(qsylvan  # based on sylvan 1.11.0
    VERSION 1.2.0
    DESCRIPTION "Q-Sylvan, a parallel decision diagram library for quantum simulation"
    HOMEPAGE_URL "https://github.com/System-Verification-Lab/q-sylvan"
    LANGUAGES C CXX
)

# The upstream Sylvan sources are taken from this release; sylvan_config.h
# reports it, independently of the Q-Sylvan version above.
set(SYLVAN_VERSION_MAJOR 1)
set(SYLVAN_VERSION_MINOR 11)
set(SYLVAN_VERSION_PATCH 0)

# Default to Release for top-level single-config builds
if(PROJECT_IS_TOP_LEVEL AND NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS Debug Release RelWithDebInfo MinSizeRel)
endif()

# Generate compile_commands.json for clangd / IDE integration
if(PROJECT_IS_TOP_LEVEL)
    set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
endif()

# Add CMake modules path
list(PREPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake/")

# ── Options ──────────────────────────────────────────────────────────────────

option(BUILD_SHARED_LIBS          "Build shared libraries"                          OFF)
option(SYLVAN_NATIVE_OPT          "Enable -march=native for compiler optimizations" OFF)
option(SYLVAN_ENABLE_PIC          "Build with -fPIC (for embedding in shared libs)" OFF)
option(SYLVAN_BUILD_EXAMPLES      "Build example tools"         ${PROJECT_IS_TOP_LEVEL})
option(SYLVAN_BUILD_DOCS          "Build documentation"                             OFF)

include(CTest)

# Q-Sylvan's C++ bindings and the QASM parser are always built; the upstream
# SYLVAN_BUILD_CPP switch is kept (fixed to ON) so the merged src/CMakeLists.txt
# and test/CMakeLists.txt can stay close to upstream.
set(SYLVAN_BUILD_CPP ON)

# ── Dependencies ─────────────────────────────────────────────────────────────

# Determine if Lace needs PIC (shared Q-Sylvan or explicit PIC request).
set(_sylvan_needs_pic OFF)
if(BUILD_SHARED_LIBS OR SYLVAN_ENABLE_PIC)
    set(_sylvan_needs_pic ON)
endif()

# First check if Lace is already a target (e.g. a parent project added it),
# then try find_package, then fall back to FetchContent.
if(NOT TARGET lace::lace AND NOT TARGET lace)
    find_package(lace 1.7.0 QUIET)
    if(NOT TARGET lace::lace AND NOT TARGET lace)
        # Normal variables override Lace's options without modifying its cache.
        set(BUILD_SHARED_LIBS OFF)
        if(_sylvan_needs_pic)
            set(LACE_ENABLE_PIC ON)
        endif()

        include(FetchContent)
        FetchContent_Declare(
            lace
            GIT_REPOSITORY https://github.com/trolando/lace.git
            GIT_TAG        v1.7.0
            GIT_SHALLOW    TRUE
        )
        FetchContent_MakeAvailable(lace)

        unset(BUILD_SHARED_LIBS)
        if(_sylvan_needs_pic)
            unset(LACE_ENABLE_PIC)
        endif()
    endif()
endif()

# Installed Lace packages use the namespaced target. Keep compatibility with a
# parent project that provides only the older unnamespaced target.
if(NOT TARGET lace::lace)
    if(TARGET lace)
        add_library(lace::lace ALIAS lace)
    else()
        message(FATAL_ERROR "Lace was found, but it did not provide a lace::lace target")
    endif()
endif()

unset(_sylvan_needs_pic)

# ── Sanitizer support (top-level builds only) ────────────────────────────────

if(PROJECT_IS_TOP_LEVEL)
    option(SYLVAN_SANITIZE_ADDRESS "Enable AddressSanitizer"           OFF)
    option(SYLVAN_SANITIZE_THREAD  "Enable ThreadSanitizer"            OFF)
    option(SYLVAN_SANITIZE_UB      "Enable UndefinedBehaviorSanitizer" OFF)
    option(SYLVAN_SANITIZE_MEMORY  "Enable MemorySanitizer"            OFF)

    if(SYLVAN_SANITIZE_THREAD AND SYLVAN_SANITIZE_ADDRESS)
        message(FATAL_ERROR "AddressSanitizer and ThreadSanitizer cannot be combined")
    endif()
    if(SYLVAN_SANITIZE_MEMORY AND SYLVAN_SANITIZE_ADDRESS)
        message(FATAL_ERROR "MemorySanitizer and AddressSanitizer cannot be combined")
    endif()
    if(SYLVAN_SANITIZE_MEMORY AND SYLVAN_SANITIZE_THREAD)
        message(FATAL_ERROR "MemorySanitizer and ThreadSanitizer cannot be combined")
    endif()

    set(_san "")
    if(SYLVAN_SANITIZE_ADDRESS)
        list(APPEND _san "address")
    endif()
    if(SYLVAN_SANITIZE_THREAD)
        list(APPEND _san "thread")
    endif()
    if(SYLVAN_SANITIZE_UB)
        list(APPEND _san "undefined")
    endif()
    if(SYLVAN_SANITIZE_MEMORY)
        list(APPEND _san "memory")
    endif()

    if(_san)
        if(NOT CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
            message(FATAL_ERROR "Sanitizers are not supported with compiler '${CMAKE_C_COMPILER_ID}'")
        endif()
        string(REPLACE ";" "," _san_list "${_san}")
        add_library(sylvan_sanitizers INTERFACE)
        target_compile_options(sylvan_sanitizers INTERFACE
            "-fsanitize=${_san_list}"
            -fno-omit-frame-pointer
        )
        target_link_options(sylvan_sanitizers INTERFACE "-fsanitize=${_san_list}")
        message(STATUS "Sanitizers enabled: ${_san_list}")
        unset(_san_list)
    endif()
    unset(_san)
endif()

# ── Hide options when consumed as a subproject ───────────────────────────────

if(NOT PROJECT_IS_TOP_LEVEL)
    mark_as_advanced(
        BUILD_SHARED_LIBS
        SYLVAN_ENABLE_PIC
        SYLVAN_BUILD_EXAMPLES
        SYLVAN_BUILD_DOCS
    )
endif()

# ── Subdirectories ───────────────────────────────────────────────────────────

# Add the Q-Sylvan library target
add_subdirectory(src)

if(TARGET sylvan_sanitizers)
    target_link_libraries(qsylvan PUBLIC "$<BUILD_INTERFACE:sylvan_sanitizers>")
endif()

# Add the QASM parser
add_subdirectory(qasm)

if(BUILD_TESTING)
    add_subdirectory(test)
endif()

if(SYLVAN_BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()

if(SYLVAN_BUILD_DOCS)
    configure_file("docs/conf.py.in" "docs/conf.py" @ONLY)
    find_package(Sphinx REQUIRED)
    Sphinx_add_targets(qsylvan
        ${CMAKE_CURRENT_BINARY_DIR}/docs
        ${CMAKE_CURRENT_SOURCE_DIR}/docs
        ${CMAKE_CURRENT_BINARY_DIR}
    )
    add_custom_target(update_gh_pages
        COMMAND "${CMAKE_COMMAND}" -P "${PROJECT_SOURCE_DIR}/cmake/UpdateGHPages.cmake"
    )
    add_dependencies(update_gh_pages qsylvan_html)
endif()
