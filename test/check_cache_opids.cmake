# Fails when two operation ids in sylvan_int.h agree.
#
# An id is OR-ed into the first key word of a cache entry and never checked,
# so two operations that share one read each other's entries whenever the
# rest of their keys coincide. That has happened twice (see the comment above
# CACHE_LIMDD_PLUS), and neither the compiler nor any functional test sees it.
#
#   cmake -DHEADER=<path to sylvan_int.h> -P check_cache_opids.cmake
file(READ "${HEADER}" text)
string(REGEX MATCHALL "static const uint64_t[ \t]+CACHE_[A-Z0-9_]+[ \t]*=[ \t]*\\(([0-9]+)LL<<40\\)" decls "${text}")
set(seen "")
set(bad 0)
foreach(d IN LISTS decls)
    string(REGEX REPLACE ".*(CACHE_[A-Z0-9_]+).*\\(([0-9]+)LL<<40\\).*" "\\1;\\2" pair "${d}")
    list(GET pair 0 name)
    list(GET pair 1 id)
    if(DEFINED owner_${id})
        message(SEND_ERROR "cache operation id ${id} is both ${owner_${id}} and ${name}")
        set(bad 1)
    else()
        set(owner_${id} ${name})
    endif()
    list(APPEND seen ${id})
endforeach()
list(LENGTH seen count)
if(count LESS 50)
    message(SEND_ERROR "found only ${count} cache operation ids in ${HEADER}; the pattern no longer matches the header")
endif()
if(NOT bad)
    message(STATUS "${count} cache operation ids, all distinct")
endif()
