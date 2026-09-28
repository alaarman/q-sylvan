#ifndef ATOMICS_H
#define ATOMICS_H

#include <errno.h>


/**
\file atomics.h
\brief Some operations for atomic data access
*/

/**
\brief An acquire load and a release store.

    These were volatile accesses, which stop the compiler from reordering but
    not the CPU, and the comment here said so: "Strong order execution has to
    be guaranteed by the CPU for this to work." x86 gives that; ARM, Apple
    Silicon included, does not. The weight tables publish an entry by writing
    its words and then its first word last, over a LOCK, and a reader waits
    for the LOCK to go and then compares. Without ordering, a reader on ARM
    could see the first word before the others, or, for qisq2, before the GMP
    limbs the words point at, compare against stale memory, find "not equal"
    and insert the same value a second time. Two indices for one number make
    a canonical form depend on the schedule: with four workers the exact
    LIMDD and BQD tests gave one state two edges, now and then, and the value
    was right each time. Release on the last write and acquire on the read
    order everything written before it, on every architecture; on x86 they
    compile to the same plain moves as before.
*/
#define atomic_read(v)      __atomic_load_n((v), __ATOMIC_ACQUIRE)
#define atomic_write(v,a)   __atomic_store_n((v), (a), __ATOMIC_RELEASE)

#define cas(a, b, c)        __sync_bool_compare_and_swap(a,b,c)
#define cas_ret(a, b, c)    __sync_val_compare_and_swap(a,b,c)
#define fetch_or(a, b)      __sync_fetch_and_or(a,b)
#define fetch_and(a, b)     __sync_fetch_and_and(a,b)
#define fetch_add(a, b)     __sync_fetch_and_add(a,b)
#define add_fetch(a, b)     __sync_add_and_fetch(a,b)
#define or_fetch(a, b)      __sync_or_and_fetch(a,b)
#define fetch_sub(a, b)     __sync_fetch_and_sub(a,b)
#define sub_fetch(a, b)     __sync_sub_and_fetch(a,b)
#define prefetch(a)         __builtin_prefetch(a)

#define mfence() { asm volatile("mfence" ::: "memory"); }

/* Compile read-write barrier */
#define compile_barrier() asm volatile("": : :"memory")

/* Pause instruction to prevent excess processor bus usage */
#define cpu_relax() asm volatile("pause\n": : :"memory")



/**
 * rwticket lock
 */

#ifndef PTHREAD_TICKET_RWLOCK
#define PTHREAD_TICKET_RWLOCK

typedef union rwticket ticket_rwlock_t;

union rwticket {
    unsigned u;
    unsigned short us;
    __extension__ struct {
        unsigned char write;
        unsigned char read;
        unsigned char users;
    } s;
};


static inline void
rwticket_init (ticket_rwlock_t *l)
{
    l->s.read = 0;
    l->s.write = 0;
    l->s.users = 0;
    l->u = 0;
    l->us = 0;
}

static inline void
rwticket_wrlock (ticket_rwlock_t *l)
{
    unsigned me = fetch_add(&l->u, (1<<16));
    unsigned char val = me >> 16;

    while (val != l->s.write) cpu_relax();
}

static inline void
rwticket_wrunlock (ticket_rwlock_t *l)
{
    ticket_rwlock_t t = *l;

    compile_barrier();

    t.s.write++;
    t.s.read++;

    *(unsigned short *) l = t.us;
}

static inline int
rwticket_wrtrylock (ticket_rwlock_t *l)
{
    unsigned me = l->s.users;
    unsigned char menew = me + 1;
    unsigned read = l->s.read << 8;
    unsigned cmp = (me << 16) + read + me;
    unsigned cmpnew = (menew << 16) + read + me;

    if (cas_ret(&l->u, cmp, cmpnew) == cmp) return 0;

    return EBUSY;
}

static inline void
rwticket_rdlock (ticket_rwlock_t *l)
{
    unsigned me = fetch_add(&l->u, (1<<16));
    unsigned char val = me >> 16;

    while (val != l->s.read) cpu_relax();
    l->s.read++;
}

static inline void
rwticket_rdunlock (ticket_rwlock_t *l)
{
    add_fetch (&l->s.write, 1);
}

static inline int
rwticket_rdtrylock (ticket_rwlock_t *l)
{
    unsigned me = l->s.users;
    unsigned write = l->s.write;
    unsigned char menew = me + 1;
    unsigned cmp = (me << 16) + (me << 8) + write;
    unsigned cmpnew = ((unsigned) menew << 16) + (menew << 8) + write;

    if (cas_ret(&l->u, cmp, cmpnew) == cmp) return 0;

    return EBUSY;
}

#endif


#endif
