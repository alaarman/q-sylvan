#ifndef QISQ2_MAP_H
#define QISQ2_MAP_H

/**
\brief Lockless non-resizing hash table implementation for fixed-length keys

@inproceedings{Laarman:2010:BMR:1998496.1998541,
  author = {Laarman, Alfons and van de Pol, Jaco and Weber, Michael},
  title = {{Boosting Multi-Core Reachability Performance with Shared Hash Tables}},
  booktitle = {Proceedings of the 2010 Conference on Formal Methods in Computer-Aided Design},
  series = {FMCAD '10},
  year = {2010},
  location = {Lugano, Switzerland},
  pages = {247--256},
  numpages = {10},
  url = {http://eprints.eemcs.utwente.nl/19281/},
  acmid = {1998541},
  publisher = {FMCAD Inc},
  address = {Austin, TX},
}
*/

#include <stdbool.h>
#include <stdint.h>
#include <string.h>
#include <gmp.h>

/**
\param Rationals (a,b,c,d) represent the number a + b sqrt(2) + c i + d i sqrt(2)
*/
typedef struct {
  mpq_t a;
  mpq_t b;
  mpq_t c;
  mpq_t d;
} qisq2_t;

extern uint32_t qisq2_hash(qisq2_t *num);

/**
\brief Create a new database.
\param len The length of the vectors to be stored here
\return the hashtable
*/
extern void *qisq2_map_create(uint64_t size, double tolerance);

/**
\brief Free the memory used by a dbs.
*/
extern void qisq2_map_free(void *dbs);

/**
\brief Free the memory used by an element in a dbs.
*/
extern void qisq2_map_free_elt(qisq2_t *num);

/**
\brief Find a vector with respect to a database and insert it if it cannot be fo
und.
\param dbs The dbs
\param vector The int vector
\retval idx The index that the vector was found or inserted at
\return 1 if the vector was present, 0 if it was added, -1 if table was full
*/
extern /* Defined here rather than in sylvan_edge_weights_qisq2.c: this library uses
 * them itself (find_or_put reduces before hashing), and qsylvan links it, so a
 * second copy up there is a duplicate symbol at link time. */
void qisq2_clear(qisq2_t *num);
void qisq2_reduce(qisq2_t *num);

int qisq2_map_find_or_put(const void *dbs, const void *v, uint64_t *ret);

extern void * qisq2_map_get(const void *dbs, const uint64_t ref);

extern uint64_t qisq2_map_count_entries(const void *dbs);

extern uint64_t qisq2_size(const void *dbs, const uint64_t ref);

extern void qisq2_map_print_bitvalues(const void *dbs, const uint64_t ref);

extern double qisq2_map_get_tolerance();

#endif // QISQ2_MAP
