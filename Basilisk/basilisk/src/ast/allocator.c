#include <stdlib.h>
#include "allocator.h"

struct _Allocator {
  void ** m;
  int n;
  long len, maxlen;
};

Allocator * new_allocator()
{
  Allocator * a = calloc (1, sizeof (Allocator));
  a->maxlen = 1 << 16;
  a->m = malloc (sizeof (void *));
  a->m[0] = calloc (a->maxlen, 1);
  a->n = 1;
  return a;
}

void * allocate (Allocator * a, long size)
{
  if (a->len + size >= a->maxlen) {
    a->n++;
    a->m = realloc (a->m, a->n*sizeof (void *));
    a->m[a->n - 1] = calloc (a->maxlen, 1);
    a->len = 0;
  }
  void * p = (void *)(((char *)a->m[a->n - 1]) + a->len);
  a->len += size;
  return p;
}

void free_allocator (Allocator * a)
{
  for (int i = 0; i < a->n; i++)
    free (a->m[i]);
  free (a->m);
  free (a);
}
