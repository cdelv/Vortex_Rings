#include <stdlib.h>
#include <string.h>
#include "stack.h"

struct _Stack {
  void * p;
  int n, size;
};

Stack * stack_new (int size)
{
  Stack * s = calloc (1, sizeof (Stack));
  s->size = size;
  return s;
}

void stack_push (Stack * s, void * p)
{
  s->n++;
  s->p = realloc (s->p, s->n*s->size);
  char * dest = ((char *)s->p) + (s->n - 1)*s->size;
  memcpy (dest, p, s->size);
}

void * stack_pop (Stack * s)
{
  if (!s->n)
    return NULL;
  return ((char *)s->p) + --s->n*s->size;
}

void * stack_index (Stack * s, int i)
{
  if (!s->n || i > s->n - 1)
    return NULL;
  return ((char *)s->p) + (s->n - i - 1)*s->size;
}

void stack_destroy (Stack * s)
{
  free (s->p);
  free (s);
}
