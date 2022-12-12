#include <stdio.h>
#include <stdlib.h>
#include "ast.h"
#include "symbols.h"

static bool check (Ast * n)
{
  int sym_error = sym_YYerror;
  #include "grammar.h"
  return false;
}

Ast * ast_check_grammar (Ast * n, bool recursive)
{
  if (n && n->child) {    
    if (!check (n)) {
      fprintf (stderr, "\ngrammatical error:\n");
      ast_print_tree (n, stderr, 0, 0, 10);
      abort ();
    }
    if (recursive)
      for (Ast ** c = n->child; *c; c++) {
	assert ((*c)->parent == n);
	ast_check_grammar (*c, true);
      }
  }
  return n;
}
