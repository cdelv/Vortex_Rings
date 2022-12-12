/**
# Interface for the Abstract Syntax Tree (AST) library

It is divided in two parts: a first part for functions which are meant
to be independent from the grammar, and a second part for (Basilisk C)
grammar-dependent functions. 

## Grammar-independent functions */

#include <stdio.h>
#include <stdbool.h>
#include <assert.h>
#include "allocator.h"
#include "stack.h"

typedef struct _Ast Ast;

struct _Ast {
  int sym;
  Ast ** child, * parent;
};

typedef struct {
  Ast ast;
  char * before, * start, * after;
  const char * file;
  int line;
  void * value;
} AstTerminal;

typedef struct {
  Ast ast;
  char * before, * after;
  const char * file;
  Allocator * alloc;
  Stack * stack;
  bool type_already_specified;
} AstRoot;

AstRoot * ast_parse            (const char * code, AstRoot * parent);
Ast *     ast_parse_expression (const char * expr, AstRoot * parent);
Ast *     ast_parse_external_declaration (const char * decl, AstRoot * parent);
void      ast_destroy          (Ast * n);
AstRoot * ast_parse_file       (FILE * fp, AstRoot * parent);
void      ast_print            (Ast * n, FILE * fp, int kind);
void      ast_print_tree       (Ast * n, FILE * fp, const char * indent,
				bool compress, int maxdepth);
void      ast_print_file_line  (Ast * n, FILE * fp);
AstRoot * ast_get_root         (const Ast * n);
void      ast_identifier_print (Ast * identifier, FILE * fp);
void      ast_stack_print      (Stack * stack, FILE * fp);

static inline Ast * ast_last_child (const Ast * n)
{
  if (!n)
    return NULL;
  Ast ** c = n->child;
  if (!c)
    return NULL;
  while (*(c + 1)) c++;
  return *c;
}

extern Ast * const ast_placeholder;

static inline Ast * ast_child (const Ast * n, int sym)
{
  if (!n)
    return NULL;
  Ast ** c = n->child;
  if (c == NULL)
    return NULL;
  while (*c && (*c == ast_placeholder || (*c)->sym != sym)) c++;
  return *c;
}

static inline Ast * ast_parent (const Ast * n, int sym)
{
  if (!n || n == ast_placeholder)
    return NULL;
  Ast * parent = n->parent;
  while (parent && parent->sym != sym)
    parent = parent->parent;
  return parent;
}

static inline AstTerminal * ast_left_terminal (const Ast * n)
{
  while (n && n != ast_placeholder && n->child) {
    Ast ** c = n->child;
    while (*c && *c == ast_placeholder) c++;
    n = *c;
  }
  return (AstTerminal *) (n != ast_placeholder ? n : NULL);
}

static inline AstTerminal * ast_right_terminal (const Ast * n)
{
  while (n && n != ast_placeholder && n->child)
    n = ast_last_child (n);
  return (AstTerminal *) (n != ast_placeholder ? n : NULL);
}
AstTerminal * ast_next_terminal (const Ast * n);

AstTerminal * ast_terminal_new (Ast * parent, int symbol, const char * start);
#define ast_terminal_new_char(p,s)			\
  (Ast *) ast_terminal_new (p, token_symbol((s)[0]), s)

Ast * ast_new_internal (Ast * parent, ...);
#define ast_new(parent,...) ast_new_internal (parent, __VA_ARGS__, -1)
Ast * ast_schema_internal (const Ast * n, ...);
#define ast_schema(n,...) ast_schema_internal (n, __VA_ARGS__, -1)
Ast * ast_find_internal (const Ast * n, ...);
#define ast_find(n,...) ast_find_internal (n, __VA_ARGS__, -1)
Ast * ast_copy_single (const Ast * n,
		       AstRoot ** dst_root, AstRoot ** src_root);
Ast * ast_copy_internal (const Ast * n, ...);
#define ast_copy(n,...) ast_copy_internal (n, __VA_ARGS__ + 0, -1)
Ast * ast_attach_internal (Ast * parent, ...);
#define ast_attach(p,...) ast_attach_internal (p, __VA_ARGS__, NULL)
Ast * ast_new_children_internal (Ast * parent, ...);
#define ast_new_children(p,...) \
  ast_new_children_internal (p, __VA_ARGS__, NULL)

void ast_set_child (Ast * parent, int index, Ast * child);
void ast_replace_child (Ast * parent, int index, Ast * child);

static inline int ast_child_index (const Ast * n)
{
  assert (n->parent);
  int index = 0;
  Ast ** c;
  for (c = n->parent->child; *c && *c != n; c++, index++);
  return *c == n ? index : - 1;
}

static inline Ast * ast_ancestor (Ast * n, int i)
{
  while (n && i)
    n = n->parent, i--;
  return n;
}

char * str_append_realloc (char * dst, ...);
#define str_append(dst, ...)						\
  do { dst = str_append_realloc (dst, __VA_ARGS__, NULL); } while(0)
char * str_prepend_realloc (char * dst, ...);
#define str_prepend(dst, ...)						\
  do { dst = str_prepend_realloc (dst, __VA_ARGS__, NULL); } while(0)

#define ast_before(n,...) str_append(ast_left_terminal (n)->before, __VA_ARGS__)
#define ast_after(n,...)  str_append(ast_right_terminal (n)->after, __VA_ARGS__)
#define ast_terminal(n) ((n)->child ? NULL : (AstTerminal *)(n))
#define ast_root(n) ((n)->parent ? NULL : (AstRoot *)(n))
char *  ast_str_append (Ast * n, char * s);

static inline void ast_hide (AstTerminal * n)
{
  for (char * s = n->start; *s != '\0'; s++)
    *s = ' ';
}

char * ast_line (AstTerminal * t);
#define ast_file_line(n, nolineno)					\
  "\"", ast_terminal((Ast *)n)->file, "\",",				\
    nolineno ? "0" : ast_line(ast_terminal((Ast *)n))
void ast_set_line (Ast * n, AstTerminal * l);
Ast * ast_flatten (Ast * n, AstTerminal * t);
AstTerminal * ast_replace (Ast * n, const char * terminal, Ast * with);

#define NN(parent,sym,...) ast_new_children (ast_new (parent, sym), __VA_ARGS__)

static inline AstTerminal * NB (Ast * parent, int sym, const char * name)
{
  AstTerminal * t = ast_terminal_new (parent, sym, name);
  AstTerminal * r = ast_left_terminal (parent);
  t->before = r->before, r->before = NULL;
  return t;
}

static inline AstTerminal * NA (Ast * parent, int sym, const char * name)
{
  AstTerminal * t = ast_terminal_new (parent, sym, name);
  AstTerminal * r = ast_right_terminal (parent);
  t->line = r->line;
  return t;
}

#define NCB(parent,sym) NB(parent, token_symbol((sym)[0]), sym)
#define NCA(parent,sym) NA(parent, token_symbol((sym)[0]), sym)

/**
## Grammar-specific functions */

void  ast_push_declaration         (Stack * stack, Ast * declaration);
Ast * ast_push_function_definition (Stack * stack, Ast * declarator);
void  ast_pop_scope                (Stack * stack, Ast * scope);

void  ast_traverse (Ast * n, Stack * stack,
		    void func (Ast *, Stack *, void *),
		    void * data);
#define foreach_item(list, index, item)					\
  for (Ast * _list = list, * item = _list && _list != ast_placeholder ?	\
	 (_list->child[1] ? _list->child[index] : _list->child[0]) : NULL; \
       (_list = _list && _list != ast_placeholder &&			\
	_list->child[1] ? _list->child[0] : NULL), item;		\
       item = _list && _list != ast_placeholder ?			\
	 (_list->child[1] ? _list->child[index] :			\
	  _list->child[0]) : NULL					\
       )

Ast * ast_identifier_declaration (Stack * stack, const char * identifier);
Ast * ast_identifier_declaration_from_to (Stack * stack, const char * identifier,
					  const Ast * start, const Ast * end);
Ast * ast_function_identifier (const Ast * function_definition);
Ast * ast_function_call_identifier (const Ast * n);

void ast_set_char (Ast * n, int c);
void ast_remove_internal (Ast * n, AstTerminal * before);
void ast_remove (Ast * n, AstTerminal * before);
void ast_erase (Ast * n);
void ast_check (Ast * n);
Ast * ast_is_typedef (const Ast * identifier);
Ast * ast_find_function (Ast * n, const char * name);
Ast * ast_list_append_list (Ast * list, Ast * list1);
Ast * ast_block_list_append (Ast * list, int item_sym, Ast * item);
Ast * ast_block_list_insert_after (Ast * insert, Ast * item);
Ast * ast_block_list_insert_before (Ast * insert, Ast * item);
Ast * ast_block_list_insert_before2 (Ast * insert, Ast * item);
Ast * ast_block_list_get_item (Ast * statement);
Ast * ast_list_append (Ast * list, int item_sym, Ast * item);
Ast * ast_list_prepend (Ast * list, int item_sym, Ast * item);
Ast * ast_list_remove (Ast * list, Ast * item);
void  ast_argument_list (Ast * expression);
Ast * ast_new_unary_expression (Ast * parent);
Ast * ast_new_constant (Ast * parent, int symbol, const char * value);
Ast * ast_new_identifier (Ast * parent, const char * name);
Ast * ast_new_member_identifier (Ast * parent, const char * name);
Ast * ast_is_unary_expression (const Ast * n);
Ast * ast_is_identifier_expression (const Ast * n);
Ast * ast_is_simple_expression (const Ast * n);
Ast * ast_get_struct_name (Ast * declaration_specifiers);
bool  ast_are_identical (const Ast * a, const Ast * b);
Ast * ast_expression_type (Ast * expr, Stack * stack, bool higher_dimension);
AstTerminal * ast_type (const Ast * identifier);
char * ast_typedef_name (const Ast * identifier);

Ast * ast_check_grammar (Ast * n, bool recursive);

bool  ast_is_stencil_function (Ast * n);
Ast * ast_stencil (Ast * n, bool parallel, bool overflow, bool nowarning);
