/**
# Automatic stencils and boundary conditions

Basilisk automatically computes, at runtime, the access pattern
(i.e. "stencils") of (basic) foreach loops (foreach(), foreach_face(),
foreach_vertex()).

This is done in practice by `qcc` which automatically adds, before
each foreach loop, a minimal version of the loop body.

The resulting access pattern is stored in the `read` and `write`
arrays associated with each field.

The `dirty` attribute is used to store the status of boundary
conditions for each field. */

attribute {
  // fixme: use a structure
  bool input, output;
  int width; // maximum stencil width/height/depth
  int dirty; // // boundary conditions status:
  // 0: all conditions applied
  // 1: nothing applied
  // 2: boundary_face applied
}

typedef struct {
  const char * fname; // name of the source file
  int line;           // line number in the source
  int first;          // is this the first time the loop is called?
  int face;           // the face component(s) being traversed
  bool vertex;        // is this a vertex traversal?
} ForeachData;

// fixme: this should be rewritten using better macros
@def foreach_stencil() {
  static ForeachData _loop = {
    S__FILE__, S_LINENO,
    1, 0, 0
  };
  if (baseblock) for (scalar s = baseblock[0], * i = baseblock;
		s.i >= 0; i++, s = *i) {
    _attribute[s.i].input = _attribute[s.i].output = false;
    _attribute[s.i].width = 0;
  }
  int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);
  Point point = {0}; NOT_UNUSED (point);
@

@def end_foreach_stencil()
  end_stencil (&_loop);
  _loop.first = 0;
}
@

@define foreach_vertex_stencil() foreach_stencil() _loop.vertex = true;
@define end_foreach_vertex_stencil() end_foreach_stencil()

@define foreach_face_stencil() foreach_stencil()
@define end_foreach_face_stencil() end_foreach_stencil()

@define foreach_visible_stencil(...) foreach_stencil()
@define end_foreach_visible_stencil(...) end_foreach_stencil()

@define _stencil_is_face_x() { _loop.face |= (1 << 0);
@define end__stencil_is_face_x() }
@define _stencil_is_face_y() { _loop.face |= (1 << 1);
@define end__stencil_is_face_y() }
@define _stencil_is_face_z() { _loop.face |= (1 << 2);
@define end__stencil_is_face_z() }

void stencil_val (Point p, scalar s, int i, int j, int k,
		  const char * file, int line, bool overflow);
void stencil_val_a (Point p, scalar s, int i, int j, int k, bool input,
		    const char * file, int line);

@def _stencil_val(a,_i,_j,_k)
  stencil_val (point, a, _i, _j, _k, S__FILE__, S_LINENO, false)
@
@def _stencil_val_o(a,_i,_j,_k)
  stencil_val (point, a, _i, _j, _k, S__FILE__, S_LINENO, true)
@
@def _stencil_val_a(a,_i,_j,_k)
  stencil_val_a (point, a, _i, _j, _k, false, S__FILE__, S_LINENO)
@
@def _stencil_val_r(a,_i,_j,_k)
  stencil_val_a (point, a, _i, _j, _k, true, S__FILE__, S_LINENO)
@

@define _stencil_fine(a,_i,_j,_k) _stencil_val(a,_i,_j,_k)
@define _stencil_fine_a(a,_i,_j,_k) _stencil_val_a(a,_i,_j,_k)
@define _stencil_fine_r(a,_i,_j,_k) _stencil_val_r(a,_i,_j,_k)

@define _stencil_coarse(a,_i,_j,_k) _stencil_val(a,_i,_j,_k)
@define _stencil_coarse_a(a,_i,_j,_k) _stencil_val_a(a,_i,_j,_k)
@define _stencil_coarse_r(a,_i,_j,_k) _stencil_val_r(a,_i,_j,_k)

@define r_assign(x)
@define _assign(x)

@define _stencil_neighbor(i,j,k)
@define _stencil_child(i,j,k)
@define _stencil_aparent(i,j,k)
@define _stencil_aparent_a(i,j,k)
@define _stencil_aparent_r(i,j,k)

@define _stencil_neighborp(i,j,k) neighborp(i,j,k)

int _stencil_nop;
@define _stencil_val_higher_dimension (_stencil_nop = 1)
@define _stencil__val_constant(a,_i,_j,_k) (_stencil_nop = 1)

typedef void _stencil_undefined;

@define o_stencil -2

/**
## Automatic boundary conditions

Boundary conditions need to be applied if `s` is dirty, or if any of
the field `d` it depends on is dirty. */

static inline bool scalar_is_dirty (scalar s)
{
  if (s.dirty)
    return true;
  scalar * depends = s.depends;
  for (scalar d in depends)
    if (d.dirty)
      return true;
  return false;
}

/**
Does the boundary conditions on `a` depend on those on `b`? */

static inline bool scalar_depends_from (scalar a, scalar b)
{
  scalar * depends = a.depends;
  for (scalar s in depends)
    if (s.i == b.i)
      return true;
  return false;
}

/**
There are two types of boundary conditions: "full" boundary
conditions, done by `boundary_internal()` and "flux" boundary
conditions (i.e. normal components on faces only) done by
`boundary_face()`. */

void boundary_internal (scalar * list, const char * fname, int line);
void (* boundary_face)  (vectorl);

/**
This function is called after the stencil access detection, just
before the (real) foreach loop is executed. This is where we use the
stencil access pattern to see whether boundary conditions need to be
applied. */

void end_stencil (ForeachData * loop)
{
  scalar * listc = NULL, * dirty = NULL;
  vectorl listf = {NULL};
  bool flux = false;
  
  /**
  We check the accesses for each field... */
  
  for (scalar s in baseblock) {
    bool write = s.output, read = s.input;
    
#ifdef foreach_layer
    if (_layer == 0 || s.block == 1)
#endif
    {

      /**
      If the field is read and dirty, we need to check if boundary
      conditions need to be applied. */
      
      if (read && scalar_is_dirty (s)) {

	/**
	If this is a face field, we check whether "full" BCs need to
	be applied, or whether "face" BCs are sufficient. */
	
	if (s.face) {
	  if (s.width > 0) // face, stencil wider than 0
	    listc = list_append (listc, s);
	  else if (!write) { // face, flux only
	    scalar sn = s.v.x.i >= 0 ? s.v.x : s;
	    foreach_dimension()
	      if (s.v.x.i == s.i) {

		/* fixme: imposing BCs on fluxes should be done by
		   boundary_face() .*/
		
		if (sn.boundary[left] || sn.boundary[right])
		  listc = list_append (listc, s);
		else if (s.dirty != 2) {
		  listf.x = list_append (listf.x, s);
		  flux = true;
		}
	      }
	  }
	}

	/**
	For dirty, centered fields BCs need to be applied if the
	stencil is wider than zero. */
	
	else if (s.width > 0)
	  listc = list_append (listc, s);
      }

      /**
      Write accesses need to be consistent with the declared field
      type (i.e. face or vertex). */
      
      if (write) {
	if (dimension > 1 && !loop->vertex && loop->first) {
	  bool vertex = true;
	  foreach_dimension()
	    if (s.d.x != -1)
	      vertex = false;
	  if (vertex)
	    fprintf (stderr,
		     "%s:%d: warning: vertex scalar '%s' should be assigned with"
		     " a foreach_vertex() loop\n",
		     loop->fname, loop->line, s.name);
	}
	if (s.face) {
	  if (loop->face == 0 && loop->first)
	    fprintf (stderr,
		     "%s:%d: warning: face vector '%s' should be assigned with"
		     " a foreach_face() loop\n",
		     loop->fname, loop->line, s.name);
	}
	else if (loop->face) {
	  if (s.v.x.i < 0) { // scalar
	    int d = 1, i = 0;
	    foreach_dimension() {
	      if (loop->face == d) {
		s.face = 2, s.v.x.i = s.i;
		s.boundary[left] = s.boundary[right] = NULL;
#if PRINTBOUNDARY
		fprintf (stderr,
			 "%s:%d: turned %s into a face vector %c-component\n",
			 loop->fname, loop->line, s.name, 'x' + i);
#endif
	      }
	      d *= 2, i++;
	    }
	    if (!s.face && loop->first)
	      fprintf (stderr,
		       "%s:%d: warning: scalar '%s' should be assigned with "
		       "a foreach_face(x|y|z) loop\n",
		       loop->fname, loop->line, s.name);
	  }
	  else { // vector
	    char * name = NULL;
	    if (s.name) {
	      name = strdup (s.name);
	      char * s = name + strlen(name) - 1;
	      while (s != name && *s != '.') s--;
	      if (s != name) *s = '\0';
	    }
	    struct { int x, y, z; } input, output;
	    vector v = s.v;
#if 1 // fixme: should not be necessary	    
	    foreach_dimension()
	      input.x = v.x.input, output.x = v.x.output;
#endif
	    init_face_vector (v, name);
#if 1 // fixme: should not be necessary	    
	    
	    foreach_dimension()
	      v.x.input = input.x, v.x.output = output.x;
#endif
#if PRINTBOUNDARY
	    fprintf (stderr, "%s:%d: turned %s into a face vector\n",
		     loop->fname, loop->line, name);
#endif
	    free (name);
	  }
	}
	else if (loop->vertex) {
	  bool vertex = true;
	  foreach_dimension()
	    if (s.d.x != -1)
	      vertex = false;
	  if (!vertex) {
	    char * name = NULL;
	    if (s.name) name = strdup (s.name); // fixme: may not be necessary
	    init_vertex_scalar (s, name);
	    foreach_dimension()
	      s.v.x.i = -1;
#if PRINTBOUNDARY
	    fprintf (stderr, "%s:%d: turned %s into a vertex scalar\n",
		     loop->fname, loop->line, name);
#endif
	    free (name);
	  }
	}

	/**
	If the field is write-accessed, we add it to the 'dirty'
	list. */
	
	dirty = list_append (dirty, s);
	for (scalar d in baseblock)
	  if (scalar_depends_from (d, s))
	    dirty = list_append (dirty, d);
      }
    }
  }

  /**
  We apply face (flux) boundary conditions. */
  
  if (flux) {
#if PRINTBOUNDARY
    int i = 0;
    foreach_dimension() {
      if (listf.x) {
	fprintf (stderr, "%s:%d: flux %c:", loop->fname, loop->line, 'x' + i);
	for (scalar s in listf.x)
	  fprintf (stderr, " %d:%s", s.i, s.name);
	fputc ('\n', stderr);
      }
      i++;
    }
#endif
    boundary_face (listf);
    foreach_dimension()
      free (listf.x);
  }
  
  /**
  We apply "full" boundary conditions. */

  if (listc) {
#if PRINTBOUNDARY
    fprintf (stderr, "%s:%d: listc:", loop->fname, loop->line);
    for (scalar s in listc)
      fprintf (stderr, " %d:%s", s.i, s.name);
    fputc ('\n', stderr);
#endif
    boundary_internal (listc, loop->fname, loop->line);
    free (listc);
  }

  /**
  We update the dirty status of fields which will be write-accessed by
  the foreach loop. */
  
  if (dirty) {
#if PRINTBOUNDARY
    fprintf (stderr, "%s:%d: dirty:", loop->fname, loop->line);
    for (scalar s in dirty)
      fprintf (stderr, " %d:%s", s.i, s.name);
    fputc ('\n', stderr);
#endif
    for (scalar s in dirty)
      s.dirty = true;
    free (dirty);
  }
}

/**
## See also

* [Stencil test case](/src/test/stencils.c)
*/
