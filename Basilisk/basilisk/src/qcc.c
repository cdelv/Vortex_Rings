/**
# The Basilisk C to C99 pre-processor

This is the front-end for the Basilisk C to C99 translator described
in [ast/README](). 

Usage:

~~~
qcc -grid=[GRID] [OPTIONS] FILE.c
~~~

A summary of the options/switches:

* `-grid=[GRID]` : specifies the grid to use
* `-MD` : generates .d dependency file
* `-tags` : generates .tags file
* `-python` : generates python wrapper code
* `-debug` : internal debugging
* `-events` : displays a trace of events on standard error
* `-catch` : catch floating point errors
* `-source` : generate C99 source file (with an underscore prefix)
* `-autolink` : use the 'autolink' pragma to link required libraries
* `-progress` : the running code will generate a 'progress' file
* `-cadna` : support for CADNA
* `-nolineno` : does not generate code containing code line numbers

All other options will be passed directly to the C compiler. */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <ctype.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <assert.h>

int dimension = 2, bghosts = 0, layers = 0;
  
int debug = 0, catch = 0, cadna = 0, nolineno = 0, events = 0, progress = 0;
int parallel = 0;
char dir[] = ".qccXXXXXX";

char * autolink = NULL;
int autolinks = 0, source = 0;
  
char * dname (const char * fname);
FILE * dopen (const char * fname, const char * mode);
  
void includes (int argc, char ** argv,
	       char ** grid, int * default_grid,
	       int * dimension, int * bg, int * layers,
	       const char * dir);

FILE * writepath (char * path, const char * mode)
{
  char * s = path;
  while ((s = strchr (s, '/'))) {
    *s = '\0';
    if (access (path, R_OK|W_OK|X_OK) && mkdir (path, 0700))
      return NULL;
    *s++ = '/';
  }
  return fopen (path, mode);
}

void cleanup (int status, const char * dir)
{
  if (!debug && dir) {
    char command[80] = "rm -r -f ";
    strcat (command, dir);
    int s = system (command); s = s;
  }
  exit (status);
}

char * dname (const char * fname)
{
  char * out = malloc (strlen (dir) + strlen (fname) + 2);
  strcpy (out, dir); strcat (out, "/"); strcat (out, fname);
  return out;
}

FILE * dopen (const char * fname, const char * mode)
{
  char * out = dname (fname);
  FILE * fout = fopen (out, mode);
  free (out);
  return fout;
}

void compdir (FILE * fin, FILE * fout, FILE * swigfp, 
	      char * swigname, char * grid)
{
  FILE * fout1 = dopen ("_endfor.c", "w");

  void endfor (FILE * fin, FILE * fout,
	       const char * grid, int dimension,
	       int nolineno, int progress, int catch, int parallel,
	       FILE * swigfp, char * swigname);

  endfor (fin, fout1, grid, dimension, nolineno, progress, catch, parallel,
	  swigfp, swigname);
  fclose (fout1);
  
  fout1 = dopen ("_endfor.c", "r");
  extern int postproc (FILE * fin, FILE * fout, char ** autolink, int nolineno);
  postproc (fout1, fout, &autolink, nolineno);
  fclose (fout1);
  fflush (fout);

  if (source && autolinks && autolink)
    printf ("%s\n", autolink);
}

int main (int argc, char ** argv)
{
  char * cc = getenv ("CC99"), command[1000], command1[1000] = "";
  if (cc == NULL)
    strcpy (command, CC99);
  else
    strcpy (command, cc);
  char * file = NULL;
  int i, dep = 0, tags = 0, swig = 0;
  for (i = 1; i < argc; i++) {
    if (!strncmp (argv[i], "-grid=", 6))
      ;
    else if (!strcmp (argv[i], "-MD"))
      dep = 1;
    else if (!strcmp (argv[i], "-tags"))
      tags = 1;
    else if (!strcmp (argv[i], "-python"))
      swig = 1;
    else if (!strcmp (argv[i], "-debug"))
      debug = 1;
    else if (!strcmp (argv[i], "-events"))
      events = 1;
    else if (!strcmp (argv[i], "-catch"))
      catch = 1;
    else if (!strcmp (argv[i], "-source"))
      source = 1;
    else if (!strcmp (argv[i], "-autolink"))
      autolinks = 1;
    else if (!strcmp (argv[i], "-progress"))
      progress = 1;
    else if (!strcmp (argv[i], "-Wall")) {
      char * s = strchr (command, ' ');
      if (s) {
	char command1[1000];
	strcpy (command1, s);
	*(s+1) = '\0';
	strcat (command, argv[i]);
	strcat (command, command1);
      }
      else
	strcat (command, argv[i]);
    }
    else if (!strcmp (argv[i], "-cadna")) {
      cadna = 1;
      char * cc = getenv ("CADNACC");
      if (cc == NULL)
	strcpy (command, CADNACC);
      else
	strcpy (command, cc);
    }
    else if (!strncmp (argv[i], "-Ddimension=", 12))
      dimension = 1 + argv[i][12] - '1';
    else if (catch && !strncmp (argv[i], "-O", 2))
      ;
    else if (!strcmp (argv[i], "-nolineno"))
      nolineno = 1;
    else if (!strcmp (argv[i], "-o")) {
      strcat (command1, " ");
      strcat (command1, argv[i++]);
      if (i < argc) {
	strcat (command1, " ");
	strcat (command1, argv[i]);
      }
    }
    else if (argv[i][0] != '-' && 
	     (tags || !strcmp (&argv[i][strlen(argv[i]) - 2], ".c"))) {
      if (file) {
	fprintf (stderr, "usage: qcc -grid=[GRID] [OPTIONS] FILE.c\n");
	return 1;
      }
      file = argv[i];
    }
    else if (!file) { 
      strcat (command, " ");
      strcat (command, argv[i]);
    }
    else {
      strcat (command1, " ");
      strcat (command1, argv[i]);
    }
  }
  if (strstr (command, "-D_MPI"))
    parallel = 1;
  char * openmp = strstr (command, "-fopenmp");
  if (openmp) {
    parallel = 1;
    if (strstr (command, "-D_MPI")) {
      fprintf (stderr,
	       "qcc: warning: OpenMP cannot be used with MPI (yet): "
	       "switching it off\n");
      int i;
      for (i = 0; i < strlen("-fopenmp"); i++)
	openmp[i] = ' ';
    }
    else if (swig) {
      fprintf (stderr,
	       "qcc: warning: OpenMP cannot be used with Python (yet): "
	       "switching it off\n");
      int i;
      for (i = 0; i < strlen("-fopenmp"); i++)
	openmp[i] = ' ';
    }
  }
  int status;
  if (debug) {
    status = system ("rm -r -f .qcc");
    strcpy (dir, ".qcc");
    status = mkdir (dir, 0700);
  }
  else
    status = (mkdtemp (dir) == NULL);
  if (status) {
    perror (dir);
    return 1;
  }
  if (file) {
    char * grid = NULL;
    int default_grid;
    includes (argc, argv, &grid, &default_grid,
	      &dimension, &bghosts, &layers,
	      dep || tags ? NULL : dir);
    FILE * swigfp = NULL;
    char swigname[80] = "";
    if (swig) {
      strcpy (swigname, file);
      char * dot = strchr (swigname, '.');
      *dot = '\0'; strcat (swigname, ".i");
      swigfp = fopen (swigname, "a");
      if (!swigfp) {
	fprintf (stderr, "qcc: could not open '%s': ", swigname);
	return 1;
      }
      *dot = '\0';
    }
    if (!dep && !tags) {
      char * basename = strdup (file), * ext = basename;
      while (*ext != '\0' && *ext != '.') ext++;
      char * cpp = malloc (strlen(basename) + strlen("-cpp") + strlen(ext) + 1);
      if (*ext == '.') {
	*ext = '\0';
	strcpy (cpp, basename);
	strcat (cpp, "-cpp");
	*ext = '.';
	strcat (cpp, ext);
      }
      else {
	strcpy (cpp, basename);
	strcat (cpp, "-cpp");
      }
      free (basename);
      FILE * fin = dopen (file, "r");
      if (!fin) {
	perror (file);
	cleanup (1, dir);
      }
      FILE * fp = dopen ("_attributes.h", "w");
      fputs ("typedef struct {\n", fp);
      fclose (fp);
      fp = dopen ("_maps.h", "w");
      fclose (fp);
      FILE * fout = dopen (cpp, "w");
      if (swig)
	fputs ("@include <Python.h>\n", fout);
      fputs ("@if _XOPEN_SOURCE < 700\n"
	     "  @undef _XOPEN_SOURCE\n"
	     "  @define _XOPEN_SOURCE 700\n"
	     "@endif\n"
	     "@if _GNU_SOURCE\n"
	     "@include <stdint.h>\n"
	     "@include <string.h>\n"
	     "@include <fenv.h>\n"
	     "@endif\n",
	     fout);
      if (catch)
	fputs ("#define TRASH 1\n"
	       "#define _CATCH last_point = point;\n", 
	       fout);
      else
	fputs ("#define _CATCH\n", fout);
      fprintf (fout, "#define dimension %d\n", dimension);
      if (bghosts)
	fprintf (fout, "#define BGHOSTS %d\n", bghosts);
      if (layers)
	fprintf (fout, "#define LAYERS 1\n");
      fputs ("#include \"common.h\"\n", fout);
      /* catch */
      if (catch)
	fputs ("void catch_fpe (void);\n", fout);
      /* grid */
      if (default_grid)
	fprintf (fout, "#include \"grid/%s.h\"\n", grid);
      char s[81];
      while (fgets (s, 81, fin))
	fputs (s, fout);
      if (swigfp)
	fputs ("#include \"python.h\"\n", fout);
      if (progress)
	fputs ("#include \"grid/progress.h\"\n", fout);
      fclose (fout);
      fclose (fin);
      fout = dopen (file, "w");
      if (!fout) {
	perror (file);
	cleanup (1, dir);
      }

      char preproc[1000], * cppcommand = getenv ("CPP99");
      strcpy (preproc, "cd ");
      strcat (preproc, dir);
      strcat (preproc, " && ");
      if (!cppcommand && strcmp (CPP99, ""))
	cppcommand = CPP99;
      if (!cppcommand) {
	if (source) {
	  /* remove -D_GNU_SOURCE flags from $CC99 */
	  char tmp[1000]; strcpy (tmp, command);
	  char * s = strtok (tmp, " \t");
	  while (s) {
	    if (strncmp (s, "-D_GNU_SOURCE", 13)) {
	      strcat (preproc, s);
	      strcat (preproc, " ");
	    }
	    s = strtok (NULL, " \t");
	  }
	}
	else
	  strcat (preproc, command);
	strcat (preproc, " -E");
      }
      else
	strcat (preproc, cppcommand);
      // remove "-pedantic option from preprocessor
      // This is mostly to avoid the preprocessor warning:
      // "ISO C99 requires at least one argument in a variadic macro"
      // note that the option is kept for the final compilation
      char * pedantic = strstr (preproc, "-pedantic");
      if (pedantic)
	for (i = 0; i < strlen ("-pedantic"); i++)
	  pedantic[i] = ' ';
      strcat (preproc, " -I");
      strcat (preproc, BASILISK);
      strcat (preproc, "/ast/std");
      strcat (preproc, " -I. -I");
      strcat (preproc, LIBDIR);
      strcat (preproc, " ");
      if (events) {
	strcat (preproc, " -DDEBUG_EVENTS=1 -DBASILISK=\"\\\"");
	strcat (preproc, BASILISK);
	strcat (preproc, "\\\"\" ");
      }
      if (nolineno)
	strcat (preproc, " -D'LINENO=0' ");
      strcat (preproc, cpp);
      free (cpp);
      if (debug) {
	fprintf (stderr, "preproc: %s\n", preproc);
	strcat (preproc, " | tee _preproc.c");
      }

      fin = popen (preproc, "r");
      if (!fin) {
	fclose (fout);
	perror (preproc);
	cleanup (1, dir);
      }

      compdir (fin, fout, swigfp, swigname, grid);
      int status = pclose (fin);
      fclose (fout);
      if (status == -1 ||
	  (WIFSIGNALED (status) && (WTERMSIG (status) == SIGINT || 
				    WTERMSIG (status) == SIGQUIT)))
	cleanup (1, dir);

      fout = dopen ("_tmp", "w");
      fin = dopen (file, "r");
      char line[1024];
      // rest of the file
      while (fgets (line, 1024, fin))
	fputs (line, fout);
      fclose (fin);
      fclose (fout);

      char src[80], dst[80];
      strcpy (src, dir); strcat (src, "/_tmp");
      if (source) {
	strcpy (dst, "_");
      }
      else {
	strcpy (dst, dir); strcat (dst, "/");
      }
      strcat (dst, file);
      rename (src, dst);

      strcat (command, " -I");
      strcat (command, LIBDIR);
      strcat (command, " ");
      strcat (command, dir);
      strcat (command, "/");
      strcat (command, file); 
      strcat (command, command1);
    }
  }
  else if (dep || tags) {
    fprintf (stderr, "usage: qcc -grid=[GRID] [OPTIONS] FILE.c\n");
    cleanup (1, dir);
  }
  else
    strcat (command, command1);
  /* compilation */
  if (!dep && !tags && !source) {
    if (autolinks && autolink)
      strcat (command, autolink);
    if (debug)
      fprintf (stderr, "command: %s\n", command);
    status = system (command);
    if (status == -1 ||
	(WIFSIGNALED (status) && (WTERMSIG (status) == SIGINT || 
				  WTERMSIG (status) == SIGQUIT)))
      cleanup (1, dir);
    cleanup (WEXITSTATUS (status), dir);
  }
  cleanup (0, dir);
  return 0;
}
