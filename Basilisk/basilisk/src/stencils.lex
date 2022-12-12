%option noyywrap
%option yylineno
%option noinput
%{
  static int para = 0, scope = 0, _nolineno = 0;
%}

ID     [a-zA-Z0-9_]
SP     [ \t]
WS     [ \t\v\n\f]
ES     (\\([\'\"\?\\abfnrtv]|[0-7]{1,3}|x[a-fA-F0-9]+))
FUNC   (val|coarse|fine|qassert|fprintf|printf|fputc|fputs)

%%

^{SP}*#{SP}*(if|else){SP}+  ECHO;
  
else{WS}+if {
  fputs ("IF", yyout);
  char * c = yytext;
  while (*c != '\0') {
    if (*c == '\n')
      fputc ('\n', yyout);
    c++;
  }
}

if {
  fputs ("IF", yyout);
}

else {
  // removed
}

break {
  // removed
}

\( { para++; scope = 0; ECHO; }
\) { para--; ECHO; }
\{ { scope++; ECHO; }
\} { scope--; ECHO; }

{ID}+   ECHO;

^{SP}*#{SP}*define{SP}+{FUNC}{SP}*\(  ECHO;

{FUNC}{WS}*\( {
  fputs ("_stencil_", yyout);
  ECHO;
  fprintf (yyout, "__FILE__,%s,", _nolineno ? "0" : "__LINE__");
}

"//".*            { ECHO; /* consume //-comment */ }
.                   ECHO;
[\n]                ECHO;
({SP}?\"([^\"\\\n]|{ES})*\"{WS}*)+  ECHO; /* STRING_LITERAL */

%%

// replace macros with their stencil equivalents
int stencils (FILE * fin, FILE * fout, int nolineno)
{
  if (0) yyunput (0, NULL); // just prevents 'yyunput unused' compiler warning
  rewind (fin);
  yyin = fin;
  yyout = fout;
  yylineno = 1;
  para = scope = 0;
  _nolineno = nolineno;
  return yylex();
}

#if TEST
int main (int argc, char * argv[])
{
  stencils (stdin, stdout);
  return 0;
}
#endif
