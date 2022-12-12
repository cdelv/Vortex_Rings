%option noyywrap

%{

static int input (void);

static void comment (void)
{
  int c;  
  while ((c = input()) != 0) {

    if (c == '*') {
      while ((c = input()) == '*')
	;
      
      if (c == '/')
	return;
      
      if (c == 0)
	break;
    }
  }
}
 
int getput(void)
{
  int c = input();
  fputc (c, yyout);
  return c;
}
 
static int line = 1, indef = 0, nolineno = 0;
static char * autolink = NULL;

#define nonspace(s) { while (strchr(" \t\v\n\f", *s)) s++; }
#define space(s) { while (!strchr(" \t\v\n\f", *s)) s++; }
 
#define YY_INPUT(buf,result,max_size)			      \
  {							      \
    int c = fgetc(yyin);				      \
    if (c == '\n') { line++; }				      \
    result = (c == EOF) ? YY_NULL : (buf[0] = c, 1);	      \
  }
  
%}

SP     [ \t]
WS     [ \t\v\n\f]
ES     (\\([\'\"\?\\abfnrtv]|[0-7]{1,3}|x[a-fA-F0-9]+))
STRING \"([^"\\\n]|{ES})*\"

%%

"/*"                                    { comment(); }
"//".*                                  { ECHO; /* consume //-comment */ }

^{SP}*#{SP}*pragma{SP}+autolink{SP}+.*$ {
  char * s = strstr (yytext, "autolink"); space (s);
  if (!autolink) {
    autolink = malloc (strlen(s) + 1);
    autolink[0] = '\0';
  }
  else
    autolink = realloc (autolink, strlen(s) + 1 + strlen (autolink));
  strcat (autolink, s);
}
	   
^{SP}*#{SP}+[0-9]+{SP}+{STRING}({SP}+[0-9]+)* {
  /* replace # 3 "machin.h" 2 4 with #line 3 "machin.h" */
  char * s = strchr(yytext, '#') + 1;
  line = atoi(s);
  if (!indef) {
    fputs ("#line", yyout);
    char * s1 = strchr(s, '"');
    s1 = strchr(s1 + 1, '"');
    *(s1 + 1) = '\0';  
    fputs (s, yyout);
  }
}

^{SP}*@{SP}*def{SP}.*                   {
  /* replace @def ... @ with #define ... \\ */
  fputs ("#define", yyout);
  fputs (strstr (yytext, "def") + 3, yyout);
  indef = 1;
}

^{SP}*@.*" Pragma(" {
  yytext = strchr(yytext, '@'); yytext++;
  char * s = strstr (yytext, "Pragma("); *s++ = '\0';
  fprintf (yyout, "#%s", yytext);
  fputs ("_Pragma(", yyout);
  register int oldc = 0, c;
  for (;;) {
    while ((c = getput()) != '\n' && c != EOF)
      oldc = c;    /* eat up text of preproc */
    if (c == EOF || oldc != '\\')
      break;
  }
}
	   
@ {
  if (indef) {
    indef = 0;
    fprintf (yyout, "\n#line %d\n", line - 1);
  }
  else
    REJECT;
}

^{SP}*@{SP}*[^ \t\n] {
  /* replace @define etc. with #define etc. */
  *strchr(yytext, '@') = '#'; ECHO;
}

S__FILE__ {
  fputs ("__FILE__", yyout);
}

S_LINENO {
  fputs (nolineno ? "0" : "__LINE__", yyout);
}
	   
{STRING}+	{ ECHO; }
\n            {
  if (indef)
    fputs ("\\\n", yyout);
  else
    ECHO;
}
.		{ ECHO; }

%%

int postproc (FILE * fin, FILE * fout, char ** autolink1, int _nolineno)
{
  if (0) yyunput (0, NULL); // just prevents 'yyunput unused' compiler warning
  yyin = fin;
  yyout = fout;
  line = 1;
  indef = 0;
  nolineno = _nolineno;
  int status = yylex();
  *autolink1 = autolink;
  return status;
}
