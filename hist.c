/*
 * Copyright (c) 2017
 *
 *     Yuan Mei and Cheng Zhang
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 * notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 * 3. Users of source code or binary forms are encouraged to take the
 * challenge of adopting a vegetarian diet (but eggs and dairy products
 * are permissible) for a full day chosen at their convenience as a token
 * of appreciation.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND CONTRIBUTORS ``AS IS'' AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHOR OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
 * OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE.
 */
/*
 * Inspired by the syntax of the TTree::Draw expression in CERN ROOT.
 * For example
 *
 * hist 'sin($1+$2^(log(if($3>1,$3,1)))):$2::exp(-$5)>>h(100, -1.0, 1.0, 50, 0.0, 2.0)' data.dat
 *            X-axis expr                 Y    W          nx  xmin xmax  ny ymin ymax
 *
 * will generate a 2D histogram according to the x and y binning
 * specifications using the data from the file.  The values of x, y, and
 * histogram weight are evaluated from the above X, Y, and W expressions,
 * respectively.  The weight is assumed to be 1 if ::W is absent.  1, 2 and 3
 * dimensional histograms are supported.  The generated histogram data
 * are printed to stdout and are suitable for plotting in gnuplot:
 *
 * plot "<hist 'sin($1)>>h(...)'" w step
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <limits.h>
#include <time.h>

struct token_s;

/** Parameters settable from commandline and some global states */
typedef struct param
{
  int normalize;       /**< normalize histogram instead of raw counts */
  int bincenter;       /**< print bin center instead of lower edge */
  const char *command; /**< histogram command string */
  const char *dfname;  /**< input data file name */
  /* to be filled by the parser from the input command */
  int dim;
  char *sx, *sy, *sz, *sw, *sh, *sxyz;
  int xn, yn, zn, colmax;
  double xmin, xmax, dx;
  double ymin, ymax, dy;
  double zmin, zmax, dz;
  struct token_s *px, *py, *pz, *pw;
} param_t;

/** set default parameters, C89 compatible */
static void param_set_default(param_t *pm)
{
  pm->normalize = 0;
  pm->bincenter = 0;
  pm->command = "sin(-$1+$2):cos($3)>>h(100, -1.0, 1.0, -23, -1.1, 1.2)";
  pm->dfname = "data.dat";
  pm->dim = 1;
  pm->colmax = 0;
  pm->sx = pm->sy = pm->sz = pm->sw = pm->sh = pm->sxyz = NULL;
  pm->px = pm->py = pm->pz = pm->pw = NULL;
}

/** free memory */
static void param_free(param_t *pm)
{
  free(pm->sxyz);
  free(pm->px); free(pm->py); free(pm->pz); free(pm->pw);
}

static void print_usage(const param_t *pm)
{
  fprintf(stderr, "Usage:\n");
  fprintf(stderr, "    hist -c -m 0 'exprX:exprY:exprZ::exprW>>h(nx,xl,xh,ny...)' dfname\n");
  fprintf(stderr, "         -c print bin center, default is bin lower edge, currently %s.\n", (pm->bincenter ? "on" : "off"));
  fprintf(stderr, "         -m normalization, 0: raw counts, 1: make the sum=1.0, 2: make the integral=1.0, currently %d.\n", pm->normalize);
  fprintf(stderr, "         exprX,Y,Z : expressions for computing X,Y,Z values\n");
  fprintf(stderr, "                     use $1 to address the leftmost data column in file.\n");
  fprintf(stderr, "         exprW expression for weight, if omitted, weight=1.\n");
  fprintf(stderr, "         Up to 3-dimensional histograms are supported.\n");
}

static int param_parse(param_t *pm, int argc, char **argv)
{
  int i, ich, iarg = 0;
  const char *p;

  param_set_default(pm);
  /* parse switches */
  for ( i = 1; i < argc; i++ ) {
    if ( argv[i][0] == '-' && strstr(argv[i], ">>") == NULL ) {
      /* going over characters in this options, s.t. "-cm" <==> "-c -m" */
      for ( ich = 1; ich > 0 && argv[i][ich] != '\0'; ich++ ) {
        switch ( argv[i][ich] ) {
        case 'c':
          pm->bincenter = 1;
          break;
        case 'm':
          p = argv[i] + ich + 1; /* string next to '-m' */
          if ( *p == '\0' && i < argc - 1 && isdigit(*argv[i+1]) ) {
            ich = -1; /* skip the remaining of this argument */
            p = argv[++i];
          }
          if ( isdigit(*p) ) {
            char *q;
            pm->normalize = (int) strtol(p, &q, 0);
            if ( ich >= 0 ) ich = q - argv[i] - 1;
          } else {
            pm->normalize++; /* such that "-mm" means "-m2" */
          }
          break;
        default:
          fprintf(stderr, "Error: unknown option [%c] in argument %d [%s]\n", argv[i][ich], i, argv[i]);
          print_usage(pm);
          return EXIT_FAILURE;
        }
      }
    } else {
      if ( iarg == 0 ) {
        pm->command = argv[i];
      } else if ( iarg == 1 ) {
        pm->dfname = argv[i];
      }
      iarg++;
    }
  }

  if ( iarg < 2
       || ( pm->normalize < 0 || pm->normalize > 2 ) ) {
    print_usage(pm);
    return EXIT_FAILURE;
  }
#ifdef DEBUG
  fprintf(stderr, "Input file:    %s\n"
          "Command:       %s\n"
          "Bin-center:    %s\n"
          "Normalization: %d\n",
          pm->dfname, pm->command, (pm->bincenter ? "on" : "off"), pm->normalize);
#endif
  return EXIT_SUCCESS;
}

/** expression token types */
enum { NULLTYPE, NUMBER, OPERATOR, FUNCTION, COLUMN, VARIABLE };

/* supported operators:
 * + and - (binary or unary), *, /, % (remainder), ^ or ** (power),
 * <, <=, >, >=, ==, !=, && (logical AND), || (logical OR) and parentheses () */
const char *ops = ",+-*/%^_<>!=&|()";

#define VARNAME_MAX 128

typedef struct token_s {
  int type;
  char s[VARNAME_MAX]; /**< operator/function name */
  double val;
  int col;             /**< column id in the data */
} token_t;

/** Get a token from the string and return a pointer after the token.
 * @param[out] t an identified token, modified.
 * @param[in] s expression string, not modified.
 * @return a pointer pointing to the character in string s right after the identified token.
 */
static char *hist_expr_token_get(token_t *t, const char *s)
{
  char *p;
  int i;

  /* skip leading spaces */
  while ( *s && isspace(*s) ) s++;
  if ( *s == '\0' ) return NULL;

  t->s[0] = '\0';
  t->val = 0;
  t->col = 0;

  if ( strchr(ops, *s) != NULL ) {
    t->type = OPERATOR;
    t->s[0] = *s;
    t->s[1] = '\0';
    p = (char *)s + 1;
    if ( *s == '*' && s[1] == '*' ) {
      t->s[0] = '^'; /* convert "**" to "^" */
      p++;
    }
  } else if ( isdigit(*s) || *s == '.' ) { /* .5 is understood as a number */
    t->type = NUMBER;
    t->val = strtod(s, &p);
  } else if ( isalpha(*s) ) {
    for ( p = (char*)s, i = 0; isalnum(*p) && i < VARNAME_MAX - 1; p++, i++ ) {
      t->s[i] = *p;
    }
    t->s[i] = '\0';
    while ( *p && isspace(*p) ) p++;
    if ( *p == '(' ) {
      t->type = FUNCTION;
      p++;
    } else {
      t->type = VARIABLE;
    }
  } else if ( *s == '$' ) {
    t->type = COLUMN;
    t->col = (int)strtol(s + 1, &p, 10);
  } else {
    fprintf(stderr, "Error: unknown token [%s]\n", s);
    exit(1);
  }
  return p;
}

#ifdef DEBUG
/** String representation of token
 * @param[in] tok
 * @param[out] s string representation of token.  Must be large enough to hold the string.
 * @return s
 */
static char *hist_expr_token2str(char *s, const token_t *tok)
{
  if ( tok->type == NUMBER ) {
    sprintf(s, "%g", tok->val);
  } else if ( tok->type == OPERATOR || tok->type == FUNCTION || tok->type == VARIABLE ) {
    sprintf(s, "%s", tok->s);
  } else if ( tok->type == COLUMN ) {
    sprintf(s, "$%d", tok->col);
  }
  return s;
}
#endif

static void hist_expr_token_copy(token_t *a, const token_t *b)
{
  a->type = b->type;
  strncpy(a->s, b->s, VARNAME_MAX);
  a->val = b->val;
  a->col = b->col;
}

/** Operator precedence table.
 * @return the precedence of the operator.
 */
static int hist_expr_preced(const char *s)
{
  static char prtable[256];
  static struct { const char *s; int p; } prlist[] = {
    { "<=", 9 }, { ">=", 9 },
    { "!=", 8 }, { "==", 8 }, { "<>", 8 },
    { "&&", 5 },
    { "||", 4 },
    { NULL, 0 } }, *prptr;

  if ( prtable['+'] == 0 ) { /* initialize the precedence table */
    prtable['('] = prtable[')'] = 16;
    prtable['_'] = prtable['!'] = 15;
    prtable['^'] = 13;
    prtable['*'] = prtable['/'] = prtable['%'] = 12;
    prtable['+'] = prtable['-'] = 11;
    prtable['<'] = prtable['>'] = 9;
    prtable[','] = 1;
  }

  if ( isalpha(*s) ) { /* functions take the highest precedence */
    return 1000;
  } else if ( s[1] == '\0' ) {
    return prtable[(int)(*s)];
  } else { /* search the list */
    for ( prptr = prlist; prptr->s != NULL; prptr++ ) {
      if ( strcmp(prptr->s, s) == 0 ) {
        return prptr->p;
      }
    }
  }
  fprintf(stderr, "Error: unknown operator [%s]\n", s);
  return 0;
}

static int hist_expr_isleftassoc(const char *s)
{
  if ( *s == '^' ) return 0;
  return 1;
}

/** Convert an expression to a stack of postfix expression.
 * The shunting-yard algorithm:
 * https://en.wikipedia.org/wiki/Shunting-yard_algorithm
 *
 * @return a malloc-ed queue stack of tokens.
 */
static token_t *hist_expr_parse2postfix(const char *s)
{
  token_t *que; /* output queue */
  token_t *ost; /* operator stack */
  token_t *pos, *top, tok[2];
  char *p;
  int n;

  n = strlen(s);

  /* initialize the output queue */
  if ( (que = calloc(n + 1, sizeof(*que))) == NULL ) exit(-1);
  pos = que;

  /* initialize the operator stack */
  if ( (ost = calloc(n + 1, sizeof(*ost))) == NULL ) exit(-1);
  top = ost;
  strcpy(ost[0].s, ""); /* special operator */

  tok[1].type = NULLTYPE; /* for the previous token */

  /* when there are tokens to be read, read a token */
  for ( p = (char*)s; (p = hist_expr_token_get(tok, p)) != NULL; ) {
    if ( tok->type == NUMBER || tok->type == COLUMN || tok->type == VARIABLE ) {
      /* if the token is a number, then push it to the output queue */
      hist_expr_token_copy(pos++, tok);
    } else if ( tok->s[0] == '(' || tok->type == FUNCTION ) {
      /* push "(" or a function onto the operator stack */
      hist_expr_token_copy(++top, tok); /* top = token */
    } else if ( tok->s[0] == ')' ) {
      /* while the operator at the top of the operator stack is not a "("
       * or a function */
      while ( top->s[0] != '(' && top->type != FUNCTION ) {
        /* pop operators from the operator stack onto the output queue */
        hist_expr_token_copy(pos++, top--);
        if ( top <= ost ) break;
      }
      if ( top->s[0] == '(' ) {
        top--; /* pop the left bracket from the stack */
      } else {
        hist_expr_token_copy(pos++, top--);
      }
    } else if ( tok->type == OPERATOR ) {
      if ( (tok->s[0] == '+' || tok->s[0] == '-') &&
           ( tok[1].type == NULLTYPE || tok[1].type == FUNCTION
             || (tok[1].type == OPERATOR && tok[1].s[0] != ')') ) ) {
        /* handle a unary operator */
        if ( tok->s[0] == '-' ) {
          tok->s[0] = '_'; /* change it to negation */
          hist_expr_token_copy(++top, tok); /* top = token */
        } /* skip the unary '+' */
      } else {
        while ( ((hist_expr_preced(top->s) > hist_expr_preced(tok->s))
                 || (hist_expr_preced(top->s) == hist_expr_preced(tok->s)
                     && hist_expr_isleftassoc(top->s)))
                && top->s[0] != '(' && top->type != FUNCTION ) {
          /* pop operators from the operator stack onto the output queue */
          hist_expr_token_copy(pos++, top--); /* pos = top */
          if ( top <= ost ) break;
        }
        /* push the read operator onto the operator stack */
        hist_expr_token_copy(++top, tok); /* top = token */
      }
    }
#ifdef DEBUG
    { /* debug routine to print out the output queue and operator stack */
      token_t *t;
      char buf[VARNAME_MAX];
      fprintf(stderr, "String: %s; Token: %s\nQueue: ", s, hist_expr_token2str(buf, tok));
      for ( t = que; t->type != NULLTYPE; t++ )
        fprintf(stderr, "%s ", hist_expr_token2str(buf, t));
      fprintf(stderr, "\nStack: ");
      for ( t = top; t > ost; t-- )
        fprintf(stderr, "%s ", hist_expr_token2str(buf, t));
      fprintf(stderr, "\n\n");
    }
#endif
    hist_expr_token_copy(tok+1, tok); /* make a copy */
  }

  /* while there are still operator tokens on the stack,
   * pop the operator onto the output queue */
  while ( top > ost )
    hist_expr_token_copy(pos++, top--);

  pos->type = NULLTYPE;
  free(ost);
  return que;
}

/* Tausworthe random number generator */
static double rand01(void)
{
  static unsigned long s1, s2, s3;

#if ( ULONG_MAX > 4294967295 )  /*  64 bit version */
#define MASK 0xffffffffUL
#define TAUSWORTHE(s, a, b, c, d) (((s & c) << d) & MASK) ^ ((((s << a) & MASK) ^ s) >> b)
#else /* 32 bit version */
#define TAUSWORTHE(s, a, b, c, d) (((s & c) << d) ^ (((s << a) ^ s) >> b) )
#endif
#define TAUS3() s1 = TAUSWORTHE(s1, 13, 19, 4294967294UL, 12); \
                s2 = TAUSWORTHE(s2,  2, 25, 4294967288UL,  4); \
                s3 = TAUSWORTHE(s3,  3, 11, 4294967280UL, 17)

  if ( s1 == 0 ) {
    s3 = (unsigned long) clock() & 0xffffffffUL;
    s1 = ( 1664525 * s3 + 1013904223 ) & 0xffffffffUL;
    s2 = ( 1664525 * s1 + 1013904223 ) & 0xffffffffUL;
    s3 = ( 1664525 * s2 + 1013904223 ) & 0xffffffffUL;
    TAUS3(); TAUS3(); TAUS3(); TAUS3(); TAUS3(); TAUS3();
  }
  TAUS3();
  return ( s1 ^ s2 ^ s3 ) / 4294967296.0;
#undef TAUSWORTHE
#undef TAUS3
}

static double dblmax(double a, double b) { return ( a > b ) ? a : b; }
static double dblmin(double a, double b) { return ( a < b ) ? a : b; }
static double iif(double c, double a, double b) { return ( c != 0 ) ? a : b; }

typedef struct {
  char s[VARNAME_MAX];
  double (*f)();
  int narg;
} funcmap_t;

typedef struct {
  char s[VARNAME_MAX];
  double val;
} varmap_t;

static funcmap_t funcmap[] = {
  {"rand", rand01, 0},
  {"rand01", rand01, 0},
  {"rand0_1", rand01, 0},
  {"fabs", fabs, 1},
  {"abs", fabs, 1},
  {"sqrt", sqrt, 1},
  {"exp", exp, 1},
  {"sinh", sinh, 1},
  {"cosh", cosh, 1},
  {"tanh", tanh, 1},
  {"log", log, 1},
  {"ln", log, 1},
  {"log10", log10, 1},
  {"pow", pow, 2},
  {"sin", sin, 1},
  {"cos", cos, 1},
  {"tan", tan, 1},
  {"asin", asin, 1},
  {"acos", acos, 1},
  {"atan", atan, 1},
  {"atan2", atan2, 2},
  {"ceil", ceil, 1},
  {"floor", floor, 1},
  {"fmod", fmod, 2},
  {"min", dblmin, 2},
  {"max", dblmax, 2},
  {"if", iif, 3},
  {"iif", iif, 3},
  {"", NULL, 0}
};

#ifndef M_PI
#define M_PI   3.14159265358979323846264338327950288
#endif

#ifndef M_E
#define M_E    2.71828182845904523536028747135266250
#endif

static varmap_t varmap[] = {
  {"pi", M_PI},
  {"Pi", M_PI},
  {"PI", M_PI},
  {"e",  M_E},
  {"E", M_E},
  {"", 0}
};

/** Evaluate the postfix expression stack.
 */
static double hist_expr_eval_postfix(const token_t *que, const double *arr, int narr)
{
  int i, n, top;
  double *st, ans;
  const token_t *pos;

  /* determine the length of the expression */
  for ( n = 0; que[n].type != NULLTYPE; n++ ) ;
  /* allocate the evaluation stack */
  if ( (st = calloc(n, sizeof(*st))) == NULL ) exit(-1);
  top = 0;

  for ( pos = que; pos->type != NULLTYPE; pos++ ) {
    if ( pos->type == NUMBER ) {
      st[top++] = pos->val;
    } else if ( pos->type == COLUMN ) {
      st[top++] = ( pos->col <= narr ) ? arr[pos->col] : 0.0;
    } else if ( pos->type == VARIABLE ) {
      for ( i = 0; varmap[i].s[0] != '\0'; i++ ) {
        if ( strcmp(varmap[i].s, pos->s) == 0 ) {
          st[top++] = varmap[i].val;
          break;
        }
      }
    } else if ( pos->type == OPERATOR ) {
      if ( strchr("_!", pos->s[0]) != NULL && pos->s[1] == '\0' ) { /* unary operators */
        if ( top < 1 ) {
          fprintf(stderr, "Error: insufficient arguments for operator [%s]\n",
                  pos->s);
          exit(1);
        }
        if ( pos->s[0] == '_' ) {
          st[top-1] = -st[top-1];
        } else if ( pos->s[0] == '!' ) {
          st[top-1] = !( st[top-1] != 0 );
        }
      } else if ( pos->s[0] == ',' ) { /* do nothing for comma */
        ;
      } else { /* binary operators */
        if ( top < 2 ) {
          fprintf(stderr, "Error: insufficient arguments for operator [%s], has %d\n",
                  pos->s, top);
          exit(1);
        }
        --top;
        if ( pos->s[0] == '+' ) {
          st[top-1] += st[top];
        } else if ( pos->s[0] == '-' ) {
          st[top-1] -= st[top];
        } else if ( pos->s[0] == '*' ) {
          st[top-1] *= st[top];
        } else if ( pos->s[0] == '/' ) {
          st[top-1] /= st[top];
        } else if ( pos->s[0] == '%' ) {
          st[top-1] = fmod(st[top-1], st[top]);
        } else if ( pos->s[0] == '^' ) {
          st[top-1] = pow(st[top-1], st[top]);
        } else if ( strcmp(pos->s, "<") == 0 ) {
          st[top-1] = ( st[top-1] < st[top] );
        } else if ( strcmp(pos->s, ">") == 0 ) {
          st[top-1] = ( st[top-1] > st[top] );
        } else if ( strcmp(pos->s, "<=") == 0 ) {
          st[top-1] = ( st[top-1] <= st[top] );
        } else if ( strcmp(pos->s, ">=") == 0 ) {
          st[top-1] = ( st[top-1] >= st[top] );
        } else if ( strcmp(pos->s, "!=") == 0 || strcmp(pos->s, "<>") == 0 ) {
          st[top-1] = ( st[top-1] != st[top] );
        } else if ( strcmp(pos->s, "==") == 0 ) {
          st[top-1] = ( st[top-1] == st[top] );
        } else if ( strcmp(pos->s, "&&") == 0 ) {
          st[top-1] = ( (st[top-1] != 0) && (st[top] != 0) );
        } else if ( strcmp(pos->s, "||") == 0 ) {
          st[top-1] = ( (st[top-1] != 0) || (st[top] != 0) );
        } else {
          fprintf(stderr, "Error: uknown operator [%s]\n", pos->s);
          exit(1);
        }
      }
    } else if ( pos->type == FUNCTION ) {
      for ( i = 0; funcmap[i].f != NULL; i++ ) {
        if ( strcmp(funcmap[i].s, pos->s) == 0 ) {
          if ( top < funcmap[i].narg ) {
            fprintf(stderr, "Error: insufficient arguments for function [%s], has %d, require %d\n",
                    funcmap[i].s, top, funcmap[i].narg);
            exit(1);
          }
          if ( funcmap[i].narg == 0 ) {
            st[top++] = (*funcmap[i].f)();
          } else if ( funcmap[i].narg == 1 ) {
            st[top-1] = (*funcmap[i].f)(st[top-1]);
          } else if ( funcmap[i].narg == 2 ) {
            --top;
            st[top-1] = (*funcmap[i].f)(st[top-1], st[top]);
          } else if ( funcmap[i].narg == 3 ) {
            top -= 2;
            st[top-1] = (*funcmap[i].f)(st[top-1], st[top], st[top+1]);
          }
          break;
        }
      }
    }
#if (defined(DEBUG) && (DEBUG >= 2))
    { /* print out the evaluation stack */
      int j;
      char s[VARNAME_MAX];
      fprintf(stderr, "%10s: ", hist_expr_token2str(s, pos));
      for ( j = 0; j < top; j++ ) fprintf(stderr, "%g ", st[j]);
      fprintf(stderr, "\n");
      //getchar();
    }
#endif
  }
  ans = st[0];
  free(st);
  return ans;
}

/** Get the maximal number of required columns in data */
static int hist_expr_get_colmax(token_t *p)
{
  int colmax = 0;

  for ( ; p->type != NULLTYPE; p++ )
    if ( p->type == COLUMN && p->col > colmax )
      colmax = p->col;
  return colmax;
}

/** Read a long line from file.
 * @param[inout] s string of the line, is allocated when s==NULL and n==0 and grown as needed.
 * @param[inout] n current size of s
 * @return s
 */
static char *file_read_long_line(char **s, size_t *n, FILE *fp)
{
  const int bufsz = 1024;
  char *p;
  size_t cnt, sz;

  if ( *s == NULL && *n == 0 ) {
    *n = bufsz;
    if ( (*s = calloc(*n, sizeof(char))) == NULL ) exit(-1);
  }
  p = *s;
  sz = *n;
  while ( 1 ) {
    if ( fgets(p, sz, fp) == NULL ) return NULL;
    cnt = strlen(*s);
    if ( (*s)[cnt-1] == '\n' ) {
      break;
    } else { /* line too long, expand the buffer */
      *n += bufsz;
      if ( (*s = realloc(*s, (*n)*sizeof(char))) == NULL ) exit(-1);
      p = *s + cnt;
      sz = bufsz;
    }
  }
  return *s;
}

/** parse the histogram command to components */
static int hist_parse_command(param_t *pm, FILE *fplog)
{
  char *p, *q;
  double harr[9];
  int i, np, n, noxyz = 0, hdim = 0;
  const char *delims = " \t\r\n,;";

  /* make a copy of the command */
  n = strlen(pm->command);
  if ( (pm->sxyz = calloc(n + 1, 1)) == NULL ) exit(-1);
  strcpy(pm->sxyz, pm->command);

  /* get the histogram string */
  if ( (p = strstr(pm->sxyz, ">>")) != NULL ) {
    *p = '\0';
    pm->sh = p + 2;
  } else if ( (p = strstr(pm->sxyz, "h(")) != NULL ) {
    noxyz = 1;
    pm->sh = p;
  } else {
    fprintf(stderr, "Error: no \'>>\' in [%s]\n", pm->command);
    return EXIT_FAILURE;
  }

  /* get histogram parameters */
  if ( (p = strchr(pm->sh, '(')) == NULL ) {
    fprintf(stderr, "Error: no ( for the histogram string %s\n", pm->sh);
    return EXIT_FAILURE;
  }
  /* remove the terminal ")" */
  if ( (q = strrchr(pm->sh, ')')) != NULL ) *q = '\0';
  /* split the arguments of h(...) into an array */
  p = strtok(p + 1, delims);
  for ( np = 0; np < 9; ) {
    token_t *t = hist_expr_parse2postfix(p);
    harr[np++] = hist_expr_eval_postfix(t, NULL, 0);
    free(t);
    if ( (p = strtok(NULL, delims)) == NULL ) break;
  }
  hdim = np / 3; /* dimension suggested by the h() function */
  /* assign the parameters of each dimension */
  if ( np >= 3 ) {
    pm->xn = (int) (harr[0] + 0.5);
    pm->xmin = harr[1];
    pm->xmax = harr[2];
    pm->dx = ( pm->xmax - pm->xmin ) / pm->xn;
    pm->zn = pm->yn = pm->xn;
    pm->zmin = pm->ymin = pm->xmin;
    pm->zmax = pm->ymax = pm->xmax;
    pm->dz = pm->dy = pm->dx;
  }
  if ( np >= 6 ) {
    pm->yn = (int) (harr[3] + 0.5);
    pm->ymin = harr[4];
    pm->ymax = harr[5];
    pm->dy = ( pm->ymax - pm->ymin ) / pm->yn;
    pm->zn = pm->yn;
    pm->zmin = pm->ymin;
    pm->zmax = pm->ymax;
    pm->dz = pm->dy;
  }
  if ( np >= 9 ) {
    pm->zn = (int) (harr[6] + 0.5);
    pm->zmin = harr[7];
    pm->zmax = harr[8];
  }

  if ( noxyz ) { /* default: plot the first column */
    pm->dim = ( hdim > 1 ) ? hdim : 1;
    pm->colmax = pm->dim;
    pm->sx = "$1";
    pm->sy = "$2";
    pm->sz = "$3";
    pm->sw = NULL;
  } else {
    /* get the weight expression */
    if ( (p = strstr(pm->sxyz, "::")) != NULL ) {
      *p = '\0';
      pm->sw = p + 2;
    } else {
#ifdef DEBUG
      fprintf(stderr, "Note: no \'::\' for weight in [%s], assuming 1.0\n", pm->sxyz);
#endif
    }

    /* get the x expression */
    pm->dim = 1;
    pm->sx = pm->sxyz;
    /* get the y expression */
    if ( (p = strchr(pm->sx, ':')) != NULL ) {
      pm->dim++;
      *p = '\0';
      pm->sy = p + 1;
      /* get the z expression */
      if ( (p = strchr(pm->sy, ':')) != NULL ) {
        pm->dim++;
        *p = '\0';
        pm->sz = p + 1;
      }
    }
  }

  /* get the postfix expressions for x, y, z */
  pm->colmax = 0;
  if ( pm->sw ) {
    pm->pw = hist_expr_parse2postfix(pm->sw);
    if ( (i = hist_expr_get_colmax(pm->pw)) > pm->colmax ) pm->colmax = i;
  }
  pm->px = hist_expr_parse2postfix(pm->sx);
  if ( (i = hist_expr_get_colmax(pm->px)) > pm->colmax ) pm->colmax = i;
  if ( pm->dim >= 2 ) {
    pm->py = hist_expr_parse2postfix(pm->sy);
    if ( (i = hist_expr_get_colmax(pm->py)) > pm->colmax ) pm->colmax = i;
    if ( pm->dim >= 3 ) {
      pm->pz = hist_expr_parse2postfix(pm->sz);
      if ( (i = hist_expr_get_colmax(pm->pz)) > pm->colmax ) pm->colmax = i;
    }
  }

  if ( fplog != NULL ) {
    if ( pm->dim == 1 ) {
      fprintf(fplog, "# hist1D: x: %s ; weight: %s ; xbins (%d, %g, %g)\n",
              pm->sx, (pm->sw ? pm->sw : "1"), pm->xn, pm->xmin, pm->xmax);
    } else if ( pm->dim == 2 ) {
      fprintf(fplog, "# hist2D: x: %s ; y: %s ; weight: %s ; xbins (%d, %g, %g) ; ybins (%d, %g, %g)\n",
              pm->sx, pm->sy, (pm->sw ? pm->sw : "1"), pm->xn, pm->xmin, pm->xmax, pm->yn, pm->ymin, pm->ymax);
    } else if ( pm->dim == 3 ) {
      fprintf(fplog, "# hist3D: x: %s ; y: %s ; z: %s ; weight: %s ; xbins (%d, %g, %g) ; ybins (%d, %g, %g) ; zbins (%d, %g, %g)\n",
              pm->sx, pm->sy, pm->sz, (pm->sw ? pm->sw : "1"), pm->xn, pm->xmin, pm->xmax,
              pm->yn, pm->ymin, pm->ymax, pm->zn, pm->zmin, pm->zmax);
    }
  }

  return EXIT_SUCCESS;
}

static int hist_from_file(param_t *pm)
{
  int i, col, ix, iy, iz, hn;
  double x, y, z, w, wtot = 0, norm, *hist = NULL;
  FILE *fp;
  double *data = NULL;
  size_t bufsz = 0, lineid = 0;
  char *buf = NULL, *p;
  const char *delims = " \t\r\n,;";

  /* allocate the data array */
  if ( (data = calloc(pm->colmax + 2, sizeof(*data))) == NULL ) exit(-1);

  hn = pm->xn;
  if ( pm->dim >= 2 ) hn *= pm->yn;
  if ( pm->dim >= 3 ) hn *= pm->zn;
  if ( (hist = calloc(hn, sizeof(*hist))) == NULL) exit(-1);
  for ( i = 0; i < hn; i++ ) hist[i] = 0;

  if ( (fp = fopen(pm->dfname, "r")) == NULL ) {
    fprintf(stderr, "Error: cannot read %s\n", pm->dfname);
    goto END;
  }

  /* read file line by line */
  while ( file_read_long_line(&buf, &bufsz, fp) ) {
    if ( buf[0] == '#' || buf[0] == '!' ) continue;
    p = strtok(buf, delims);
    data[0] = ++lineid;
    for ( col = 1; col <= pm->colmax && p != NULL; col++ ) {
      data[col] = atof(p);
      p = strtok(NULL, delims);
      //fprintf(stderr, "%g ", data[col]);
    }

    /* compute the index */
    x = hist_expr_eval_postfix(pm->px, data, col);
    ix = ( x >= pm->xmin ) ? (int) (( x - pm->xmin ) / pm->dx) : -1;
    i = ix;
    if ( ix < pm->xn ) {
      if ( pm->dim >= 2 ) {
        i *= pm->yn;
        y = hist_expr_eval_postfix(pm->py, data, col);
        iy = ( y >= pm->ymin ) ? (int) (( y - pm->ymin ) / pm->dy) : -1;
        if ( iy < pm->yn ) {
          i += iy;
          if ( pm->dim >= 3 ) {
            i *= pm->zn;
            z = hist_expr_eval_postfix(pm->pz, data, col);
            iz = ( z >= pm->zmin ) ? (int) (( z - pm->zmin ) / pm->dz) : -1;
            if ( iz < pm->zn ) {
              i += iz;
            } else {
              i = -1;
            }
          }
        } else {
          i = -1;
        }
      }
    } else {
      i = -1;
    }
    if ( i >= 0 ) {
      w = ( pm->pw != NULL ) ? hist_expr_eval_postfix(pm->pw, data, col) : 1.0;
      //fprintf(stderr, " |  %d %g\n", i, w);
      hist[i] += w;
      wtot += w;
    }
  }

  /* print the histogram */
  w = 0.0; if ( pm->bincenter ) w = 0.5;
  if ( pm->normalize == 0 ) {
    norm = 1.0;
  } else if ( pm->normalize == 1 ) {
    norm = 1.0/wtot;
  } else {
    double dvol = pm->dx;
    if ( pm->dim >= 2 ) dvol *= pm->dy;
    if ( pm->dim >= 3 ) dvol *= pm->dz;
    norm = 1.0/(wtot*dvol);
  }
  if ( pm->dim == 1 ) { /* 1D histogram */
    for ( i = 0 ; i < pm->xn; i++ ) {
      printf("%g %g\n", pm->xmin + (i + w) * pm->dx, hist[i]*norm);
    }
  } else if ( pm->dim == 2 ) { /* 2D histogram */
    for ( i = 0, ix = 0; ix < pm->xn; ix++ ) {
      for ( iy = 0; iy < pm->yn; iy++, i++ ) {
        printf("%g %g %g\n", pm->xmin + (ix + w) * pm->dx,
               pm->ymin + (iy + w) * pm->dy, hist[i]*norm);
      }
      printf("\n");
    }
  } else if ( pm->dim == 3 ) { /* 3D histogram */
    for ( i = 0, ix = 0; ix < pm->xn; ix++ ) {
      for ( iy = 0; iy < pm->yn; iy++, i++ ) {
        for ( iz = 0; iz < pm->zn; iz++, i++ ) {
          printf("%g %g %g %g\n", pm->xmin + (ix + w) * pm->dx,
                 pm->ymin + (iy + w) * pm->dy, pm->zmin + (iz + w) * pm->dz,
                 hist[i]*norm);
        }
      }
      printf("\n");
    }
  }

  fclose(fp);
END:
  free(hist);
  free(buf);
  free(data);
  return 0;
}

#if 0
/** Generate a random test file, if not exist */
static void gen_rand_file(const char *fn) {
  FILE *fp;
  int i, j;
  double x, y;
  if ( (fp = fopen(fn, "r")) != NULL ) {
    fclose(fp);
    return;
  }
  // srand(clock());
  if ((fp = fopen(fn, "w")) != NULL) {
    for ( i = 0; i < 10000; i++ ) {
      for ( x = 0, j = 0; j < 10; j++ )
        x += 1.0*rand()/RAND_MAX;
      x /= 10;
      for ( y = 0, j = 0; j < 10; j++ )
        y += 1.0*rand()/RAND_MAX;
      y /= 10;
      fprintf(fp, "%g %g %g\n", x*2.19, x*y*6.19, y*y*9.19);
    }
    fclose(fp);
  }
}
#endif

int main(int argc, char **argv)
{
  param_t pm;

  if ( param_parse(&pm, argc, argv) != EXIT_SUCCESS )
    return EXIT_FAILURE;
  /* get parameters from the histogram command */
  if ( hist_parse_command(&pm, stdout) != EXIT_SUCCESS )
    return EXIT_FAILURE;
  hist_from_file(&pm);
  param_free(&pm);

  return EXIT_SUCCESS;
}

/* Local Variables:  */
/* c-basic-offset: 2 */
/* End:              */
