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
 * Simple usage: the command
 *
 * hist '$2' data.dat
 *
 * will generate a 1D histogram using the second column of the data file.
 * If the '$2' is omitted, the first column is used.
 *
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
 * plot "<hist 'sin($1)>>h(...)' data.dat" w step
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <float.h>
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
  char *sx, *sy, *sz, *sw, *sh, *sxyz, *cmdcopy;
  int xn, yn, zn, colmax;
  double xmin, xmax, dx;
  double ymin, ymax, dy;
  double zmin, zmax, dz;
  struct token_s *px, *py, *pz, *pw;
  double x, y, z, w, *data;
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
  pm->sx = "$1";
  pm->sy = "$2";
  pm->sz = "$3";
  pm->sw = pm->sh = pm->sxyz = pm->cmdcopy = NULL;
  pm->xn = pm->yn = pm->zn = 0;
  pm->px = pm->py = pm->pz = pm->pw = NULL;
  pm->data = NULL;
}

/** free memory */
static void param_free(param_t *pm)
{
  free(pm->cmdcopy);
  free(pm->px); free(pm->py); free(pm->pz); free(pm->pw);
  free(pm->data);
}

static void print_usage(const param_t *pm)
{
  fprintf(stderr, "Usage:\n");
  fprintf(stderr, "    hist -c -m 0 \'exprX:exprY:exprZ::exprW>>h(nx,xmin,xmax,ny...)\' dfname\n");
  fprintf(stderr, "         -c print bin center, default is bin lower edge, currently %s.\n", (pm->bincenter ? "on" : "off"));
  fprintf(stderr, "         -m normalization, 0: raw counts, 1: make the sum=1.0, 2: make the integral=1.0, currently %d.\n", pm->normalize);
  fprintf(stderr, "         exprX,Y,Z : expressions for computing X,Y,Z values\n");
  fprintf(stderr, "                     use $1 to address the leftmost data column in file.\n");
  fprintf(stderr, "         exprW : expression for weight, if omitted, weight=1.\n");
  fprintf(stderr, "         If the \'>>h(...)\' part is omitted, the histogram range is guessed\n");
  fprintf(stderr, "         in a preliminary scan, which can be inefficient for large data files.\n");
  fprintf(stderr, "         Up to 3-dimensional histograms are supported.\n");
  fprintf(stderr, "\nFor example,\n");
  fprintf(stderr, "\n    hist dfname\n\n");
  fprintf(stderr, "will output a 1D histogram using the first column of the data file \'dfname\' on an automatic grid.\n");
  fprintf(stderr, "\n    hist \'$2:$3**2\' dfname\n\n");
  fprintf(stderr, "will output a 2D histogram using the second column and the square of the third column.\n");
  fprintf(stderr, "\n    hist \'$2>>h(10,-1,1)\' dfname\n\n");
  fprintf(stderr, "will output a 1D histogram of 10 bins from -1 to 1 using the second column.\n");
}

static int param_parse(param_t *pm, int argc, char **argv)
{
  int i, ich, iarg = 0;
  const char *p;

  param_set_default(pm);
  /* parse switches */
  for ( i = 1; i < argc; i++ ) {
    if ( argv[i][0] == '-' && isalpha(argv[i][1]) ) {
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

  if ( iarg == 0
       || ( pm->normalize < 0 || pm->normalize > 2 ) ) {
    print_usage(pm);
    return EXIT_FAILURE;
  }
  if ( iarg == 1 ) {
    pm->dfname = pm->command;
    pm->command = NULL;
  }
#ifdef DEBUG
  fprintf(stderr, "Input file:    %s\n"
          "Command:       %s\n"
          "Bin-center:    %s\n"
          "Normalization: %d\n",
          pm->dfname, (pm->command ? pm->command : ""),
          (pm->bincenter ? "on" : "off"), pm->normalize);
#endif
  return EXIT_SUCCESS;
}

#define VARNAME_MAX 128

static double dblneg(double a) { return -a; }
static double dbladd(double a, double b) { return a + b; }
static double dblsub(double a, double b) { return a - b; }
static double dblmul(double a, double b) { return a * b; }
static double dbldiv(double a, double b) { return a / b; }
static double dblnot(double a) { return a == 0; }
static double dbllt(double a, double b) { return a < b; }
static double dblle(double a, double b) { return a <= b; }
static double dblgt(double a, double b) { return a > b; }
static double dblge(double a, double b) { return a >= b; }
static double dbleq(double a, double b) { return a == b; }
static double dblneq(double a, double b) { return a != b; }
static double dbland(double a, double b) { return (a != 0) && (b != 0); }
static double dblor(double a, double b)  { return (a != 0) || (b != 0); }
static double dblmax(double a, double b) { return ( a > b ) ? a : b; }
static double dblmin(double a, double b) { return ( a < b ) ? a : b; }
static double iif(double c, double a, double b) { return ( c != 0 ) ? a : b; }

/** Tausworthe random number generator for a number in [0, 1) */
static double rand01(void)
{
  static unsigned long s1, s2, s3;

#define TAUSWORTHE(s, a, b, c, d) (((s & c) << d) & 0xffffffffUL) ^ ((((s << a) & 0xffffffffUL) ^ s) >> b)
#define TAUS3() s1 = TAUSWORTHE(s1, 13, 19, 4294967294UL, 12);  \
  s2 = TAUSWORTHE(s2,  2, 25, 4294967288UL,  4);                \
  s3 = TAUSWORTHE(s3,  3, 11, 4294967280UL, 17)

  if ( s1 == 0 ) {
    s3 = (unsigned long) time(NULL) & 0xffffffffUL;
    s1 = ( 219619919 * s3 + 195984 ) & 0xffffffffUL;
    s2 = ( 219619919 * s1 + 195984 ) & 0xffffffffUL;
    s3 = ( 219619919 * s2 + 195984 ) & 0xffffffffUL;
    TAUS3(); TAUS3(); TAUS3(); TAUS3(); TAUS3(); TAUS3();
  }
  TAUS3();
  return ( s1 ^ s2 ^ s3 ) / 4294967296.0;
}

/** expression token types */
enum { NULLTYPE, NUMBER, OPERATOR, FUNCTION, COLUMN, VARIABLE };

typedef struct {
  const char *s;
  double (*f)();  /**< corresponding function */
  int nops;       /**< number of operands */
  int preced;     /**< precedence */
  int rightassoc; /**< right associative */
} opinfo_t;

/** supported operators: */
opinfo_t hist_expr_oplist[] = {
  { "",      NULL, 0,  0, 0 }, /* root of the operator stack */
  { "(",     NULL, 0, 16, 0 }, { ")",     NULL, 0, 16, 0 },
  { "_",   dblneg, 1, 15, 0 }, /* negative */
  { "!",   dblnot, 1, 15, 0 }, { "not", dblnot, 1, 15, 0 },
  { "^",      pow, 2, 13, 1 }, { "**",     pow, 2, 13, 1 }, /* power */
  { "*",   dblmul, 2, 12, 0 }, { "/",   dbldiv, 2, 12, 0 }, { "%",     fmod, 2, 12, 0 },
  { "+",   dbladd, 2, 11, 0 }, { "-",   dblsub, 2, 11, 0 },
  { "<",    dbllt, 2,  9, 0 }, { ">",    dblgt, 2,  9, 0 },
  { "lt",   dbllt, 2,  9, 0 }, { "gt",   dblgt, 2,  9, 0 },
  { "<=",   dblle, 2,  9, 0 }, { ">=",   dblge, 2,  9, 0 },
  { "le",   dblle, 2,  9, 0 }, { "ge",   dblge, 2,  9, 0 },
  { "==",   dbleq, 2,  8, 0 }, { "eq",   dbleq, 2,  8, 0 },
  { "!=",  dblneq, 2,  8, 0 }, { "<>",  dblneq, 2,  8, 0 }, { "neq", dblneq, 2,  8, 0 },
  { "&&",  dbland, 2,  5, 0 }, { "and", dbland, 2,  5, 0 },
  { "||",   dblor, 2,  4, 0 }, { "or",   dblor, 2,  4, 0 },
  { ",",     NULL, 0,  1, 0 },
  { NULL,    NULL, 0,  0, 0 }
};

typedef struct {
  const char *s;
  double (*f)();
  int narg;
} funcmap_t;

funcmap_t hist_expr_funclist[] = {
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
  {NULL, NULL, 0}
};

typedef struct token_s {
  int type;
  char s[VARNAME_MAX]; /**< operator/function name */
  double val;
  int col;             /**< column id in the data */
  opinfo_t *op;        /**< pointer to operator information */
  funcmap_t *func;     /**< pointer to function information */
} token_t;

/** Get a token from the string and return a pointer after the token.
 * @param[out] t an identified token, modified.
 * @param[in] s expression string, not modified.
 * @return a pointer pointing to the character in string s right after the identified token.
 */
static const char *hist_expr_token_get(token_t *t, const char *s)
{
  char *p;
  int i;
  size_t oplen;
  opinfo_t *op;

  /* skip leading spaces */
  while ( *s && isspace(*s) ) s++;
  if ( *s == '\0' ) return NULL;

  t->s[0] = '\0';
  t->val = 0;
  t->col = 0;

  if ( isdigit(*s) || (*s == '.' && isdigit(s[1])) ) { /* .5 is understood as a number */
    t->type = NUMBER;
    t->val = strtod(s, &p);
  } else if ( *s == '$' ) {
    t->type = COLUMN;
    t->col = (int)strtol(s + 1, &p, 10);
  } else {
    /* scan for operators */
    for ( oplen = 3; oplen > 0; oplen-- ) {
      for ( op = hist_expr_oplist + 1; op->s != NULL; op++ ) {
        if ( strlen(op->s) == oplen && strncmp(op->s, s, oplen) == 0 ) {
          t->type = OPERATOR;
          t->op = op;
          strcpy(t->s, op->s);
          return s + oplen;
        }
      }
    }
    if ( isalpha(*s) ) {
      for ( i = 0; isalnum(s[i]) && i < VARNAME_MAX - 1; i++ )
        t->s[i] = s[i];
      t->s[i] = '\0';
      t->type = VARIABLE;
      return s + i;
    } else {
      fprintf(stderr, "Error: unknown token [%s]\n", s);
      exit(1);
    }
  }
  return (const char *) p;
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

/** a = b */
static void hist_expr_token_copy(token_t *a, const token_t *b)
{
  a->type = b->type;
  strncpy(a->s, b->s, VARNAME_MAX);
  a->val = b->val;
  a->col = b->col;
  a->op = b->op;
  a->func = b->func;
}

/** set token `a' as operator `s' */
static void hist_expr_set_op(token_t *a, const char *s)
{
  opinfo_t *op;
  a->type = OPERATOR;
  for ( op = hist_expr_oplist; op->s != NULL; op++ ) {
    if ( strcmp(op->s, s) == 0 ) {
      strcpy(a->s, s);
      a->op = op;
      break;
    }
  }
}

/** convert a variable to a function */
static void hist_expr_make_func(token_t *t)
{
  funcmap_t *func;

  t->type = FUNCTION;
  for ( func = hist_expr_funclist; func->s != NULL; func++ ) {
    if ( strcmp(func->s, t->s) == 0 ) {
      t->func = func;
      return;
    }
  }
  fprintf(stderr, "Error: unknown function [%s]\n", t->s);
  exit(1);
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
  const char *p;
  int n;

  n = strlen(s);

  /* initialize the output queue */
  if ( (que = calloc(n + 1, sizeof(*que))) == NULL ) exit(-1);
  pos = que;

  /* initialize the operator stack */
  if ( (ost = calloc(n + 1, sizeof(*ost))) == NULL ) exit(-1);
  top = ost;
  /* initialize the operator stack with a null operator */
  hist_expr_set_op(top, "");
  hist_expr_token_copy(&tok[1], top); /* tok[1] is the previous token */

  /* when there are tokens to be read, read a token */
  for ( p = s; (p = hist_expr_token_get(tok, p)) != NULL; ) {
    if ( tok->type == NUMBER || tok->type == COLUMN || tok->type == VARIABLE ) {
      /* if the token is a number, then push it to the output queue */
      hist_expr_token_copy(pos++, tok);
    } else if ( tok->s[0] == '(' ) {
      if ( tok[1].type == VARIABLE ) {
        /* convert the preceding variable to a function */
        hist_expr_make_func(--pos);
        /* move the function onto the operator stack */
        hist_expr_token_copy(++top, pos);
      }
      /* push the "(" onto the operator stack */
      hist_expr_token_copy(++top, tok); /* top = token */
    } else if ( tok->s[0] == ')' ) {
      /* while the operator at the top of the operator stack is not a "(" */
      while ( top->s[0] != '(' ) {
        /* pop operators from the operator stack onto the output queue */
        hist_expr_token_copy(pos++, top--);
        if ( top <= ost ) break;
      }
      if ( top->s[0] == '(' ) {
        top--; /* pop the left bracket from the stack */
        if ( top->type == FUNCTION ) {
          /* pop the function from the operator stack onto the output queue */
          hist_expr_token_copy(pos++, top--);
        }
      } else {
        hist_expr_token_copy(pos++, top--);
      }
    } else if ( tok->type == OPERATOR ) {
      if ( (tok->s[0] == '+' || tok->s[0] == '-') && tok->s[1] == '\0'
           && ( tok[1].type == OPERATOR && tok[1].s[0] != ')' ) ) {
        /* handle a unary operator */
        if ( tok->s[0] == '-' ) {
          hist_expr_set_op(tok, "_"); /* change it to negation */
          hist_expr_token_copy(++top, tok); /* top = token */
        } /* skip the unary '+' */
      } else {
        while ( ((top->op->preced > tok->op->preced)
                 || (top->op->preced == tok->op->preced && !top->op->rightassoc))
                && top->s[0] != '(' ) {
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
      for ( t = que; t < pos && t->type != NULLTYPE; t++ )
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

typedef struct {
  char s[VARNAME_MAX];
  double val;
} varmap_t;

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

/** Evaluate the postfix expression stack. */
static double hist_expr_eval_postfix(const token_t *que, const double *arr, int narr)
{
  int i, n, top, narg;
  double *st, ans, (*f)();
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
    } else if ( pos->type == OPERATOR || pos->type == FUNCTION ) {
      if ( pos->type == OPERATOR ) {
        narg = pos->op->nops;
        f = pos->op->f;
      } else {
        narg = pos->func->narg;
        f = pos->func->f;
      }

      if ( f != NULL ) { /* NOTE: operators like "," have no function */
        if ( top < narg ) {
          fprintf(stderr, "Error: insufficient arguments for [%s], has %d, require %d\n",
                  pos->s, top, narg);
          exit(1);
        }
        top -= narg - 1;
        if ( narg == 0 ) {
          st[top-1] = (*f)();
        } else if ( narg == 1 ) {
          st[top-1] = (*f)(st[top-1]);
        } else if ( narg == 2 ) {
          st[top-1] = (*f)(st[top-1], st[top]);
        } else if ( narg == 3 ) {
          st[top-1] = (*f)(st[top-1], st[top], st[top+1]);
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
  ans = ( top > 0 ) ? st[top-1] : 0.0;
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
static int hist_parse_command(param_t *pm)
{
  char *p, *q;
  double harr[9];
  int i, np, n, hdim = 0;
  const char *delims = " \t\r\n,;";

  pm->sxyz = pm->sh = pm->cmdcopy = NULL;
  if ( pm->command != NULL ) {
    /* make a copy of the command */
    n = strlen(pm->command);
    if ( (pm->cmdcopy = calloc(n + 1, 1)) == NULL ) exit(-1);
    strcpy(pm->cmdcopy, pm->command);

    /* split the command to an xyz string and a histogram string */
    if ( (p = strstr(pm->cmdcopy, ">>")) != NULL ) {
      *p = '\0';
      pm->sxyz = pm->cmdcopy;
      pm->sh = p + 2;
    } else if ( (p = strstr(pm->cmdcopy, "h(")) != NULL ) {
      pm->sh = p;
    } else {
      pm->sxyz = pm->cmdcopy;
    }
  }

  if ( pm->sh != NULL ) { /* get the histogram parameters */
    if ( (p = strchr(pm->sh, '(')) == NULL ) {
      fprintf(stderr, "Error: no ( for the histogram string %s\n", pm->sh);
      return EXIT_FAILURE;
    }
    /* remove the terminal ")" */
    if ( (q = strrchr(pm->sh, ')')) != NULL ) *q = '\0';
    /* split the arguments of h(...) into an array */
    p = strtok(p + 1, delims);
    for ( np = 0; p != NULL && np < 9; ) {
      token_t *t = hist_expr_parse2postfix(p);
      harr[np++] = hist_expr_eval_postfix(t, NULL, 0);
      free(t);
      p = strtok(NULL, delims);
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
      if ( np >= 6 ) {
        pm->yn = (int) (harr[3] + 0.5);
        pm->ymin = harr[4];
        pm->ymax = harr[5];
        pm->dy = ( pm->ymax - pm->ymin ) / pm->yn;
        pm->zn = pm->yn;
        pm->zmin = pm->ymin;
        pm->zmax = pm->ymax;
        pm->dz = pm->dy;
        if ( np >= 9 ) {
          pm->zn = (int) (harr[6] + 0.5);
          pm->zmin = harr[7];
          pm->zmax = harr[8];
        }
      }
    }
  }

  if ( pm->sxyz == NULL ) { /* default: plot the first column */
    pm->dim = ( hdim > 1 ) ? hdim : 1;
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

  /* allocate the data array */
  if ( (pm->data = calloc(pm->colmax + 2, sizeof(*pm->data))) == NULL ) exit(-1);

  return EXIT_SUCCESS;
}

/** get the x, y, z, and weight from the line buffer */
static int hist_get_xyzw(char *buf, param_t *pm, size_t *lineid)
{
  char *p;
  int col;
  const char *delims = " \t\r\n,;";

  for ( p = buf; *p && isspace(*p); p++ ) ; /* skip leading spaces */
  if ( *p == '#' || *p == '!' ) return 1; /* skip a comment line */

  /* split the line to an array of data */
  p = strtok(p, delims);
  pm->data[0] = ++(*lineid);
  for ( col = 1; col <= pm->colmax && p != NULL; col++ ) {
    pm->data[col] = atof(p);
    p = strtok(NULL, delims);
    //fprintf(stderr, "%g ", pm->data[col]);
  }

  pm->x = hist_expr_eval_postfix(pm->px, pm->data, col);
  if ( pm->dim >= 2 ) {
    pm->y = hist_expr_eval_postfix(pm->py, pm->data, col);
    if ( pm->dim >= 3 ) {
      pm->z = hist_expr_eval_postfix(pm->pz, pm->data, col);
    }
  }
  pm->w = ( pm->pw != NULL ) ? hist_expr_eval_postfix(pm->pw, pm->data, col) : 1.0;
  return 0;
}

/** produce a suitable grid for (*xmin, *xmax)
 * @return the number of bins */
static int hist_make_grid(double *xmin, double *xmax, double *dx)
{
  double del = *xmax - *xmin, xave;
  int n;

  if ( del <= 0 ) {
    xave = (*xmin + *xmax)/2;
    del = fabs(xave);
    if ( del < 1e-6 ) del = 1;
    *xmin = xave - del/2;
    *xmax = xave + del/2;
  }
  *dx = pow(10, floor(log10(del)) - 1);
  n = (int) (del / *dx);
  if ( n > 50 ) {
    *dx *= 5;
  } else if ( n > 20 ) {
    *dx *= 2;
  }
  *xmin = floor(*xmin / *dx) * (*dx);
  *xmax = (floor(*xmax / *dx) + 1) * (*dx);
  return (int) ( ( *xmax - *xmin ) / *dx + 0.5 );
}

/** guess the histogram range */
static int hist_guess_range(param_t *pm)
{
  FILE *fp;
  size_t bufsz = 0, lineid = 0;
  char *buf = NULL;

  if ( (fp = fopen(pm->dfname, "r")) == NULL ) {
    fprintf(stderr, "Error: cannot read %s\n", pm->dfname);
    return -1;
  }

  pm->xmin = pm->ymin = pm->zmin = DBL_MAX;
  pm->xmax = pm->ymax = pm->zmax = -DBL_MAX;
  /* read file line by line */
  while ( file_read_long_line(&buf, &bufsz, fp) ) {
    /* get the x, y, z, and weight from the line */
    if ( hist_get_xyzw(buf, pm, &lineid) != 0 ) continue;
    pm->xmin = dblmin(pm->xmin, pm->x);
    pm->xmax = dblmax(pm->xmax, pm->x);
    if ( pm->dim >= 2 ) {
      pm->ymin = dblmin(pm->ymin, pm->y);
      pm->ymax = dblmax(pm->ymax, pm->y);
      if ( pm->dim >= 3 ) {
        pm->zmin = dblmin(pm->zmin, pm->z);
        pm->zmax = dblmax(pm->zmax, pm->z);
      }
    }
  }
  fclose(fp);

  if ( lineid > 0 ) {
    pm->xn = hist_make_grid(&pm->xmin, &pm->xmax, &pm->dx);
#ifdef DEBUG
    fprintf(stderr, "Auto x range: (%g, %g) bin %g, %d bins\n",
            pm->xmin, pm->xmax, pm->dx, pm->xn);
#endif
    if ( pm->dim >= 2 ) {
      pm->yn = hist_make_grid(&pm->ymin, &pm->ymax, &pm->dy);
#ifdef DEBUG
      fprintf(stderr, "Auto y range: (%g, %g) bin %g, %d bins\n",
              pm->ymin, pm->ymax, pm->dy, pm->yn);
#endif
      if ( pm->dim >= 3 ) {
        pm->zn = hist_make_grid(&pm->zmin, &pm->zmax, &pm->dz);
#ifdef DEBUG
        fprintf(stderr, "Auto z range: (%g, %g) bin %g, %d bins\n",
                pm->zmin, pm->zmax, pm->dz, pm->zn);
#endif
      }
    }
  }

  free(buf);
  return lineid <= 0;
}

static int hist_from_file(param_t *pm)
{
  int i, ix, iy, iz, hn;
  double wtot = 0, del, norm, *hist;
  FILE *fp;
  size_t bufsz = 0, lineid = 0;
  char *buf = NULL;

  hn = pm->xn;
  if ( pm->dim >= 2 ) {
    hn *= pm->yn;
    if ( pm->dim >= 3 ) hn *= pm->zn;
  }
  if ( (hist = calloc(hn, sizeof(*hist))) == NULL) exit(-1);
  for ( i = 0; i < hn; i++ ) hist[i] = 0;

  if ( pm->dim == 1 ) {
    fprintf(stdout, "# hist1D: x: %s ; weight: %s ; xbins (%d, %g, %g)\n",
            pm->sx, (pm->sw ? pm->sw : "1"), pm->xn, pm->xmin, pm->xmax);
  } else if ( pm->dim == 2 ) {
    fprintf(stdout, "# hist2D: x: %s ; y: %s ; weight: %s ; xbins (%d, %g, %g) ; ybins (%d, %g, %g)\n",
            pm->sx, pm->sy, (pm->sw ? pm->sw : "1"), pm->xn, pm->xmin, pm->xmax, pm->yn, pm->ymin, pm->ymax);
  } else if ( pm->dim == 3 ) {
    fprintf(stdout, "# hist3D: x: %s ; y: %s ; z: %s ; weight: %s ; xbins (%d, %g, %g) ; ybins (%d, %g, %g) ; zbins (%d, %g, %g)\n",
            pm->sx, pm->sy, pm->sz, (pm->sw ? pm->sw : "1"), pm->xn, pm->xmin, pm->xmax,
            pm->yn, pm->ymin, pm->ymax, pm->zn, pm->zmin, pm->zmax);
  }

  if ( (fp = fopen(pm->dfname, "r")) == NULL ) {
    fprintf(stderr, "Error: cannot read %s\n", pm->dfname);
    goto END;
  }

  /* read file line by line */
  while ( file_read_long_line(&buf, &bufsz, fp) ) {
    /* get the x, y, z, and weight from the line */
    if ( hist_get_xyzw(buf, pm, &lineid) != 0 ) continue;

    /* compute the histogram index */
    ix = ( pm->x >= pm->xmin ) ? (int) (( pm->x - pm->xmin ) / pm->dx) : -1;
    i = ( ix < pm->xn ) ? ix : -1;
    if ( i >= 0 && pm->dim >= 2 ) {
      iy = ( pm->y >= pm->ymin ) ? (int) (( pm->y - pm->ymin ) / pm->dy) : -1;
      i = ( iy < pm->yn ) ? i * pm->yn + iy : -1;
      if ( i >= 0 && pm->dim >= 3 ) {
        iz = ( pm->z >= pm->zmin ) ? (int) (( pm->z - pm->zmin ) / pm->dz) : -1;
        i = ( iz < pm->zn ) ? i * pm->zn + iz : -1;
      }
    }
    if ( i >= 0 ) {
      //fprintf(stderr, " |  %d %g\n", i, w);
      hist[i] += pm->w;
      wtot += pm->w;
    }
  }

  /* print the histogram */
  del = ( pm->bincenter ) ? 0.5 : 0.0;
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
      printf("%g %g\n", pm->xmin + (i + del) * pm->dx, hist[i]*norm);
    }
  } else if ( pm->dim == 2 ) { /* 2D histogram */
    for ( i = 0, ix = 0; ix < pm->xn; ix++ ) {
      for ( iy = 0; iy < pm->yn; iy++, i++ ) {
        printf("%g %g %g\n", pm->xmin + (ix + del) * pm->dx,
               pm->ymin + (iy + del) * pm->dy, hist[i]*norm);
      }
      printf("\n");
    }
  } else if ( pm->dim == 3 ) { /* 3D histogram */
    for ( i = 0, ix = 0; ix < pm->xn; ix++ ) {
      for ( iy = 0; iy < pm->yn; iy++, i++ ) {
        for ( iz = 0; iz < pm->zn; iz++, i++ ) {
          printf("%g %g %g %g\n", pm->xmin + (ix + del) * pm->dx,
                 pm->ymin + (iy + del) * pm->dy, pm->zmin + (iz + del) * pm->dz,
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
  // srand(time(NULL));
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
  if ( hist_parse_command(&pm) != EXIT_SUCCESS )
    return EXIT_FAILURE;
  if ( pm.xn <= 0 && hist_guess_range(&pm) != 0 )
    return EXIT_FAILURE;
  hist_from_file(&pm);
  param_free(&pm);

  return EXIT_SUCCESS;
}

/* Local Variables:  */
/* c-basic-offset: 2 */
/* End:              */
