/*
 * Copyright (c) 2017
 *
 *     Cheng Zhang and Yuan Mei
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
 *
 * hist 'sin($1+$2^(log($3))):$2::1/$5>>h(100, -1.0, 1.0, 50, 0.0, 2.0)' data.dat
 *            X-axis expr      Y    W      nx  xmin xmax  ny ymin ymax
 *
 * will generate a 2D histogram according to the x,y,z binning
 * specification using data from the file, transformed by the X,Y,Z
 * and W (weight) expression, then filled into the histogram with
 * weight.  Weight is assumed to be 1 if ::W is absent.  1, 2 and 3
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
#include <unistd.h>

/** Parameters settable from commandline */
typedef struct param
{
  int normalize;  /**< normalize histogram instead of raw counts */
  int bincenter;  /**< print bin center instead of lower edge */
} param_t;

param_t param_default = {
  .normalize = 0,
  .bincenter = 0
};

static void print_usage(const param_t *pm)
{
  printf("Usage:\n");
  printf("    hist -c -m 0 'exprX:exprY:exprZ::exprW>>h(nx,xl,xh,ny...)' dfname\n");
  printf("         -c print bin center, default is bin lower edge.\n");
  printf("         -m [0]:counts, 1:normalize to sum=1.0, 2:integral=1.0.\n");
  printf("         exprX,Y,Z : expressions for computing X,Y,Z values\n");
  printf("                     use $1 to address the leftmost data column in file.\n");
  printf("         exprW expression for weight, if omitted, weight=1.\n");
  printf("         Up to 3-dimensional histograms are supported.\n");
}

/** expression token types */
enum { NULLTYPE, NUMBER, OPERATOR, FUNCTION, COLUMN, VARIABLE };

static const char *ops = ",+-*/%^_()";

#define VARNAME_MAX 128

typedef struct {
  int type;
  char s[VARNAME_MAX]; /* operator/function name */
  double val;
  int col;             /* column id in the data */
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
    p = (char*)s + 1;
    if ( *s == '*' && s[1] == '*' ) {
      t->s[0] = '^'; /* convert "**" to "^" */
      p++;
    }
  } else if ( isdigit(*s) ) {
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
    t->col = (int) strtol(s + 1, &p, 10);
  } else {
    fprintf(stderr, "Error: unknown token [%s]\n", s);
    exit(1);
  }
  return p;
}

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

  if ( prtable['+'] == 0 ) { /* initialize the precedence table */
    prtable['('] = prtable[')'] = 16;
    prtable['_'] = 15;
    prtable['^'] = 13;
    prtable['*'] = prtable['/'] = prtable['%'] = 12;
    prtable['+'] = prtable['-'] = 11;
    prtable[','] = 1;
  }

  /* functions take the highest precedence */
  return isalpha(*s) ? 1000 : prtable[(int)(*s)];
}

static int hist_expr_isleftassoc(const char *s)
{
  if ( *s == '^' ) return 0;
  return 1;
}

/** Convert an expression to a stack of postfix expression.
 * Shunting Yard algorithm:
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

  /* where there are tokens to be read, read a token */
  for ( p = (char*)s; (p = hist_expr_token_get(tok, p)) != NULL; ) {
    if ( tok->type == NUMBER || tok->type == COLUMN || tok->type == VARIABLE ) {
      hist_expr_token_copy(pos++, tok);
    } else if ( tok->s[0] == ',' ) {
      /* do nothing for the comma */
    } else if ( tok->s[0] == '(' || tok->type == FUNCTION ) {
      hist_expr_token_copy(++top, tok); /* top = token */
    } else if ( tok->s[0] == ')' ) {
      /* while the operator at the top of the operator stack is not a left bracket
       * and is not a function */
      while ( top->s[0] != '(' && !isalpha(top->s[0]) ) {
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
                && top->s[0] != '(' && !isalpha(top->s[0]) ) {
          /* pop operators from the operator stack onto the output queue */
          hist_expr_token_copy(pos++, top--); /* pos = top */
          if ( top <= ost ) break;
        }
        //printf("pushing tok %s\n", tok->s);
        hist_expr_token_copy(++top, tok); /* top = token */
      }
    }
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

static double max(double a, double b) { return (a > b) ? a : b; }
static double min(double a, double b) { return (a < b) ? a : b; }
static double iif(double c, double a, double b) { return (c == 0) ? a : b; }

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
  {"min", min, 2},
  {"max", max, 2},
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
static double hist_expr_eval_postfix(const token_t *que, const double *arr)
{
  int i, n, top;
  double *st, ans;
  token_t *pos;

  /* determine the length of the expression */
  for ( n = 0; que[n].type != NULLTYPE; n++ ) ;
  /* allocate the evaluation stack */
  if ((st = calloc(n, sizeof(*st))) == NULL) exit(-1);
  top = 0;

  for ( pos = (token_t*)que; pos->type != NULLTYPE; pos++ ) {
    if ( pos->type == NUMBER ) {
      st[top++] = pos->val;
    } else if ( pos->type == COLUMN ) {
      st[top++] = arr[pos->col - 1]; /* array index off-by-one */
    } else if ( pos->type == VARIABLE ) {
      for ( i = 0; varmap[i].s[0] != '\0'; i++ ) {
        if ( strcmp(varmap[i].s, pos->s) == 0 ) {
          st[top++] = varmap[i].val;
          break;
        }
      }
    } else if ( pos->type == OPERATOR ) {
      if ( pos->s[0] == '_' ) { /* unary operators */
        st[top-1] = -st[top-1];
      } else { /* binary operators */
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
        }
      }
    } else if ( pos->type == FUNCTION ) {
      for ( i = 0; funcmap[i].f != NULL; i++ ) {
        if ( strcmp(funcmap[i].s, pos->s) == 0 ) {
          if ( funcmap[i].narg == 1 ) {
            st[top-1] = (*funcmap[i].f)(st[top-1]);
          } else if ( funcmap[i].narg == 2 ) {
            --top;
            st[top-1] = (*funcmap[i].f)(st[top-1], st[top]);
          }
          break;
        }
      }
    }
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

static int hist_from_file(const char *command, const char *dfname, const param_t *pm)
{
  token_t *px, *py = NULL, *pz = NULL, *pw;
  char *sx, *sy = NULL, *sz = NULL, *sw = NULL;
  const char *p, *q;
  int dim = 1, i, n, colmax, ix, iy, iz, xn, yn, zn, hn;
  double x, xmin, xmax, dx;
  double y, ymin, ymax, dy;
  double z, zmin, zmax, dz;
  double w, wtot = 0, norm, *hist = NULL;
  FILE *fp;
  double *data = NULL;
  size_t bufsz = 0;
  char *buf = NULL;
  const char *delims = " \t\r\n,";

  /* get the x expression */
  if ( (p = strpbrk(command, ":>")) == NULL ) {
    fprintf(stderr, "Syntax error: no '>' in [%s]\n", command);
    return -1;
  }
  n = p - command;
  if ( (sx = calloc(n + 1, 1)) == NULL ) exit(-1);
  strncpy(sx, command, n + 1);
  sx[n] = '\0';
  px = hist_expr_parse2postfix(sx);
  colmax = hist_expr_get_colmax(px);
  p++;

  /* get the y expression */
  if ( *p != ':' && *p != '>' && (q = strpbrk(p, ":>")) != NULL ) {
    dim++;
    n = q - p;
    if ( (sy = calloc(n + 1, 1)) == NULL ) exit(-1);
    strncpy(sy, p, n + 1);
    sy[n] = '\0';
    py = hist_expr_parse2postfix(sy);
    if ( (i = hist_expr_get_colmax(py)) > colmax ) colmax = i;
    p = (*q == ':') ? q + 1 : q;

    /* get the z expression */
    if ( *p != ':' && *p != '>' && (q = strpbrk(p, ":>")) != NULL ) {
      dim++;
      n = q - p;
      if ( (sz = calloc(n + 1, 1)) == NULL ) exit(-1);
      strncpy(sz, p, n + 1);
      sz[n] = '\0';
      pz = hist_expr_parse2postfix(sz);
      if ( (i = hist_expr_get_colmax(pz)) > colmax ) colmax = i;
      p = (*q == ':') ? q + 1 : q;
    }
  }

  /* get the weight */
  if ( (p = strstr(command, "::")) != NULL ) {
    p += 2;
    q = strstr(p, ">>");
    n = q - p;
    if ( (sw = calloc(n + 1, 1)) == NULL ) exit(-1);
    strncpy(sw, p, n + 1);
    sw[n] = '\0';
    pw = hist_expr_parse2postfix(sw);
    if ( (i = hist_expr_get_colmax(pw)) > colmax ) colmax = i;
  } else {
    // fprintf(stderr, "no :: for weight in [%s|%s], assuming 1.0\n", command, p);
    sw = NULL;
    pw = NULL;
    p = command;
  }

  /* allocate the data array */
  if ( (data = calloc(colmax + 1, sizeof(*data))) == NULL ) exit(-1);

  /* get histogram parameters */
  if ( (p = strstr(p, "h(")) == NULL ) {
    fprintf(stderr, "Syntax error: no histogram specification h(...) in [%s]\n", command);
    goto END;
  }
  if (3 != sscanf(p, " h ( %d, %lf, %lf%n", &xn, &xmin, &xmax, &i)) {
    fprintf(stderr, "Error: invalid histogram x-axis specification [%s]\n", p);
    goto END;
  }
  p += i;
  dx = (xmax - xmin) / xn;
  zn = yn = xn;
  zmin = ymin = xmin;
  zmax = ymax = xmax;
  dz = dy = dx;

  if ( dim >= 2 ) { /* scan histogram y parameters */
    while ( *p && isspace(*p) ) p++;
    if ( *p != ')' ) {
      p++; /* skip a , or ; */
      if (3 != sscanf(p, "%d, %lf, %lf%n", &yn, &ymin, &ymax, &i)) {
        fprintf(stderr, "Error: invalid histogram y-axis specification [%s]\n", p);
        goto END;
      }
      p += i;
      dy = (ymax - ymin) / yn;
      zn = yn;
      zmin = ymin;
      zmax = ymax;
      dz = dy;
      if ( dim >= 3 ) { /* scan histogram z parameters */
        while ( *p && isspace(*p) ) p++;
        if ( *p != ')' ) {
          p++; /* skip a , or ; */
          if (3 != sscanf(p, "%d, %lf, %lf%n", &zn, &zmin, &zmax, &i)) {
            fprintf(stderr, "Error: invalid histogram z-axis specification [%s]\n", p);
            goto END;
          }
          p += i;
          dz = (zmax - zmin) / zn;
        }
      }
    }
  }

  if ( dim == 1 ) {
    fprintf(stdout, "# hist1D: x: %s, weight: %s, xbins (%d, %g, %g)\n",
        sx, (sw ? sw : "1"), xn, xmin, xmax);
  } else if ( dim == 2 ) {
    fprintf(stdout, "# hist2D: x: %s, y: %s, weight: %s, xbins (%d, %g, %g); ybins (%d, %g, %g)\n",
        sx, sy, (sw ? sw : "1"), xn, xmin, xmax, yn, ymin, ymax);
  } else if ( dim == 3 ) {
    fprintf(stdout, "# hist3D: x: %s, y: %s, z: %s, weight: %s, xbins (%d, %g, %g); ybins (%d, %g, %g); zbins (%d, %g, %g)\n",
        sx, sy, sz, (sw ? sw : "1"), xn, xmin, xmax, yn, ymin, ymax, zn, zmin, zmax);
  }

  hn = xn;
  if ( dim >= 2 ) hn *= yn;
  if ( dim >= 3 ) hn *= zn;
  if ( (hist = calloc(hn, sizeof(*hist))) == NULL) exit(-1);
  for ( i = 0; i < hn; i++ ) hist[i] = 0;

  if ((fp = fopen(dfname, "r")) == NULL ) {
    fprintf(stderr, "Error: cannot read %s\n", dfname);
    goto END;
  }

  /* read file line by line */
  while ( file_read_long_line(&buf, &bufsz, fp) ) {
    if (buf[0] == '#' || buf[0] == '!') continue;
    p = strtok(buf, delims);
    for ( i = 0; i < colmax && p != NULL; i++ ) {
      data[i] = atof(p);
      p = strtok(NULL, delims);
      //fprintf(stderr, "%g ", data[i]);
    }

    /* compute the index */
    x = hist_expr_eval_postfix(px, data);
    ix = ( x >= xmin ) ? ( x - xmin ) / dx : -1;
    i = -1;
    if ( ix < xn ) {
      i = ix;
      if ( dim >= 2 ) {
        i *= yn;
        y = hist_expr_eval_postfix(py, data);
        iy = ( y >= ymin ) ? ( y - ymin ) / dy : -1;
        if ( iy < yn ) {
          i += iy;
          if ( dim >= 3 ) {
            i *= zn;
            z = hist_expr_eval_postfix(pz, data);
            iz = ( z >= zmin ) ? ( z - zmin ) / dz : -1;
            if ( iz < zn ) {
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
      w = (pw != NULL) ? hist_expr_eval_postfix(pw, data) : 1.0;
      //fprintf(stderr, " |  %d %g\n", i, w);
      hist[i] += w;
      wtot += w;
    }
  }

  /* print the histogram */
  w = 0.0; if ( pm->bincenter ) w = 0.5;
  if ( dim == 1 ) { /* 1D histogram */
    if ( pm->normalize == 0 ) {
      norm = 1.0;
    } else if ( pm->normalize == 1 ) {
      norm = 1.0/wtot;
    } else {
      norm = 1.0/(wtot*dx);
    }
    for (i=0 ; i < xn; i++ ) {
      printf("%g %g\n", xmin + (i + w) * dx, hist[i]*norm);
    }
  } else if ( dim == 2 ) { /* 2D histogram */
    if ( pm->normalize == 0 ) {
      norm = 1.0;
    } else if ( pm->normalize == 1 ) {
      norm = 1.0/wtot;
    } else {
      norm = 1.0/(wtot*dx*dy);
    }
    for ( i = 0, ix = 0; ix < xn; ix++ ) {
      for ( iy = 0; iy < yn; iy++, i++ ) {
        printf("%g %g %g\n", xmin + (ix + w) * dx,
            ymin + (iy + w) * dy, hist[i]*norm);
      }
      printf("\n");
    }
  } else if ( dim == 3 ) { /* 3D histogram */
    if ( pm->normalize == 0 ) {
      norm = 1.0;
    } else if ( pm->normalize == 1 ) {
      norm = 1.0/wtot;
    } else {
      norm = 1.0/(wtot*dx*dy*dz);
    }
    for ( i = 0, ix = 0; ix < xn; ix++ ) {
      for ( iy = 0; iy < yn; iy++, i++ ) {
        for ( iz = 0; iz < zn; iz++, i++ ) {
          printf("%g %g %g %g\n", xmin + (ix + w) * dx,
              ymin + (iy + w) * dy, zmin + (iz + w) * dz,
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
  free(sx); free(sy); free(sz); free(sw);
  free(px); free(py); free(pz); free(pw);
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
  int optC = 0;
  param_t pm;
  //const char *command = "sin(-$1+$2)::1.0>>h(100, -1.0, 1.0)";
  const char *command = "sin(-$1+$2):cos($3)::1.0>>h(100, -1.0, 1.0)";
  const char *dfname = "data.dat";

  memcpy(&pm, &param_default, sizeof(pm));
  /* parse switches */
  while((optC = getopt(argc, argv, "cm:")) != -1) {
    switch(optC) {
    case 'c':
      pm.bincenter = 1;
      break;
    case 'm':
      pm.normalize = strtol(optarg, NULL, 0);
      break;
    default:
      print_usage(&pm);
      return EXIT_FAILURE;
      break;
    }
  }
  argc -= optind;
  argv += optind;

  if(argc<2) {
    print_usage(&pm);
    return EXIT_FAILURE;
  }
  command = argv[0];
  dfname  = argv[1];
  hist_from_file(command, dfname, &pm);

  return EXIT_SUCCESS;
}

/* Local Variables:  */
/* c-basic-offset: 2 */
/* End:              */
