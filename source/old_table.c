/*------------------------------------------------------*
 | Author: Maurizio Loreti, aka MLO or (HAM) I3NOO      |
 | Work:   University of Padova - Department of Physics |
 |         Via F. Marzolo, 8 - 35131 PADOVA - Italy     |
 | Phone:  (39)(49) 827-7216     FAX: (39)(49) 827-7102 |
 | EMail:  loreti@padova.infn.it                        |
 | WWW:    http://wwwcdf.pd.infn.it/~loreti/mlo.html    |
 *------------------------------------------------------*

 $Id: old_table.c,v 1.1 2005/03/01 10:06:08 loreti Exp $

 Builds up the LaTeX code to print out the statistical tables in the
 appendix F of "The Book".  These statistical tables are built using
 portable code taken, with modifications, from: Press, Flannery,
 Teukolsky & Vetterling - "Numerical Recipes in C" - 1995 - Cambridge
 University Press - ISBN 0-521-43108-5 .

 ------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define N_T      10        /* Nr. of percentiles in the Student's tables  */
#define STUD_ACC 0.5e-6    /* Accuracy for inversion in the above tables  */
#define N_C      13        /* Nr. of percentiles in the chi square tables */
#define N_NC     30        /* Nr. of ndf entries in the above tables      */
#define N_F       2        /* Nr. of percentiles in the Fisher's tables   */
#define N_FM     12        /* Nr. of M values in the Fisher's tables      */
#define N_FN     30        /* Nr. of N values in the Fisher's tables      */
#define CHI_ACC  0.5e-6    /* Accuracy for inversion in the above tables  */
#define MAXIT    100       /* Iterations (gamma functions)                */
#define ITMAX    100       /* Iterations (beta functions)                 */
#define EPS      3.0E-7    /* Accuracy for the gamma/beta functions       */
#define FPMIN    1.0E-30   /* Underflow guard                             */
#define ZRMAXIT  60        /* Maximum number of iterations in zriddr      */
#define UNUSED   -1.11e30  /* A big negative number                       */

#define SIGN(a,b) ((b) > 0.0 ? fabs(a) : -fabs(a))

#define erroutl(a) errout((a), __FILE__, __LINE__)

/**
 | Local global variables
**/

static int ndf;             /* Number of degrees of freedom */
static int ndg;             /* Second n.d.f. (Fisher)       */
static double target;       /* Probability level to achieve */

static void errout(const char *, const char *, int);

static double istud1(double);
static double pchi1(double);
static double pF1(double);
static void   t_chi(void);
static void   t_gauss(void);
static void   t_stud(void);
static void   t_fisher(void);

static double istud(double, int);
static double pchi(double, int);
static double pF(double, int, int);
static double betacf(double, double, double);
static double betai(double, double, double);
static double gammp(double, double);
static double gammln(double);
static void   gcf(double *, double, double, double *);
static void   gser(double *, double, double, double *);

static double zriddr(double (*func)(double), double, double, double);

int main()
{
  time_t now = time(0);

  puts("% Computer generated file - do not edit!");
  printf("%% Generation date: %s", asctime(localtime(&now)));
  puts("%");
  puts("");

  t_gauss();
  t_stud();
  puts("\\begin{landscape}");
  t_chi();
  t_fisher();
  puts("\\end{landscape}");
  puts("");
  puts("%");
  puts("% End-of-file (tabout.tex)");
  puts("\\endinput");
  return 0;
}

static void errout(
  const char *message,
  const char *file,
  int         line
) {
  fprintf(stderr, "Fatal error - %s\n", message);
  fprintf(stderr, "Error message from file %s, line %d.\n",
          file, line);
  exit(EXIT_FAILURE);
}

static double istud1(
  double x
) {

  /**
   | Auxiliary function passed to "zriddr" in order to find the
   | abscissa "x" where the Student's cumulative distribution function
   | reaches the value "target"; "zriddr" needs a function of a double
   | returning double, and finds its zeros.
  **/

  return istud(x, ndf) - target;
}

static double pchi1(
  double x
) {

  /**
   | Auxiliary function for "zriddr" (chi square).
  **/

  return pchi(x, ndf) - target;
}

static double pF1(
  double x
) {

  /**
   | Auxiliary function for "zriddr" (Fisher).
  **/

  return pF(x, ndf, ndg) - target;
}

static void t_chi(void)
{

  /**
   | Percentile tables of the chi square function; for the probability
   | levels stored in the array c[N_C], finds the zero of the function
   | (integral from zero to x of the chi square function) minus
   | (probability level).
  **/

  double c[N_C] = {
    0.9990, 0.9950, 0.9900, 0.9500, 0.9000,
    0.7500, 0.5000, 0.2500, 0.1000, 0.0500,
    0.0100, 0.0050, 0.0010
  };
  int n[N_NC] = {
     1,  2,  3,  4,  5,  6,  7,  8,  9, 10,
    11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
    21, 22, 23, 24, 25, 26, 27, 28, 29, 30
  };
  double res[N_C][N_NC], xmax;
  int i, j;

  puts("  \\begin{center}");
  puts("    \\footnotesize");
  puts("    \\renewcommand{\\arraystretch}{0.89}");
  puts("    \\begin{tabular}{rrrrrrrrrrrrrr}");
  puts("      & \\multicolumn{13}{c}{Probabilit\\`a (in percentuale)}"
       " \\\\[\\belowrulesep]");
  puts("      \\cmidrule{2-14}");
  puts("      $N$");
  for (i=0; i<N_C; i++) {
    printf("        & \\multicolumn{1}{c}{%.1f}\n", c[i]*100.0);
  }
  puts("        \\\\");
  puts("      \\midrule");

  xmax = 150.0;
  for (j = N_NC;  j-- > 0; xmax = res[0][j]) {
    double xmin;

    ndf = n[j];
    xmin = 0.0;
    for (i = 0;  i < N_C;  i++) {
      target = c[i];
      res[i][j] = zriddr(pchi1, xmin, xmax, CHI_ACC);
      xmax = res[i][j];
    }
  }
  for (j = 0;  j < N_NC;  ) {
    printf("      %d", n[j]);
    for (i = 0;  i < N_C;  i++) {
      printf(" & %.2f", res[i][j]);
      if (i == 5) printf("\n       ");
    }
    puts(" \\\\");
    if (++j % 10 == 0   &&   j != N_NC) {
      puts("      \\cmidrule(lr){2-14}");
    }
  }
  puts("      \\bottomrule");
  puts("      \\addlinespace");
  puts("      \\multicolumn{14}{c}"
       "{Percentili della distribuzione del $\\chi^2$}");
  puts("    \\end{tabular}");
  puts("  \\end{center}");
}

static void t_gauss(void)
{

  /**
   | Tables of the Gauss function: ordinates y = f(x), integral from
   | -x to +x, integral from -infinity to x; all of that, for x going
   | from 0 to 4 with step 0.01
  **/

  int j;
  char *roman[] = {
    "I", "II", "III", "IV", "V"
  };

  for (j = 0;  j < 5;  j++) {
    int i, imin, imax;
    imin = 80 * j;
    imax = imin + 40;

    puts("\\begin{center}");
    puts("  \\footnotesize");
    puts("  \\begin{tabular}{cc}");
    puts("    \\hbox{");
    puts("      \\begin{tabular}{cccc}");
    puts("        \\cmidrule[\\heavyrulewidth]{2-4}");
    puts("        $x$ & $y$ & $I_1$ & $I_2$ \\\\");
    puts("        \\midrule");

    for (i = imin;  i < imax;  i++) {
      double x, y, xs2, i1, i2;

      if (i > imin   &&   i%10 == 0)
        puts("        \\cmidrule(lr){2-4}");

      x = i / 100.0;
      xs2 = x / 1.414213562;
      y = 0.3989422803 * exp(-0.5 * x * x);
      i1 = erf(xs2);
      i2 = 0.5 + 0.5 * erf(xs2);
      printf("        %4.2f & %7.5f & %7.5f & %7.5f \\\\\n",
             x, y, i1, i2);
    }

    puts("        \\bottomrule");
    puts("      \\end{tabular}");
    puts("    } & \\hbox {");
    puts("      \\begin{tabular}{cccc}");
    puts("        \\cmidrule[\\heavyrulewidth]{2-4}");
    puts("        $x$ & $y$ & $I_1$ & $I_2$ \\\\");
    puts("        \\midrule");

    imin = imax;
    imax = imin + 40;
    for (i = imin;  i < imax;  i++) {
      double x, y, xs2, i1, i2;

      if (i > imin   &&   i%10 == 0)
        puts("        \\cmidrule(lr){2-4}");

      x = i / 100.0;
      xs2 = x / 1.414213562;
      y = 0.3989422803 * exp(-0.5 * x * x);
      i1 = erf(xs2);
      i2 = 0.5 + 0.5 * erf(xs2);
      printf("        %4.2f & %7.5f & %7.5f & %7.5f \\\\\n",
             x, y, i1, i2);
    }

    puts("        \\bottomrule");
    puts("      \\end{tabular}");
    puts("    } \\\\");
    puts("    \\addlinespace");
    printf("    \\multicolumn{2}{c}{Tabelle della distribuzione "
           "normale (");
    printf("%s)}\n", roman[j]);
    puts("  \\end{tabular}");
    puts("\\end{center}");
    puts("\\clearpage");
  }
}

static void t_stud(void)
{

  /**
   | Percentile tables of the Student's function; for the probability
   | levels in the array t[N_T], finds the zero of the function
   | (integral from -infinity to x of the Student's function) minus
   | (probability level).
  **/

  double t[N_T] = {
    0.9990, 0.9980, 0.9950, 0.9900, 0.9800,
    0.9500, 0.9000, 0.8000, 0.7500, 0.6000
  };
  double res[N_T] = {320.0};
  int i;

  puts("\\begin{center}");
  puts("  \\footnotesize");
  puts("  \\addtolength{\\tabcolsep}{-1.6pt}");
  puts("  \\begin{tabular}{rrrrrrrrrrr}");
  puts("    & \\multicolumn{10}{c}{Probabilit\\`a (in percentuale)} "
       "\\\\[\\belowrulesep]");
  puts("    \\cmidrule(l){2-11}");
  puts("    $N$");
  for (i=0; i<N_T; i++) {
    printf("      & \\multicolumn{1}{c}{%4.1f}\n", t[i]*100.0);
  }
  puts("      \\\\");
  puts("    \\midrule");

  for (ndf = 1;  ndf < 41;  ndf++) {
    double xmin, xmax;

    if (ndf > 1   &&   ndf%10 == 1) puts("    \\cmidrule(lr){2-11}");

    xmin = 0.0;
    xmax = res[0];
    printf("    %d", ndf);
    for (i = 0;  i < N_T;  i++) {
      target = t[i];
      res[i] = zriddr(istud1, xmin, xmax, STUD_ACC);
      xmax = res[i];
      printf(" & %.3f", xmax);
      if (i == 5) printf("\n     ");
    }
    puts(" \\\\");
  }
  puts("    \\bottomrule");
  puts("    \\addlinespace");
  puts("    \\multicolumn{11}{c}{Percentili della distribuzione "
       "di Student}");
  puts("  \\end{tabular}");
  puts("\\end{center}");
  puts("\\clearpage");
}

static void t_fisher(void)
{

  /**
   | Percentile tables of the Fisher function; for the probability
   | levels in the array f[N_F], finds the zero of the function
   | (integral from 0 to x of the Fisher's function) minus
   | (probability level).
  **/

  double f[N_F] = {0.95, 0.99};
  int m[N_FM] = { 1,  2,  3,  4,  5,  6,  8, 12, 24, 36,
                 48,  1000000};
  int n[N_FN] = { 1,  2,  3,  4,  5,  6,  7,  8,  9, 10,
                 11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
                 21, 22, 23, 24, 25, 30, 40, 60,120, 1000000};
  int kf, km, kn;
  double x, xmax;

  for (kf = 0;  kf < N_F;  kf++) {
    puts("  \\clearpage");
    puts("  \\begin{center}");
    puts("    \\footnotesize");
    puts("    \\renewcommand{\\arraystretch}{0.9}");
    if (kf > 0) {
      puts(   "    \\addtolength{\\tabcolsep}{-0.12pt}");
    }
    puts("    \\begin{tabular}{rrrrrrrrrrrrr}");
    puts("      & \\multicolumn{12}{c}{$M$} \\\\");
    puts("      \\cmidrule(lr){2-13}");
    puts("      \\multicolumn{1}{c}{$N$}");

    for (km = 0;   km < N_FM - 1;   km++) {
      printf("        & \\multicolumn{1}{c}{%d}\n", m[km]);
    }
    puts("        & \\multicolumn{1}{c}{$\\infty$} \\\\");
    puts("      \\midrule");

    xmax = 7000.0;
    target = 1.0 - f[kf];
    for (kn = 0;  kn < N_FN;  ) {
      double new_xmax = 0.0;

      ndg = n[kn];
      if (kn < N_FN - 1) {
        printf("      %d", ndg);
      } else {
        printf("      $\\infty$");
      }

      for (km = 0;  km < N_FM;  ) {
        ndf = m[km];
        x = zriddr(pF1, 0.0, xmax, CHI_ACC);
        printf(" & %.2f", x);
        if (x > new_xmax) new_xmax = x;
        if ((++km % 5) == 0) printf("\n       ");
      }
      puts(" \\\\");
      if ((++kn % 10) == 0) {
        if (kn < N_FN) {
          puts("      \\cmidrule(lr){2-13}");
        } else {
          puts("      \\bottomrule");
        }
      }
      xmax = new_xmax;
    }
    puts("      \\addlinespace");
    puts("      \\multicolumn{13}{c}%");
    printf("        {Percentili della distribuzione di Fisher "
           "per $P=%.2f$}\n", f[kf]);
    puts("    \\end{tabular}");
    puts("  \\end{center}");
  }
}

/*--------------------------------------------------------------------------*

 STATISTICAL FUNCTIONS
 =====================

 Gauss Distribution related functions: erf (defined in <math.h> and
 implemented in the Standard Library, as of the C-99 Standard (ISO/IEC
 9899-1999).

 Defining f(x) = 2 / sqrt(pi) * exp(-x**2) :
 - erf(x) is the integral from 0 to x of f(x): erf(-x) = -erf(x),
   erf(0) = 0, erf(+infinity) = 1.

 The Gauss normalized function is defined as:
 N(x) = 1 / (sqrt(2 * pi)) * exp(-x**2 / 2) ;
 - the integral of N(x) from -x to x is erf(x / sqrt(2)), with x > 0;
 - the integral from -infinity to x is 0.5 + 0.5 * erf(x / sqrt(2)) .

 Gamma related functions:
 - gammp(a,x): Incomplete Gamma Function P(a,x)
 - betai(a,b,x): Incomplete Beta Function Ix(a,b)
 - betacf(a,b,x): continue fraction expansion for Incomplete Beta
   Function
 - gser(s,a,b,g) evaluates the Incomplete Gamma Function P(a,x) by its
   series representation in *s; Returns also Gamma(a) as *g
 - gcf(s,a,b,g) evaluates the Incomplete Gamma Function Q(a,x) by its
   continued fraction representation in *s and Gamma(a) in *g
 - gammln(x): the value of log(Gamma(x)), for x>0

 Student's function:
 - istud(x,n) is the integral from -infinity to x of the Student's
   function with n degrees of freedom

 Chi square function:
 - pchi(x,n) is the integral from 0 to x of the chi square
   distribution function with n degrees of freedom

 Fisher's F function:
 - pF(x,n1,n2) is the integral from 0 to x of the F function with
   n1 and n2 degrees of freedom

 *--------------------------------------------------------------------------*/

static double istud(
  double x,
  int    n
) {
  double y = n / (n + x * x);
  return 1 - 0.5 * betai(0.5 * n,  0.5,  y);
}

static double pchi(
  double x,
  int    n
) {
  return gammp((n / 2.0),  (x / 2.0));
}

static double pF(
  double x,
  int    n1,
  int    n2
) {
  return betai((n2 / 2.0),  (n1 / 2.0),  (n2 / (n2 + n1 * x)));
}

static double betacf(
  double a,
  double b,
  double x
) {
  int m, m2;
  double aa, c, d, del, h, qab, qam, qap;

  qab = a + b;
  qap = a + 1.0;
  qam = a - 1.0;
  c = 1.0;
  d = 1.0 - qab * x / qap;
  if (fabs(d) < FPMIN) d = FPMIN;
  d = 1.0 / d;
  h = d;
  for (m = 1;  m <= MAXIT;  m++) {
    m2 = m + m;
    aa = m * (b - m) * x / ((qam + m2) * (a + m2));
    d = 1.0 + aa * d;
    if (fabs(d) < FPMIN) d = FPMIN;
    c = 1.0 + aa / c;
    if (fabs(c) < FPMIN) c = FPMIN;
    d = 1.0 / d;
    h *= (d * c);
    aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
    d = 1.0 + aa * d;
    if (fabs(d) < FPMIN) d = FPMIN;
    c = 1.0 + aa / c;
    if (fabs(c) < FPMIN) c = FPMIN;
    d = 1.0 / d;
    del = d * c;
    h *= del;
    if (fabs(del - 1.0) < EPS) break;
  }
  if (m > MAXIT) {
    erroutl("a or b too big, or MAXIT too small in BETACF");
  }
  return h;
}

static double betai(
  double a,
  double b,
  double x
) {
  double bt;

  if (x < 0.0   ||   x > 1.0) erroutl("Bad x in routine BETAI");
  if (x == 0.0   ||   x == 1.0) {
    bt = 0.0;
  } else {
    bt = exp(gammln(a + b) - gammln(a) -gammln(b) +
         a * log(x) + b * log(1.0 - x));
  }
  if (x   <   (a + 1.0) / (a + b + 2.0)) {
    return (bt * betacf(a, b, x) / a);
  } else {
    return (1.0 - bt * betacf(b, a, 1.0 - x) / b);
  }
}

static double gammp(
  double a,
  double x
) {
  double gamser, gammcf, gln;

  if (x < 0.0   ||   a <= 0.0) {
    erroutl("Invalid arguments in routine GAMMP");
  }
  if (x < (a + 1.0)) {
    gser(&gamser, a, x, &gln);
    return gamser;
  } else {
    gcf(&gammcf, a, x, &gln);
    return (1.0 - gammcf);
  }
}

static double gammln(
  double xx
) {
  double x, y, tmp, ser;
  int j;
  static double cof[6] = {
    76.18009172947146, -86.50532032941677, 24.01409824083091,
    -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5
  };

  y = x = xx;
  tmp = x + 5.5;
  tmp -= (x + 0.5) * log(tmp);
  ser = 1.000000000190015;
  for (j = 0;  j < 6;  j++) {
    y += 1.0;
    ser += cof[j] / y;
  }
  return -tmp + log(2.5066282746310005 * ser / x);
}

static void gcf(
  double *gammcf,
  double  a,
  double  x,
  double *gln
) {
  int i;
  double an, b, c, d, del, h;

  *gln = gammln(a);
  b = x + 1.0 - a;
  c = 1.0 / FPMIN;
  h = d = 1.0 / b;
  for (i = 1;  i <= ITMAX;  i++) {
    an = -i * (i - a);
    b += 2.0;
    d = an * d + b;
    if (fabs(d) < FPMIN) d = FPMIN;
    c = b + an / c;
    if (fabs(c) < FPMIN) c = FPMIN;
    d = 1.0 / d;
    del = d * c;
    h *= del;
    if (fabs(del - 1.0) < EPS) break;
  }
  if (i > ITMAX ) {
    erroutl("a too large, ITMAX too small in routine GCF");
  }
  *gammcf = exp(-x + a * log(x) - (*gln)) * h;
}

static void gser(
  double *gamser,
  double  a,
  double  x,
  double *gln
) {
  int n;
  double sum, del, ap;

  *gln = gammln(a);
  if (x <= 0.0) {
    if (x < 0.0) erroutl("x less than 0 in routine GSER");
    *gamser = 0.0;
    return;
  } else {
    ap = a;
    del = sum = 1.0 / a;
    for (n = 1;   n <= ITMAX;   n++) {
      ap += 1.0;
      del *= (x / ap);
      sum += del;
      if (fabs(del) < (fabs(sum) * EPS)) {
        *gamser = sum * exp(-x + a * log(x) - (*gln));
        return;
      }
    }
    erroutl("a too large, ITMAX too small in routine GSER");
    return;
  }
}

/*--------------------------------------------------------------------------*

 Find function zeros with the Ridder's method; taken (and modified)
 from: Press, Flannery, Teukolsky & Vetterling - "Numerical Recipes in
 C" - 1995 - Cambridge University Press - ISBN 0-521-43108-5

 *--------------------------------------------------------------------------*/

static double zriddr(
  double (*func)(double),
  double x1,
  double x2,
  double xacc
) {
  int j;
  double ans, fh, fl, fm, fnew, s, xh, xl, xm, xnew;

  fl = (*func)(x1);
  fh = (*func)(x2);
  if ((fl > 0.0  &&  fh < 0.0)   ||   (fl < 0.0  &&  fh > 0.0)) {
    xl = x1;
    xh = x2;
    ans = UNUSED;
    for (j = 1;  j <= ZRMAXIT;  j++) {
      xm = 0.5 * (xl + xh);
      fm = (*func)(xm);
      s = sqrt(fm * fm - fl * fh);
      if (s == 0.0) return ans;
      xnew = xm + (xm - xl) * ((fl > fh ? 1.0 : -1.0) * fm / s);
      if (fabs(xnew - ans) <= xacc) return ans;
      ans = xnew;
      fnew = (*func)(ans);
      if (fnew == 0.0) return ans;
      if (SIGN(fm, fnew) != fm) {
        xl = xm;
        fl = fm;
        xh = ans;
        fh = fnew;
      } else if (SIGN(fl, fnew) != fl) {
        xh = ans;
        fh = fnew;
      } else if (SIGN(fh, fnew) != fh) {
        xl = ans;
        fl = fnew;
      } else {
        erroutl("can't happen!");
      }
      if (fabs(xh - xl) <= xacc) return ans;
    }
    erroutl("more than ZRMAXIT iterations in ZRIDDR");
  } else {
    if (fl == 0.0) return x1;
    if (fh == 0.0) return x2;
    erroutl("root must be bracketed, in ZRIDDR");
  }
  erroutl("can't happen!");
  return 0.0;
}
