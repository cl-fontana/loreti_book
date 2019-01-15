/*------------------------------------------------------*
 | Author: Maurizio Loreti, aka MLO or (HAM) I3NOO      |
 | Work:   University of Padova - Department of Physics |
 |         Via F. Marzolo, 8 - 35131 PADOVA - Italy     |
 | Phone:  +39 (049) 827-7216   FAX: +39 (049) 827-7102 |
 | EMail:  loreti@pd.infn.it                            |
 | WWW:    http://www.pd.infn.it/~loreti/mlo.html       |
 *------------------------------------------------------*

 $Id: table.c,v 1.1 2005/03/01 10:06:08 loreti Exp $

 Builds up the LaTeX code that will typeset the statistical tables in
 the Appendix F of "The Book".  These statistical tables are built
 using the mathematical constants and procedures of the GNU Scientific
 Library (GSL): see http://www.gnu.org/software/gsl/ for more info.

 *------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include <gsl/gsl_cdf.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_sf_erf.h>

#define FSOLVER_TYPE    gsl_root_fsolver_brent
#define FISH_EPS        3.0E-7
#define FISH_ITER       1000

static void t_chi(void);
static void t_fisher(void);
static void t_gauss(void);
static void t_stud(void);

static double fishAux(double x, void *params);

typedef struct sFishParams {    /* Auxiliary structure used for */
  int    ndf;                   /*   the Fisher's distribution: */
  int    ndg;                   /*     degrees of freedom, and  */
  double target;                /*       confidence level.      */
} fishParams;

int main()
{
  time_t now = time(0);

  puts("% Computer generated file - do not edit!");
  printf("%% Generation date: %s", asctime(localtime(&now)));
  puts("%");
  puts("");
  t_gauss();                    /* Gauss function tables */
  t_stud();                     /* Student's function    */
  puts("\\begin{landscape}");
  t_chi();                      /* Chi-square function   */
  t_fisher();                   /* Fisher's distribution */
  puts("\\end{landscape}");
  puts("");
  puts("%");
  puts("% End-of-file (tabout.tex)");
  puts("\\endinput");
  return 0;
}

static void t_gauss(void)
{

  /**
   | Tables of the unit Gauss function (mean 0 and variance 1)
   | =========================================================
   |
   | For every x from 0 to 4 with step 0.01, prints:
   | - the Probability Distribution Function (PDF) ordinates;
   | - the integral from -x to +x of the PDF;
   | - the integral from -infinity to x of the PDF,
   | using the GSL version, gsl_sf_erf(), of the so-called "error
   | function" erf(x).
   |
   | Defining f(x) = 2 / sqrt(pi) * exp(-x**2) :
   | - erf(x) is the integral from 0 to x of f(x);
   | - erf(-x) = -erf(x),   erf(0)=0,   and erf(+infinity) = 1.
   |
   | The unit Gauss function: N(x) = 1 / (sqrt(2*pi)) * exp(-x**2/2) ;
   | - the integral of N(x) from -x to x (x>0) is erf(x/sqrt(2)) ;
   | - the integral from -infinity to x is 0.5 + 0.5 * erf(x/sqrt(2)).
  **/

  int           j;
  const double  factor  = M_SQRT2 * M_SQRTPI;
  const char   *roman[] = { "I", "II", "III", "IV", "V" };
  const int     nPages  = sizeof(roman) / sizeof(roman[0]);

  for (j = 0;   j < nPages;   j++) {
    int i, imin, imax;

    puts("\\begin{center}");
    puts("  \\footnotesize");
    puts("  \\begin{tabular}{cc}");
    puts("    \\hbox{");
    puts("      \\begin{tabular}{cccc}");
    puts("        \\cmidrule[\\heavyrulewidth]{2-4}");
    puts("        $x$ & $y$ & $I_1$ & $I_2$ \\\\");
    puts("        \\midrule");

    /**
     | Left column
    **/

    imin = 80 * j;
    imax = imin + 40;
    for (i = imin;   i < imax;   i++) {
      double x, y, xs2, i1, i2;

      if (i > imin   &&   i % 10 == 0) {
        puts("        \\cmidrule(lr){2-4}");
      }

      x   = i / 100.0;
      xs2 = x / M_SQRT2;
      y   = exp(-0.5 * x * x) / factor;
      i1  = gsl_sf_erf(xs2);
      i2  = 0.5 + 0.5 * i1;

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

    /**
     | Right column
    **/

    imin = imax;
    imax = imin + 40;
    for (i = imin;   i < imax;   i++) {
      double x, y, xs2, i1, i2;

      if (i > imin   &&   i % 10 == 0) {
        puts("        \\cmidrule(lr){2-4}");
      }

      x   = i / 100.0;
      xs2 = x / M_SQRT2;
      y   = exp(-0.5 * x * x) / factor;
      i1  = gsl_sf_erf(xs2);
      i2  = 0.5 + 0.5 * i1;

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
   | Percentile tables of the Student's function
   | ===========================================
   |
   | For the probability levels in the array t[N_T], inverts the
   | cumulative Student's distribution function; repeating all for a
   | number of degrees of freedom from 1 to 40.
  **/

  const double t[] = {
    0.9990, 0.9980, 0.9950, 0.9900, 0.9800,
    0.9500, 0.9000, 0.8000, 0.7500, 0.6000
  };
  const int N_T = sizeof(t) / sizeof(t[0]);

  int i, ndf;

  puts("\\begin{center}");
  puts("  \\footnotesize");
  puts("  \\addtolength{\\tabcolsep}{-1.6pt}");
  puts("  \\begin{tabular}{rrrrrrrrrrr}");
  puts("    & \\multicolumn{10}{c}{Probabilit\\`a (in "
       "percentuale)} \\\\[\\belowrulesep]");
  puts("    \\cmidrule(l){2-11}");
  puts("    $N$");

  for (i = 0;  i < N_T;  i++) {
    printf("      & \\multicolumn{1}{c}{%4.1f}\n",
           t[i] * 100.0);
  }
  puts("      \\\\");
  puts("    \\midrule");

  for (ndf = 1;   ndf < 41;   ndf++) {
    if (ndf > 1   &&   ndf % 10 == 1) {
      puts("    \\cmidrule(lr){2-11}");
    }
    printf("    %d", ndf);

    for (i = 0;   i < N_T;   i++) {
      double x;

      x = gsl_cdf_tdist_Pinv(t[i], ndf);
      printf(" & %.3f", x);
      if (i == 5) printf("\n     ");
    }
    puts(" \\\\");
  }

  puts("    \\bottomrule");
  puts("    \\addlinespace");
  puts("    \\multicolumn{11}{c}{Percentili della distribuzione di Student}");
  puts("  \\end{tabular}");
  puts("\\end{center}");
  puts("\\clearpage");
}

static void t_chi(void)
{

  /**
   | Percentile tables of the chi-squared function
   | =============================================
   |
   | For the probability levels in the array c[N_C], inverts the
   | cumulative chi-squared distribution function; repeating all for a
   | number of degrees of freedom from 1 to 30.
  **/

  const double c[] = {
    0.9990, 0.9950, 0.9900, 0.9500, 0.9000,
    0.7500, 0.5000, 0.2500, 0.1000, 0.0500,
    0.0100, 0.0050, 0.0010
  };
  const int N_C = sizeof(c) / sizeof(c[0]);
  int       i, j;

  puts("  \\begin{center}");
  puts("    \\footnotesize");
  puts("    \\renewcommand{\\arraystretch}{0.89}");
  puts("    \\begin{tabular}{rrrrrrrrrrrrrr}");
  puts("      & \\multicolumn{13}{c}{Probabilit\\`a (in "
       "percentuale)} \\\\[\\belowrulesep]");
  puts("      \\cmidrule{2-14}");
  puts("      $N$");

  for (i = 0;   i < N_C;   i++) {
    printf("        & \\multicolumn{1}{c}{%.1f}\n",
           c[i] * 100.0);
  }
  puts("        \\\\");
  puts("      \\midrule");

  for (j = 0;   j < 30;   ) {
    printf("      %d", ++j);
    for (i = 0;   i < N_C;   i++) {
      double x;

      x = gsl_cdf_chisq_Pinv(c[i], j);
      printf(" & %.2f", x);
      if (i == 5) printf("\n       ");
    }
    puts(" \\\\");
    if (j % 10 == 0    &&   j != 30) {
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

static void t_fisher(void)
{

  /**
   | Percentile tables of the Fisher's function
   | ==========================================
   |
   | For the probability levels in the array f[N_F], inverts the
   | cumulative Fisher's distribution function; repeating all for the
   | degrees of freedom in m[N_FM] and n[N_FN].
   |
   | The GSL library does not include procedures that directly invert
   | the cumulative Fisher's distribution function; an auxiliary
   | function "Fisher's PDF minus wanted confidence level" is defined,
   | and its zeros are found using the GSL "one-dimensional root
   | finding" package.
  **/

  const double f[] = {
    0.95, 0.99
  };
  const int N_F = sizeof(f) / sizeof(f[0]);

  const int m[] = {
    1,  2,  3,  4,  5,  6,  8, 12, 24, 36, 48,
    1000000
  };
  const int N_FM = sizeof(m) / sizeof(m[0]);

  const int n[] = {
    1,  2,  3,  4,  5,  6,  7,  8,  9, 10,
    11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
    21, 22, 23, 24, 25, 30, 40, 60,120, 1000000
  };
  const int N_FN = sizeof(n) / sizeof(n[0]);

  int                          kf, km, kn;
  double                       xmax;
  fishParams                   fp;
  const gsl_root_fsolver_type *pGRT = FSOLVER_TYPE;
  gsl_root_fsolver            *pGRS;
  gsl_function                 fishFun;

  /**
   | Definition and initialization of the GSL root finding environment
  **/

  if ((pGRS = gsl_root_fsolver_alloc(pGRT)) == NULL) {
    fputs("Error return from gsl_root_fsolver_alloc\n", stderr);
    exit(EXIT_FAILURE);
  }

  fishFun.function = fishAux;
  fishFun.params   = &fp;

  for (kf = 0;   kf < N_F;   kf++) {
    puts("  \\clearpage");
    puts("  \\begin{center}");
    puts("    \\footnotesize");
    puts("    \\renewcommand{\\arraystretch}{0.91}");
    if (kf == 1) {
      puts("    \\addtolength{\\tabcolsep}{-0.14pt}");
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

    xmax      = 7000.0;
    fp.target = f[kf];

    for (kn = 0;   kn < N_FN;   ) {
      double new_xmax = 0.0;

      fp.ndg = n[kn];
      if (kn < N_FN-1) {
        printf("      %d", fp.ndg);
      } else {
        printf("      $\\infty$");
      }

      for (km = 0;   km < N_FM;   ) {
        int    iter = 0;
        int    status;
        double x, y;

        fp.ndf = m[km];
        gsl_root_fsolver_set(pGRS, &fishFun, 0.0, xmax);

        do {
          iter++;
          gsl_root_fsolver_iterate(pGRS);
          x      = gsl_root_fsolver_root(pGRS);
          y      = fishAux(x, &fp);
          status = gsl_root_test_residual(y, FISH_EPS);
        } while (status == GSL_CONTINUE   &&   iter < FISH_ITER);

        if (status != GSL_SUCCESS) {
          fprintf(stderr, "t_fisher: too many iterations for N, "
                  "M = %d %d\n", fp.ndf, fp.ndg);
          gsl_root_fsolver_free(pGRS);
          exit(EXIT_FAILURE);
        }

        printf(" & %.2f", x);
        if (x > new_xmax) new_xmax = x;
        if ((++km % 5) == 0) printf("\n       ");
      } /* Loop over M (columns) */

      puts(" \\\\");
      if ((++kn % 10) == 0) {
        if (kn < N_FN) {
          puts("      \\cmidrule(lr){2-13}");
        } else {
          puts("      \\bottomrule");
        }
      }
      xmax = new_xmax;
    } /* Loop over N (rows)  */

    puts("      \\addlinespace");
    puts("      \\multicolumn{13}{c}%");
    printf("        {Percentili della distribuzione di Fisher "
           "per $P=%.2f$}\n", f[kf]);
    puts("    \\end{tabular}");
    puts("  \\end{center}");
  } /* Loop over f (confidence levels) */

  gsl_root_fsolver_free(pGRS);
}

static double fishAux(
  double  x,
  void   *params
) {
  fishParams *pFP = (fishParams *) params;
  double      nu1 = pFP->ndf;
  double      nu2 = pFP->ndg;

  return gsl_cdf_fdist_P(x, nu1, nu2) - pFP->target;
}
