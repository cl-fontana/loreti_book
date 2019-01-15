/*------------------------------------------------------*
 | Author: Maurizio Loreti, aka MLO or (HAM) I3NOO      |
 | Work:   University of Padova - Department of Physics |
 |         Via F. Marzolo, 8 - 35131 PADOVA - Italy     |
 | Phone:  +39 (049) 827-7216   FAX: +39 (049) 827-7102 |
 | EMail:  loreti@pd.infn.it                            |
 | WWW:    http://www.pd.infn.it/~loreti/mlo.html       |
 *------------------------------------------------------*

 This C snippet generates a KUMAC macro suitable to obtain an old
 figure (now removed) of the chapter 32 ("Statistics") of the PDG
 paper (see at the URL http://pdg.lbl.gov/).  Compile/link with the
 following command:

 gcc -std=c99 -pedantic -W -Wall -O2 -o chi chi.c -lgsl -lgslcblas -lm

 The above command assumes you have installed the GNU Scientific
 Library (GSL): for more, see under http://www.gnu.org/software/gsl/ .
 The output in to the standard output stream; redirect to a file as in

 ./chi > chi.kumac

 THIS CODE IS PUBLIC DOMAIN.

 *------------------------------------------------------------------*
 $Id: chirid.c,v 1.1 2005/03/10 12:15:29 loreti Exp $
 *------------------------------------------------------------------*/

#include <stdio.h>
#include <math.h>
#include <gsl/gsl_cdf.h>

/**
 | The function to be plotted will be computed with values from 1 to
 | NDFMAX (inclusive) of the number of degrees of freedom; and for the
 | confidence levels defined in the array conf[0] ... conf[nconf-1] .
**/

#define NDFMAX 50

double conf[] = {0.01, 0.05, 0.10, 0.20, 0.30,
                 0.50, 0.68, 0.90, 0.95, 0.997};
int nconf     = sizeof(conf) / sizeof(conf[0]);

int main() {

  /**
   | xmin, xmax: lower and upper limit of X-axis in the plot
   | ymin, ymax: the same, for the Y-axis (cumulative distribution
   |             function)
  **/

  double xmin = 0.0;
  double xmax = NDFMAX;
  double ymin = 0.0;
  double ymax = 2.5;
  int    i;

  /**
   | Preamble
  **/

  printf("| Generated on %s, from %s\n\n", __DATE__, __FILE__);

  puts("vect/del *");
  puts("opt *");
  puts("set *");
  puts("igset *\n");

  puts("opt nbox");
  puts("opt tic");
  puts("opt grid");
  puts("set *fon -130");
  puts("set *siz 0.4");
  puts("igset txfp -130");
  puts("set lwid 3\n");

  puts("fort/file 66 chirid.eps");
  puts("meta 66 -113\n");

  printf("sigma x=array(%d,1#%d)\n\n", NDFMAX, NDFMAX);

  for (i = 0;   i < nconf;   ++i) {

    /**
     | Loop over the confidence levels; "prob" is the current value.
    **/

    double prob = conf[i];
    double chirid[NDFMAX];
    int    ndf;

    for (ndf = 1;   ndf <= NDFMAX;   ++ndf) {

      /**
       | Loop over the degrees of freedom; "ndf" is the current value.
       |
       | Finds the value X for which the integral, from X to +\infty,
       | of the chi-square probability distribution function with
       | "ndf" degrees of freedom is equal to "prob"; then stores in
       | the array "chirid" the corresponding reduced value, X/ndf.
      **/

      double X        = gsl_cdf_chisq_Qinv(prob, ndf);
      chirid[ndf - 1] = X / ndf;
    }
    printf("vect/create y%d(%d) r", i, NDFMAX);

    for (ndf = 0;   ndf < NDFMAX;   ++ndf) {
      if ((ndf % 5) == 0) puts(" _");
      printf("  %f", chirid[ndf]);
    }
    puts("\n");
  }

  /**
   | Plotting commands
  **/

  puts("set ndvx 50205");
  puts("set ndvy 505");
  printf("hplot/null %f %f %f %f\n\n", xmin, xmax, ymin, ymax);

  for (i = 0;   i < nconf;   ++i) {
    printf("graph/prim/graph %d x y%d 'c'\n", NDFMAX, i);
  }

  puts("\nigset tang 90");
  puts("igset chhe 0.5");
  puts("itx -3.5 2.15 '[h]^2!/N'");
  puts("igset tang");
  puts("itx 46 -0.18 'N'");
  puts("igset chhe 0.3\n");

  puts("itx 11 2.27 '1%'");
  puts("itx 11 1.81 '5%'");
  puts("itx 11 1.58 '10%'");
  puts("itx 11 1.34 '20%'");
  puts("itx 11 1.10 '30%'");
  puts("itx 11 0.88 '50%'");
  puts("itx 11 0.70 '68%'");
  puts("itx 11 0.55 '90%'");
  puts("itx 11 0.36 '95%'");
  puts("itx 11 0.15 '99.7%'");

  puts("\nclose 66");
  return 0;
}
