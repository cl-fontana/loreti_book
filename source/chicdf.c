/*------------------------------------------------------*
 | Author: Maurizio Loreti, aka MLO or (HAM) I3NOO      |
 | Work:   University of Padova - Department of Physics |
 |         Via F. Marzolo, 8 - 35131 PADOVA - Italy     |
 | Phone:  +39 (049) 827-7216   FAX: +39 (049) 827-7102 |
 | EMail:  loreti@pd.infn.it                            |
 | WWW:    http://www.pd.infn.it/~loreti/mlo.html       |
 *------------------------------------------------------*

 This C snippet generates a KUMAC macro suitable to obtain the figure
 32.1 of the chapter 32 ("Statistics") of the PDG paper (see at the
 URL http://pdg.lbl.gov/).  Compile/link with the following command:

 gcc -std=c99 -pedantic -W -Wall -O2 -o chi chi.c -lgsl -lgslcblas -lm

 The above command assumes you have installed the GNU Scientific
 Library (GSL): for more, see under http://www.gnu.org/software/gsl/ .
 The output in to the standard output stream; redirect to a file as in

 ./chi > chi.kumac

 THIS CODE IS PUBLIC DOMAIN.

 *------------------------------------------------------------------*
 $Id: chicdf.c,v 1.2 2005/03/10 12:15:29 loreti Exp $
 *------------------------------------------------------------------*/

#include <stdio.h>
#include <math.h>
#include <gsl/gsl_cdf.h>

/**
 | The function to be plotted will be computed in NPT points equally
 | spaced along the X-axis.
**/

#define NPT 100

int main() {

  /**
   | The function is the chi-square cumulative distribution function,
   | for all the degrees of freedom in the array "ndf"; "nndf" is
   | their number.
  **/

  int ndf[] = {1, 2, 4, 7, 12, 20, 35, 60};
  int nndf  = sizeof(ndf) / sizeof(ndf[0]);

  /**
   | xmin, xmax: lower and upper limit of X-axis in the plot
   | ymin, ymax: the same, for the Y-axis (cumulative distribution
   |             function)
   | up: setting xmax and ymax as upper axis limits, paw does not
   |     put the labels for these values on the plot.  For that
   |     reason, paw will be told that the upper limits are xmax*up
   |     and ymax*up.
   | xval: an array to save the x values for the points in all the
   |       plotted functions.
  **/

  double xmin = 0.1;
  double xmax = 100.0;
  double ymin = 0.001;
  double ymax = 1.0;
  double up   = 1.001;
  double xval[NPT+1];

  int    i, j;

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
  puts("set xlab 1.6");
  puts("set xval 0.6");
  puts("igset txfp -130");
  puts("opt logx");
  puts("opt logy");
  puts("set lwid 3\n");

  puts("fort/file 66 chicdf.eps");
  puts("meta 66 -113\n");

  /**
   | X values, i.e. NPT points equally spaced between xmin and xmax;
   | that will be computed, stored in xval[] for later use, and
   | written to stdout to define a paw vector.  Since the scale is
   | logarithmic, a small hack is required...
  **/

  printf("vect/create x(%d) r", NPT+1);

  for (i = 0;   i <= NPT;   ++i) {
    double logxmin = log(xmin);
    double x = logxmin +
               (log(xmax) - logxmin) * ((double) i / NPT);
    xval[i] = x = exp(x);

    if ((i % 5) == 0) puts(" _");
    printf("  %f", x);
  }
  puts("\n");

  for (j = 0;   j < nndf;   ++j) {
    int ny;
    int n = ndf[j];
    double yval[NPT+1];

    /**
     | "n" will loop on all the wanted degrees of freedom; for every
     | "n", the cumulative PDF will be computed for every x and
     | written to stdout as the components of a paw vector.  It is
     | assumed that all the plotted points are less than or equal to
     | ymax; and only the first point below ymin is written.
    **/

    for (ny = i = 0;   i <= NPT;   ++i) {
      double y = gsl_cdf_chisq_Q(xval[i], n);

      yval[ny++] = y;
      if (y < ymin) break;
    }

    printf("vect/create y%d(%d) r", n, NPT+1);

    for (i = 0;   i < ny;   ++i) {
      if ((i % 5) == 0) puts(" _");
      printf("  %f", yval[i]);
    }
    puts("\n");
  }

  /**
   | Plotting commands
  **/

  printf("hplot/null %f %f %f %f\n", xmin, xmax*up, ymin, ymax*up);

  for (j = 0;   j < nndf;   ++j) {
    printf("graph/prim/graph %d x y%d 'l'\n", NPT+1, ndf[j]);
  }

  puts("\nigset tang 90");
  puts("igset chhe 0.5");
  puts("itx 0.07 0.5 'P'");
  puts("igset tang");
  puts("itx 50 6e-4 'x'");
  puts("igset chhe 0.3\n");

  puts("itx 1.4 0.1 'N = 1'");
  puts("itx 3.8 0.1 '2'");
  puts("itx 6.5 0.1 '4'");
  puts("itx 10 0.1 '7'");
  puts("itx 14 0.1 '12'");
  puts("itx 21 0.1 '20'");
  puts("itx 34 0.1 '35'");
  puts("itx 55 0.1 '60'");

  puts("\nclose 66");
  return 0;
}
