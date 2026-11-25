#include "main.h"

// One piece of a trace estimate. A piece can combine several independent estimates
// (e.g. a connected G_{i,j}, averaged over the inner timeslices t'):
//   value = sum_k w_k est_k,   var of estimate = sum_k w_k^2 s_k^2 / N_k.
#define TRACE_REPORT_MAX_PIECES 512
#define TRACE_REPORT_LABEL 96

typedef struct {
  char label[TRACE_REPORT_LABEL];
  complex_double value;
  int n_est, n_samples;           // estimates combined, samples per estimate (last one)
  double vs_re, vs_im;            // sum of per-sample variances (averaged when printed)
  double vm_re, vm_im;            // variance of the estimate of the piece
  int undefined;                  // some estimate had a single sample: variance unknown
} trace_piece;

static trace_piece pieces[TRACE_REPORT_MAX_PIECES];
static int num_pieces = 0;

void trace_report_reset( void ){
  num_pieces = 0;
}

// var_re, var_im: per-sample (unbiased) variances of the real and imaginary parts,
// negative if undefined (a single sample). weight: factor of this estimate in the piece.
void trace_report_add( const char *label, complex_double value, int n_samples,
                       double var_re, double var_im, double weight ){
  int k;
  for( k=0; k<num_pieces; k++ )
    if( strcmp( pieces[k].label, label ) == 0 ) break;
  if( k == num_pieces ){
    if( num_pieces == TRACE_REPORT_MAX_PIECES ){
      printf("trace_report: more than %d pieces, %s not recorded\n", TRACE_REPORT_MAX_PIECES, label);
      return;
    }
    memset( &pieces[k], 0, sizeof(trace_piece) );
    snprintf( pieces[k].label, TRACE_REPORT_LABEL, "%s", label );
    num_pieces++;
  }
  trace_piece *p = &pieces[k];
  p->value += weight*value;
  p->n_est++;
  p->n_samples = n_samples;
  if( var_re < 0 || var_im < 0 || n_samples < 2 ){
    p->undefined = 1;
    return;
  }
  p->vs_re += var_re;
  p->vs_im += var_im;
  p->vm_re += weight*weight*var_re/n_samples;
  p->vm_im += weight*weight*var_im/n_samples;
}

// total = |.|^2 variance (real + imaginary); the standard error of the total is its square root
static void print_variances( const char *what, double v_re, double v_im ){
  printf("    %-22s re %.6e   im %.6e   total %.6e\n", what, v_re, v_im, v_re+v_im);
}

static void print_std_errors( double v_re, double v_im ){
  printf("    %-22s re %.6e   im %.6e   total %.6e\n", "standard error", sqrt(v_re), sqrt(v_im), sqrt(v_re+v_im));
}

void trace_report_print( const char *title ){
  complex_double total = 0;
  double tv_re = 0, tv_im = 0;
  int undefined = 0;

  printf("\n==================== Trace summary: %s ====================\n", title);
  for( int k=0; k<num_pieces; k++ ){
    trace_piece *p = &pieces[k];
    total += p->value;
    tv_re += p->vm_re; tv_im += p->vm_im;
    undefined |= p->undefined;

    printf("\n%s\n", p->label);
    printf("    trace                  %.10e %c i%.10e\n", creal(p->value),
           (cimag(p->value) < 0) ? '-' : '+', fabs(cimag(p->value)));
    if( p->n_est == 1 ) printf("    samples                %d\n", p->n_samples);
    else                printf("    samples                %d estimates x %d\n", p->n_est, p->n_samples);
    if( p->undefined ){
      printf("    variance               n/a (single sample)\n");
      continue;
    }
    print_variances( (p->n_est == 1) ? "variance per sample" : "var per sample (avg)",
                     p->vs_re/p->n_est, p->vs_im/p->n_est );
    print_variances( "variance of estimate", p->vm_re, p->vm_im );
    print_std_errors( p->vm_re, p->vm_im );
  }

  printf("\nTOTAL (%d piece%s)\n", num_pieces, (num_pieces == 1) ? "" : "s");
  printf("    trace                  %.10e %c i%.10e\n", creal(total),
         (cimag(total) < 0) ? '-' : '+', fabs(cimag(total)));
  if( undefined ) printf("    (pieces with a single sample are left out of the variance)\n");
  print_variances( "variance", tv_re, tv_im );
  print_std_errors( tv_re, tv_im );
  printf("==========================================================================\n\n");
  fflush(0);
}
