#ifndef TRACE_REPORT_HEADER
  #define TRACE_REPORT_HEADER

  // Summary of the pieces of a trace estimate (one row per piece) and of the total.
  // Pieces are recorded by rank 0, master thread only.
  void trace_report_reset( void );
  void trace_report_add( const char *label, complex_double value, int n_samples,
                         double var_re, double var_im, double weight );
  void trace_report_print( const char *title );

#endif
