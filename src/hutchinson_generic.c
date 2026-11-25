#include "main.h"


struct estimate {
  int counter; //required number of estimates.
  complex_PRECISION estimate;
};

struct sample {
  // required number of estimates
  int sample_size;
  // accumulated trace
  complex_PRECISION acc_trace;
  // number of noise vectors and per-sample (unbiased) variances of the real and
  // imaginary parts of the samples; negative if undefined (one sample)
  int n_samples;
  double var_re, var_im;
};

// Record one estimate as (part of) a piece of the trace summary, with factor weight.
static void report_piece_PRECISION( const char *label, struct sample estimate, double weight,
                                    struct Thread *threading ){
  START_MASTER(threading)
  if( g.my_rank == 0 )
    trace_report_add( label, (complex_double)( estimate.acc_trace / estimate.sample_size ),
                      estimate.n_samples, estimate.var_re, estimate.var_im, weight );
  END_MASTER(threading)
}


void hutchinson_diver_PRECISION_init( level_struct *l, struct Thread *threading ) {
  hutchinson_PRECISION_struct* h = &(l->h_PRECISION);

  h->max_iters = NULL;
  h->min_iters = NULL;

  // MLMC
  h->mlmc_b1 = NULL;
  h->mlmc_b2 =NULL;
  h->mlmc_testing =NULL;
  h->rademacher_vector =NULL;
  h->rademacher_buffer =NULL;
  h->hadamard_vector =NULL;
  h->probing_vector =NULL;

  h->finest_level = l;
  h->lx_i = NULL;
  h->lx_j = NULL;

  SYNC_MASTER_TO_ALL(threading)
}


void hutchinson_diver_PRECISION_alloc( level_struct *l, struct Thread *threading ) {
  int i;
  hutchinson_PRECISION_struct* h = &(l->h_PRECISION) ;

  PUBLIC_MALLOC( h->max_iters, int, g.num_levels );
  PUBLIC_MALLOC( h->min_iters, int, g.num_levels );

  // For MLMC
  PUBLIC_MALLOC( h->mlmc_b1, complex_PRECISION, l->inner_vector_size );
  PUBLIC_MALLOC( h->mlmc_b2, complex_PRECISION, l->inner_vector_size );
  PUBLIC_MALLOC( h->mlmc_testing, complex_PRECISION, l->inner_vector_size );
  PUBLIC_MALLOC( h->rademacher_vector, complex_PRECISION, l->inner_vector_size );
  PUBLIC_MALLOC( h->rademacher_buffer, complex_PRECISION, l->inner_vector_size );
  PUBLIC_MALLOC( h->hadamard_vector, complex_PRECISION, l->inner_vector_size );
  PUBLIC_MALLOC( h->probing_vector, complex_PRECISION, l->inner_vector_size );

  for ( i=0;i<g.num_levels;i++ ) {
    h->max_iters[i] = g.trace_max_iters[i];
    h->min_iters[i] = g.trace_min_iters[i];
  }
}


void hutchinson_diver_PRECISION_free( level_struct *l, struct Thread *threading ) {
  hutchinson_PRECISION_struct* h = &(l->h_PRECISION) ;

  PUBLIC_FREE( h->max_iters, int, g.num_levels );
  PUBLIC_FREE( h->min_iters, int, g.num_levels );

  PUBLIC_FREE( h->mlmc_b1, complex_PRECISION, l->inner_vector_size );
  PUBLIC_FREE( h->mlmc_b2, complex_PRECISION, l->inner_vector_size );
  PUBLIC_FREE( h->mlmc_testing, complex_PRECISION, l->inner_vector_size );
  PUBLIC_FREE( h->rademacher_vector, complex_PRECISION, l->inner_vector_size );
  PUBLIC_FREE( h->rademacher_buffer, complex_PRECISION, l->inner_vector_size );
  PUBLIC_FREE( h->hadamard_vector, complex_PRECISION, l->inner_vector_size );
  PUBLIC_FREE( h->probing_vector, complex_PRECISION, l->inner_vector_size );
}




complex_PRECISION hutchinson_driver_PRECISION( level_struct *l, struct Thread *threading ){
  complex_PRECISION trace = 0.0;
  struct sample estimate;
  hutchinson_PRECISION_struct* h = &(l->h_PRECISION);
  level_struct* lx = l;

  // set the pointer to the finest-level Hutchinson estimator
  h->hutch_compute_one_sample = hutchinson_plain_PRECISION;

  estimate = probing_hutchinson_blind_PRECISION( lx, h, 0, 0, threading );
  trace += estimate.acc_trace / estimate.sample_size;
  report_piece_PRECISION( "plain Hutchinson (level 1)", estimate, 1.0, threading );

  return trace;
}

// Level on whose grid the stochastic vectors live for a given estimator term:
// - type 1 : vectors on l->next_level (split intermediate terms)
// - otherwise : vectors on level l (4D and 3D MGMLMC terms live on their own level)
level_struct* trace_vector_level_PRECISION( level_struct *l, hutchinson_PRECISION_struct* h, int type ){
  if( type==1 ) return l->next_level;
  if( type!=0 ) error0("Unknown value for type of Rademacher vector in relation to level of creation\n");
  return l;
}

// Build in h->probing_vector the probing vector of the current g.coloring_count and
// g.dilution_count on the grid of level vl: plain dilution (probing 0), multiplier-based
// coloring (probing 1) or Hadamard vector (probing 2). dil_level selects the dilution scheme.
void probe_create_PRECISION( level_struct *vl, hutchinson_PRECISION_struct* h, int dil_level, struct Thread *threading ){
  START_MASTER(threading)
  int mode = level_probing( vl->depth );
  if( mode == 0 )
    vector_PRECISION_dilution( h->probing_vector, 0, vl->inner_vector_size, vl, dil_level );
  else if( mode == 1 )
    vector_PRECISION_probing( h->probing_vector, 0, vl->inner_vector_size, vl, dil_level );
  else if( mode == 2 )
    vector_PRECISION_hadamard( h->probing_vector, 0, vl->inner_vector_size, vl, dil_level );
  else
    error0("Unknown probing selection %d\n", g.probing);
  END_MASTER(threading)
  SYNC_MASTER_TO_ALL(threading)
}

// Fill h->rademacher_vector with Z2 noise on the grid of level vl
void rademacher_create_PRECISION( level_struct *vl, hutchinson_PRECISION_struct* h, struct Thread *threading ){
  START_MASTER(threading)
  vector_PRECISION_define_random_rademacher( h->rademacher_vector, 0, vl->inner_vector_size, vl );
  END_MASTER(threading)
  SYNC_MASTER_TO_ALL(threading)
}

// FIXME : hacked this function a bit to avoid compiler warnings
gmres_PRECISION_struct* get_p_struct_PRECISION( level_struct* l ){
  if( l->depth==0 ){
    if ( strcmp("PRECISION","float")==0 ) { return &(l->p_PRECISION); }
    else { return (gmres_PRECISION_struct*)&(g.p); }
  }
  else{ return &(l->p_PRECISION); }
}


// Inner product on the grid of level l for the trace estimators. On coarse grids,
// processes can be merged: idle processes take no part in that level's reduction
// and receive the result from rank 0 (active on every level), so that every rank
// accumulates the same sample.
complex_PRECISION trace_inner_product_PRECISION( vector_PRECISION a, vector_PRECISION b, int start, int end,
                                                 level_struct *l, struct Thread *threading ){
  static complex_PRECISION shared;
  complex_PRECISION aux = 0;

  if( l->depth == 0 || g.num_processes == 1 )
    return global_inner_product_PRECISION( a, b, start, end, l, threading );

  if( !l->idle ) aux = global_inner_product_PRECISION( a, b, start, end, l, threading );
  START_MASTER(threading)
  if( g.my_rank == 0 && l->idle ) error0("trace_inner_product: rank 0 idle at depth %d\n", l->depth);
  shared = aux;
  MPI_Bcast( &shared, 1, MPI_COMPLEX_PRECISION, 0, g.comm_cart );
  END_MASTER(threading)
  SYNC_MASTER_TO_ALL(threading)
  aux = shared;
  SYNC_CORES(threading)
  return aux;
}

int apply_solver_PRECISION( level_struct* l, struct Thread *threading ){
  int nr_iters = 0;
  double buff1=0, buff2=0;

  if( l->idle ) return 0;   // processes merged away on this level take no part in its solves

  gmres_PRECISION_struct* p = get_p_struct_PRECISION( l );

  // force a zero initial guess: prevents a previous (possibly bad) solution
  // from contaminating this solve
  {
    int start0, end0;
    compute_core_start_end( 0, l->inner_vector_size, &start0, &end0, l, threading );
    vector_PRECISION_define( p->x, 0, start0, end0, l );
    SYNC_CORES(threading)
  }

  p->print_iters = 1;

  buff1 = p->tol;
  p->tol = g.tol;
  START_MASTER(threading);
  if( l->level==0 ){
    buff2 = g.coarse_tol;
    g.coarse_tol = g.tol;
  }
  END_MASTER(threading);
  SYNC_MASTER_TO_ALL(threading);

  if ( l->level > 0 ) {
    nr_iters = fgmres_PRECISION( p, l, threading );
  } else {
#ifdef GCRODR
    // NOTE : something that shouldn't be happening here happens, namely the RHS is changed
    //        by the function coarse_solve_odd_even_PRECISION(...). So, we back it up and restore
    //        it as necessary

    int start,end;
    //compute_core_start_end( l->p_PRECISION.v_start, l->p_PRECISION.v_end, &start, &end, l, threading );
    compute_core_start_end( 0, l->inner_vector_size, &start, &end, l, threading );
    vector_PRECISION_copy( l->p_PRECISION.rhs_bk, l->p_PRECISION.b, start, end, l );

    START_MASTER(threading)
    l->p_PRECISION.was_there_stagnation = 0;
    END_MASTER(threading)
    SYNC_MASTER_TO_ALL(threading)

    while( 1 ) {
      coarse_solve_odd_even_PRECISION( &(l->p_PRECISION), &(l->oe_op_PRECISION), l, threading );
      if ( l->p_PRECISION.was_there_stagnation==0 ) { break; }
      else if ( l->p_PRECISION.was_there_stagnation==1 && l->p_PRECISION.gcrodr_PRECISION.CU_usable==1 ) {
        // in case there was stagnation, we need to rebuild the coarsest-level data
        double time_bk = g.coarsest_time;
        coarsest_level_resets_PRECISION( l, threading );
        START_MASTER(threading)
        l->p_PRECISION.was_there_stagnation = 0;
        g.coarsest_time = time_bk;
        END_MASTER(threading)
        SYNC_MASTER_TO_ALL(threading)
        vector_PRECISION_copy( l->p_PRECISION.b, l->p_PRECISION.rhs_bk, start, end, l );
      }
      else {
        // in this case, there was stagnation but no deflation/recycling subspace is being used
        break;
      }
    }
#else
    coarse_solve_odd_even_PRECISION( &(l->p_PRECISION), &(l->oe_op_PRECISION), l, threading );
#endif
  }

  START_MASTER(threading);
  p->tol = buff1;
  if( l->level==0 ){
    g.coarse_tol = buff2;
  }
  END_MASTER(threading);
  SYNC_MASTER_TO_ALL(threading);

  p->print_iters = 0;

  return nr_iters;
}

void print_variance_pc_PRECISION(level_struct *l, int iters){

  double var_avg = 0.0;
  double var_var = 0.0;

  for(int i = 0; i < iters; i++)
    var_avg += g.var_pc[l->depth][i];

  var_avg = var_avg/iters;

  for(int i = 0; i < iters; i++){
    var_var += (g.var_pc[l->depth][i] - var_avg)*(g.var_pc[l->depth][i] - var_avg);
  }

  var_var = var_var/iters;

  double interval;

  interval = sqrt(var_var/iters);

  if(g.probing != 0) printf("\nVariance of color %d - dilution index %d = %f +- %f\n", g.coloring_count, g.dilution_count, var_avg, interval);
  if(g.probing == 0) printf("\nVariance of dilution index %d = %f +- %f\n", g.dilution_count, var_avg, interval);

  for(int i = 0; i < g.trace_max_iters[l->depth]; i ++)
    g.var_pc[l->depth][i] = 0.0;

}

// Stochastic estimation of the trace of the operator applied by h->hutch_compute_one_sample.
//
// - type 0 : noise vectors on level l (or on the finest grid for the 3D MGMLMC kernels),
//   type 1 : noise vectors on l->next_level (see trace_vector_level_PRECISION)
// - dil_level : level whose dilution scheme (g.dilution[dil_level]) is used
//
// For every sample ONE Z2 vector z0 is drawn and shared by all probing vectors
// (colors x dilution pieces) of that sample: v = z0 .* probe. The sample is
//   probing 0 (plain) and 1 (multiplier): sum over colors and dilution pieces
//   probing 2 (hierarchical)            : sum over dilution pieces, average over
//                                         the Hadamard vectors.
// The returned acc_trace/sample_size is the estimate of the trace.
struct sample probing_hutchinson_blind_PRECISION( level_struct *l, hutchinson_PRECISION_struct* h, int type,
                                                  int dil_level, struct Thread *threading ){
  int i, j, start, end;
  complex_PRECISION one_sample=0.0, variance=0.0, trace=0.0;
  struct sample estimate;

  level_struct *vl = trace_vector_level_PRECISION( l, h, type );
  int num_colors = g.num_colors[vl->depth];
  int hp         = ( level_probing( vl->depth ) == 2 );   // hierarchical probing on this grid
  int num_dil    = g.dilution[dil_level];
  int max_iters  = h->max_iters[l->depth];
  // index of g.variances where the variance of this estimator is stored
  int var_idx    = ( g.trace_op_type == 7 ) ? h->lx_i->depth*g.num_levels + h->lx_j->depth : l->depth;

  complex_PRECISION* samples = (complex_PRECISION*) malloc( max_iters*sizeof(complex_PRECISION) );
  complex_PRECISION* traces  = (complex_PRECISION*) malloc( num_colors*sizeof(complex_PRECISION) );
  for( i=0; i<max_iters; i++ ) samples[i] = 0.0;

  compute_core_start_end( 0, vl->inner_vector_size, &start, &end, vl, threading );

  estimate.acc_trace = 0.0;
  double t0 = MPI_Wtime();

  for( i=0; i<max_iters; i++ ){
    // 1. one Z2 vector z0 per sample, kept in h->rademacher_buffer
    rademacher_create_PRECISION( vl, h, threading );
    vector_PRECISION_copy( h->rademacher_buffer, h->rademacher_vector, start, end, vl );
    SYNC_CORES(threading)

    for( j=0; j<num_colors; j++ ) traces[j] = 0.0;

    for( g.coloring_count=0; g.coloring_count<num_colors; g.coloring_count++ ){
      for( g.dilution_count=1; g.dilution_count<num_dil+1; g.dilution_count++ ){

        // 2. v = z0 .* probe
        probe_create_PRECISION( vl, h, dil_level, threading );
        hadamard_PRECISION_product( h->rademacher_vector, h->probing_vector, start, end, vl );
        SYNC_CORES(threading)

        // 3. apply the operator and take the inner product
        one_sample = h->hutch_compute_one_sample( -1, l, h, threading );

        START_MASTER(threading)
        if(g.my_rank==0) {
          printf("[%d, color %d, dof %d out of %d, trace: %e %c i%e]\n",
                 i, g.coloring_count+1, g.dilution_count, num_dil, creal(one_sample),
                 (cimag(one_sample) < 0) ? '-' : '+', fabs(cimag(one_sample)));
          fflush(0);
        }
        END_MASTER(threading)

        traces[g.coloring_count] += one_sample;
        estimate.acc_trace += one_sample;

        // restore z0 for the next probing vector
        vector_PRECISION_copy( h->rademacher_vector, h->rademacher_buffer, start, end, vl );
        SYNC_CORES(threading)
      }

      START_MASTER(threading)
      if(g.my_rank==0){
        printf("\n[It %d - Trace of h=%d: %e %c i%e]\n",
               i, g.coloring_count+1, creal(traces[g.coloring_count]),
               (cimag(traces[g.coloring_count]) < 0) ? '-' : '+', fabs(cimag(traces[g.coloring_count])));
        fflush(0);
      }
      END_MASTER(threading)
    }

    // 4. sample of this noise vector
    for( j=0; j<num_colors; j++ )
      samples[i] += traces[j];
    if( hp ) samples[i] = samples[i]/num_colors;

    START_MASTER(threading)
    if(g.my_rank==0) printf("\n[Trace of it %d = %e %c i%e]\n", i, creal(samples[i]),
                            (cimag(samples[i]) < 0) ? '-' : '+', fabs(cimag(samples[i])));
    END_MASTER(threading)

    // 5. running mean and unbiased variance over the samples so far
    if( i != 0 ){
      trace = 0.0;
      for( j=0; j<=i; j++ ) trace += samples[j];
      trace = trace/(i+1);

      variance = 0.0;
      for( j=0; j<=i; j++ ) variance += conj(samples[j] - trace)*(samples[j] - trace);
      variance = variance/i;

      START_MASTER(threading)
      if(g.my_rank==0) {
        printf("\n[%d, trace: %e %c i%e, variance: %e]\n", i, creal(trace),
               (cimag(trace) < 0) ? '-' : '+', fabs(cimag(trace)), creal(variance));
        g.var_pc[l->depth][i] = creal(variance);
      }
      END_MASTER(threading)
    }
  }

  START_MASTER(threading)
  if(g.my_rank==0) {
    g.variances[var_idx] += creal(variance);
    printf("\nTime for sample computation (Avg.): \t %f\n\n", (MPI_Wtime()-t0)/max_iters);
  }
  END_MASTER(threading)

  // acc_trace sums ALL quadratures: normalize to the per-sample convention above
  estimate.sample_size = hp ? max_iters*num_colors : max_iters;

  // per-sample variances of the real and imaginary parts, for the trace summary
  estimate.n_samples = max_iters;
  estimate.var_re = -1.0; estimate.var_im = -1.0;
  if( max_iters > 1 ){
    complex_PRECISION mean = 0.0;
    for( j=0; j<max_iters; j++ ) mean += samples[j];
    mean = mean/max_iters;
    estimate.var_re = 0.0; estimate.var_im = 0.0;
    for( j=0; j<max_iters; j++ ){
      estimate.var_re += ( creal(samples[j]) - creal(mean) )*( creal(samples[j]) - creal(mean) );
      estimate.var_im += ( cimag(samples[j]) - cimag(mean) )*( cimag(samples[j]) - cimag(mean) );
    }
    estimate.var_re /= ( max_iters - 1 );
    estimate.var_im /= ( max_iters - 1 );
  }

  free(samples);
  free(traces);

  return estimate;
}

// this is the driver for plain Hutchinson
complex_PRECISION gamma_3D_hutchinson_driver_PRECISION( level_struct *l, struct Thread *threading ){
  complex_PRECISION trace = 0.0;
  struct sample estimate;
  hutchinson_PRECISION_struct* h = &(l->h_PRECISION);
  level_struct* lx = l;

  // set the pointer to the finest-level Hutchinson estimator
  h->hutch_compute_one_sample = gamma_3D_hutchinson_plain_PRECISION;

  estimate = probing_hutchinson_blind_PRECISION( lx, h, 0, 0, threading );
  trace += estimate.acc_trace / estimate.sample_size;
  report_piece_PRECISION( "plain Hutchinson (level 1)", estimate, 1.0, threading );

  // if deflation vectors are available
  //if(g.trace_deflation_type[l->depth] != 0){
  //  trace += hutchinson_deflated_direct_term_PRECISION( l, h, threading );
  //}

  return trace;
}

// this is the driver for plain Hutchinson, BUT for the connected operator (-> connected diagram)
complex_PRECISION g5_3D_connected_hutchinson_driver_PRECISION( level_struct *l, struct Thread *threading ){
  complex_PRECISION trace = 0.0;
  struct sample estimate;
  hutchinson_PRECISION_struct* h = &(l->h_PRECISION);
  level_struct* lx = l;

  // set the pointer to the finest-level Hutchinson estimator
  h->hutch_compute_one_sample = g5_3D_connected_hutchinson_plain_PRECISION;
  for ( g.time_slice_inner_connected=0; g.time_slice_inner_connected<g.global_lattice[0][0]; g.time_slice_inner_connected++ ) {
    
    // colors and dilution are looped inside the blind routine, with one shared z0 per sample
    estimate = probing_hutchinson_blind_PRECISION( lx, h, 0, 0, threading );
    trace += estimate.acc_trace / estimate.sample_size;
    report_piece_PRECISION( "connected plain Hutchinson (sum over t')", estimate, 1.0, threading );
    // estimate = hutchinson_blind_PRECISION( lx, h, 0, threading );
    // trace += estimate.acc_trace/estimate.sample_size;
  }

  // if deflation vectors are available
  //if(g.trace_deflation_type[l->depth] != 0){
  //  trace += hutchinson_deflated_direct_term_PRECISION( l, h, threading );
  //}

  return trace;
}


complex_PRECISION g5_3D_connected_hutchinson_plain_PRECISION( int type_appl, level_struct *l, hutchinson_PRECISION_struct* h, struct Thread *threading ){
  {
    int start, end;
    gmres_PRECISION_struct* p = get_p_struct_PRECISION( l );
    compute_core_start_end( 0, l->inner_vector_size, &start, &end, l, threading );

    if ( type_appl==-1 ) {
      //Apply \Pi_{t’+t}
      int bufft = g.time_slice;
      // TODO : check : is this assuming periodic in time ?
      g.time_slice = g.time_slice + g.time_slice_inner_connected;
      g.time_slice = g.time_slice%g.global_lattice[0][0];
      vector_PRECISION_ghg( h->rademacher_vector, 0, l->inner_vector_size, l );
      g.time_slice = bufft;
      vector_PRECISION_copy( p->b, h->rademacher_vector, start, end, l );
    } else {
      //vector_PRECISION_copy( p->b, l->powerit_PRECISION.vecs[type_appl], start, end, l );
    }

    // Apply Gamma_5
    gamma5_PRECISION( p->b, p->b, l, threading );
  }

  {
    apply_solver_PRECISION( l, threading );
  }

  {
    int start, end;
    gmres_PRECISION_struct* p = get_p_struct_PRECISION( l );
    compute_core_start_end( 0, l->inner_vector_size, &start, &end, l, threading );

    // apply \Pi_{t'}
    //vector_PRECISION_ghg(  p->x, 0, l->inner_vector_size, l );
    int bufft = g.time_slice;
    // TODO : check : is this assuming periodic in time ?
    g.time_slice = g.time_slice_inner_connected;
    //g.time_slice = g.time_slice%g.global_lattice[0][0];
    //vector_PRECISION_ghg( h->rademacher_vector, 0, l->inner_vector_size, l );
    vector_PRECISION_ghg( p->x, 0, l->inner_vector_size, l );
    g.time_slice = bufft;

    // apply G5
    gamma5_PRECISION( p->x, p->x, l, threading );

    // inverse again
    vector_PRECISION_copy( p->b, p->x, start, end, l );
    apply_solver_PRECISION( l, threading );
  }

  // subtract the results and perform dot product
  {
    int start, end;
    complex_PRECISION aux = 0;
    gmres_PRECISION_struct* p = get_p_struct_PRECISION( l );
    compute_core_start_end( 0, l->inner_vector_size, &start, &end, l, threading );

    if ( type_appl==-1 ) {
    //if(g.trace_deflation_type[l->depth] != 0){
    //  hutchinson_deflate_vector_PRECISION(p->x, l, threading);
    //}
      aux = trace_inner_product_PRECISION( h->rademacher_vector, p->x, 0, l->inner_vector_size, l, threading );
    } else {
      //vector_PRECISION_copy(l->powerit_PRECISION.vecs_buff1, p->x, start, end, l);
      //aux = global_inner_product_PRECISION( l->powerit_PRECISION.vecs[type_appl], l->powerit_PRECISION.vecs_buff1, p->v_start, p->v_end, l, threading );
    }

    return aux;
  }
}

complex_PRECISION gamma_3D_hutchinson_plain_PRECISION( int type_appl, level_struct *l, hutchinson_PRECISION_struct* h, struct Thread *threading ){
  {
    int start, end;
    gmres_PRECISION_struct* p = get_p_struct_PRECISION( l );
    compute_core_start_end( 0, l->inner_vector_size, &start, &end, l, threading );

    if ( type_appl==-1 ) {
      
      SYNC_CORES(threading)
      START_MASTER(threading)
      vector_PRECISION_ghg( h->rademacher_vector, 0, l->inner_vector_size, l );
      END_MASTER(threading)
      SYNC_MASTER_TO_ALL(threading)

      } else {
      //vector_PRECISION_copy( p->b, l->powerit_PRECISION.vecs[type_appl], start, end, l );
    }

    // Apply Gamma out-of-place: b = gamma * z.
    // NEVER apply spin-off-diagonal gammas in place (out==in overwrites
    // spin components before their partners are read).
    if(g.gamma_idx == 5)      gamma5_PRECISION( p->b, h->rademacher_vector, l, threading );
    else if(g.gamma_idx == 0) gamma0_PRECISION( p->b, h->rademacher_vector, l, threading );
    else if(g.gamma_idx == 1) gamma1_PRECISION( p->b, h->rademacher_vector, l, threading );
    else if(g.gamma_idx == 2) gamma2_PRECISION( p->b, h->rademacher_vector, l, threading );
    else if(g.gamma_idx == 3) gamma3_PRECISION( p->b, h->rademacher_vector, l, threading );
    else vector_PRECISION_copy( p->b, h->rademacher_vector, start, end, l );  // identity

  }

  {
    apply_solver_PRECISION( l, threading );
  }

  // subtract the results and perform dot product
  {
    int start, end;
    complex_PRECISION aux = 0;
    gmres_PRECISION_struct* p = get_p_struct_PRECISION( l );
    compute_core_start_end( 0, l->inner_vector_size, &start, &end, l, threading );

    if ( type_appl==-1 ) {
    //if(g.trace_deflation_type[l->depth] != 0){
    //  hutchinson_deflate_vector_PRECISION(p->x, l, threading);
    //}
      aux = trace_inner_product_PRECISION( h->rademacher_vector, p->x, 0, l->inner_vector_size, l, threading );
    } else {
      //vector_PRECISION_copy(l->powerit_PRECISION.vecs_buff1, p->x, start, end, l);
      //aux = global_inner_product_PRECISION( l->powerit_PRECISION.vecs[type_appl], l->powerit_PRECISION.vecs_buff1, p->v_start, p->v_end, l, threading );
    }

    return aux;
  }
}


// apply the interpolation
void apply_P_PRECISION( vector_PRECISION out, vector_PRECISION in, level_struct* l, struct Thread *threading ){
  if( l->idle ) return;   // l is the fine level; idle processes do not take part
  if( l->depth==0 ){
    interpolate3_PRECISION( l->sbuf_PRECISION[0], in, l, threading );
    trans_back_PRECISION( (vector_double)out, l->sbuf_PRECISION[0], l->s_PRECISION.op.translation_table, l, threading );
  }
  else{
    interpolate3_PRECISION( out, in, l, threading );
  }
}


// apply the restriction
void apply_R_PRECISION( vector_PRECISION out, vector_PRECISION in, level_struct* l, struct Thread *threading ){
  if( l->idle ) return;   // l is the fine level; idle processes do not take part
  if( l->depth==0 ){
    trans_PRECISION( l->sbuf_PRECISION[0], (vector_double)in, l->s_PRECISION.op.translation_table, l, threading );
    restrict_PRECISION( out, l->sbuf_PRECISION[0], l, threading );
  }
  else{
    restrict_PRECISION( out, in, l, threading );
  }
}





























// out = Gamma in on the finest grid, Gamma chosen by g.gamma_idx
// (0..3: gamma_mu, 5: gamma_5, anything else: identity). out != in.
void apply_gamma_idx_PRECISION( vector_PRECISION out, vector_PRECISION in, level_struct *l, struct Thread *threading ){
  if(g.gamma_idx == 5)      gamma5_PRECISION( out, in, l, threading );
  else if(g.gamma_idx == 0) gamma0_PRECISION( out, in, l, threading );
  else if(g.gamma_idx == 1) gamma1_PRECISION( out, in, l, threading );
  else if(g.gamma_idx == 2) gamma2_PRECISION( out, in, l, threading );
  else if(g.gamma_idx == 3) gamma3_PRECISION( out, in, l, threading );
  else {
    int start, end;
    compute_core_start_end( 0, l->inner_vector_size, &start, &end, l, threading );
    vector_PRECISION_copy( out, in, start, end, l );
  }
  SYNC_CORES(threading)
}

// out = M_l in, with M_l = R_hat_l Gamma P_t P_hat_l on the grid of level l:
// prolongate in to the finest grid, keep the finest timeslice g.time_slice,
// apply Gamma (g.gamma_idx), restrict back to level l. No solves.
// in and out live on level l; uses h->mlmc_b1 and h->mlmc_b2 as buffers
// (out must be neither of them).
void apply_M_timeslice_PRECISION( vector_PRECISION out, vector_PRECISION in, level_struct *l,
                                  hutchinson_PRECISION_struct* h, struct Thread *threading ){
  level_struct *finest_l = h->finest_level;
  level_struct *levels[g.num_levels];
  int nlev = 0, l_index = -1, start, end;

  for( level_struct *lt = finest_l; lt != NULL; lt = lt->next_level ){
    if( lt == l ) l_index = nlev;
    levels[nlev++] = lt;
  }
  if( l_index < 0 ) error0("apply_M_timeslice: level not found\n");

  // prolongate from level l to the finest grid
  compute_core_start_end( 0, l->inner_vector_size, &start, &end, l, threading );
  vector_PRECISION_copy( h->mlmc_b1, in, start, end, l );
  SYNC_CORES(threading)
  for( int i = l_index; i > 0; i-- ){
    level_struct *fine = levels[i-1];
    apply_P_PRECISION( h->mlmc_b2, h->mlmc_b1, fine, threading );
    compute_core_start_end( 0, fine->inner_vector_size, &start, &end, fine, threading );
    vector_PRECISION_copy( h->mlmc_b1, h->mlmc_b2, start, end, fine );
    SYNC_CORES(threading)
  }

  // P_t and Gamma on the finest grid (both site-local, so they commute)
  vector_PRECISION_timeslice( h->mlmc_b1, g.time_slice, finest_l, threading );
  apply_gamma_idx_PRECISION( h->mlmc_b2, h->mlmc_b1, finest_l, threading );

  // restrict back to level l
  for( int i = 0; i < l_index; i++ ){
    level_struct *fine = levels[i], *coarse = levels[i+1];
    apply_R_PRECISION( h->mlmc_b1, h->mlmc_b2, fine, threading );
    compute_core_start_end( 0, coarse->inner_vector_size, &start, &end, coarse, threading );
    vector_PRECISION_copy( h->mlmc_b2, h->mlmc_b1, start, end, coarse );
    SYNC_CORES(threading)
  }

  compute_core_start_end( 0, l->inner_vector_size, &start, &end, l, threading );
  vector_PRECISION_copy( out, h->mlmc_b2, start, end, l );
  SYNC_CORES(threading)
}

// 3D MGMLMC, level-l term (not the coarsest), on the grid of level l:
//   tr[ (A_l^-1 - P_l A_{l+1}^-1 R_l) M_l ],  M_l = R_hat_l Gamma P_t P_hat_l,
// which equals tr[ (P_hat_l A_l^-1 R_hat_l - P_hat_{l+1} A_{l+1}^-1 R_hat_{l+1}) Gamma P_t ]
// by cyclicity. M_l is site-local and lives on the coarse timeslice t_c that contains
// g.time_slice, so z is projected onto t_c of level l. Gamma from "gamma index:".
complex_PRECISION gamma_3D_mlmc_difference_PRECISION( int type_appl, level_struct *l, hutchinson_PRECISION_struct* h, struct Thread *threading ){
  int start, end;
  gmres_PRECISION_struct* p = get_p_struct_PRECISION( l );
  vector_PRECISION v = h->rademacher_vector;   // z (times probing vector) on level l
  complex_PRECISION aux = 0;

  if ( type_appl != -1 ) error0("gamma_3D_mlmc_difference: only type_appl = -1 is supported\n");

  vector_PRECISION_timeslice( v, level_time_slice( l->depth ), l, threading );

  // b = M_l v
  apply_M_timeslice_PRECISION( h->mlmc_testing, v, l, h, threading );

  // first term: x1 = A_l^-1 b, kept in h->mlmc_b1
  compute_core_start_end( 0, l->inner_vector_size, &start, &end, l, threading );
  vector_PRECISION_copy( p->b, h->mlmc_testing, start, end, l );
  apply_solver_PRECISION( l, threading );
  vector_PRECISION_copy( h->mlmc_b1, p->x, start, end, l );
  SYNC_CORES(threading)

  // second term: x2 = P_l A_{l+1}^-1 R_l b, in h->mlmc_b2
  apply_R_PRECISION( l->next_level->p_PRECISION.b, h->mlmc_testing, l, threading );
  apply_solver_PRECISION( l->next_level, threading );
  apply_P_PRECISION( h->mlmc_b2, l->next_level->p_PRECISION.x, l, threading );

  // v^dagger (x1 - x2) on level l
  vector_PRECISION_minus( h->mlmc_b1, h->mlmc_b1, h->mlmc_b2, start, end, l );
  SYNC_CORES(threading)
  aux = trace_inner_product_PRECISION( v, h->mlmc_b1, 0, l->inner_vector_size, l, threading );

  return aux;
}

// 3D MGMLMC, coarsest-level term: tr[ A_L^-1 M_L ] on the grid of level L.
complex_PRECISION gamma_3D_mlmc_coarsest_PRECISION( int type_appl, level_struct *l, hutchinson_PRECISION_struct* h, struct Thread *threading ){
  int start, end;
  gmres_PRECISION_struct* p = get_p_struct_PRECISION( l );
  vector_PRECISION v = h->rademacher_vector;   // z (times probing vector) on level L

  if ( type_appl != -1 ) error0("gamma_3D_mlmc_coarsest: only type_appl = -1 is supported\n");

  vector_PRECISION_timeslice( v, level_time_slice( l->depth ), l, threading );

  apply_M_timeslice_PRECISION( h->mlmc_testing, v, l, h, threading );

  compute_core_start_end( 0, l->inner_vector_size, &start, &end, l, threading );
  vector_PRECISION_copy( p->b, h->mlmc_testing, start, end, l );
  SYNC_CORES(threading)
  apply_solver_PRECISION( l, threading );

  return trace_inner_product_PRECISION( v, p->x, 0, l->inner_vector_size, l, threading );
}


complex_PRECISION g5_3D_connected_mlmc_difference_PRECISION( int type_appl, level_struct *l, hutchinson_PRECISION_struct* h, struct Thread *threading ){
  // store from fine to coarsest lvl in an array
  level_struct *finest_l = h->finest_level;
  level_struct *levels[g.num_levels];
  (void)levels;
  int lvl_nr = 0, l_index = -1;  (void)l_index;

  complex_PRECISION aux = 0;

  if(g.my_rank == 0) {
    printf("\n\n------------------function at level %d ------------------\n\n", l->depth);fflush(0);
  }

  for (level_struct *l_tmp = finest_l; l_tmp != NULL; l_tmp = l_tmp->next_level) {
    levels[lvl_nr] = l_tmp;
    if(g.my_rank == 0) printf("Storing depth %d \t Function called at depth %d\n", l_tmp->depth, l->depth);
    if (l_tmp == l) l_index = lvl_nr;
    lvl_nr++;
  }

  // Preallocate first and second term vectors at the finest
  vector_PRECISION first_term = NULL;
  PUBLIC_MALLOC( first_term, complex_PRECISION, finest_l->inner_vector_size );
  vector_PRECISION second_term = NULL;
  PUBLIC_MALLOC( second_term, complex_PRECISION, finest_l->inner_vector_size );

  // CONNECTED STEP #1. apply \Pi_{t+t'}
  // CONNECTED STEP #2. apply G5

  {
    int start, end;
    //gmres_PRECISION_struct* p = get_p_struct_PRECISION( l );
    compute_core_start_end( 0, finest_l->inner_vector_size, &start, &end, finest_l, threading );

    // Project rademacher vector into t+t'
    int bufft = g.time_slice;
    // TODO : check : is this assuming periodic in time ?
    g.time_slice = g.time_slice + g.time_slice_inner_connected;
    g.time_slice = g.time_slice%g.global_lattice[0][0];
    vector_PRECISION_ghg( h->rademacher_vector, 0, finest_l->inner_vector_size, finest_l );
    g.time_slice = bufft;
    //vector_PRECISION_ghg(  h->rademacher_vector, 0, finest_l->inner_vector_size, finest_l );

    // Apply gamma_5 AT FINEST level and save in buffer
    gamma5_PRECISION( h->mlmc_testing, h->rademacher_vector, finest_l, threading );
  }

  // CONNECTED STEP #3. apply level difference operator

  // LEVEL DIFFERENCE : FIRST TERM : result stored in first_term
  {
    int start, end;
    gmres_PRECISION_struct* p = get_p_struct_PRECISION( l );
    compute_core_start_end( 0, finest_l->inner_vector_size, &start, &end, finest_l, threading );

    // copy buffer to be restricted below
    if ( type_appl==-1 ) {
      vector_PRECISION_copy( h->mlmc_b1, h->mlmc_testing, start, end, finest_l );
    } else {
    // vector_PRECISION_copy( p->b, l->powerit_PRECISION.vecs[type_appl], start, end, l );
    }

    // Restrict rhs from finest to current level
    int start_tmp, end_tmp, start_tmp_c, end_tmp_c;

    for (int i = 0; i < l_index; ++i) {
      level_struct *fine   = levels[i];
      level_struct *coarse = levels[i + 1];

      compute_core_start_end(0, fine->inner_vector_size, &start_tmp, &end_tmp, fine, threading);
      compute_core_start_end(0, coarse->inner_vector_size, &start_tmp_c, &end_tmp_c, coarse, threading);
      // restrict from fine to coarse
      apply_R_PRECISION(h->mlmc_b2, h->mlmc_b1, fine, threading);
      vector_PRECISION_copy(h->mlmc_b1, h->mlmc_b2, start_tmp_c, end_tmp_c, coarse);

      if (g.my_rank == 0)
        printf("TERM 1: Restricting from depth %d to %d, \tfunction called at depth %d \n", fine->depth, coarse->depth, l->depth);
    }

    // Solve at current level
    compute_core_start_end( 0, l->inner_vector_size, &start, &end, l, threading );
    vector_PRECISION_copy( p->b, h->mlmc_b1, start, end, l);
    apply_solver_PRECISION( l, threading );

    // Prolongate solution to finest level
    vector_PRECISION_copy(h->mlmc_b1, p->x, start, end, l);
    for (int i = l_index; i > 0; --i) {
      level_struct *coarse = levels[i];
      level_struct *fine   = levels[i - 1];  /* finer level (towards finest) */

      compute_core_start_end(0, fine->inner_vector_size, &start_tmp, &end_tmp, fine, threading);
      apply_P_PRECISION(h->mlmc_b2, h->mlmc_b1, fine, threading);  /* prolongate from coarse to fine */
      vector_PRECISION_copy(h->mlmc_b1, h->mlmc_b2, start_tmp, end_tmp, fine);  /* copy result to fine level */

      if(g.my_rank == 0) {
        printf("TERM 1: Prolongating from depth %d to %d, \tfunction called at depth %d\n\n", coarse->depth, fine->depth, l->depth);
      }
    }

    // Copy result in first term vector at finest
    compute_core_start_end(0, finest_l->inner_vector_size, &start_tmp, &end_tmp, finest_l, threading);
    // copy for dot product
    vector_PRECISION_copy(first_term, h->mlmc_b1, start_tmp, end_tmp, finest_l);
  }
  // LEVEL DIFFERENCE : SECOND TERM : result stored in h->mlmc_b2
  // 1. Restrict (iteratively from finest)
  // 2. invert
  // 3. Prolongate (iteratively to finest)
  {
    int start, end, start_tmp, end_tmp, start_tmp_c, end_tmp_c;
    compute_core_start_end( 0, finest_l->inner_vector_size, &start, &end, finest_l, threading );

    if ( type_appl==-1 ) {
      // Copy gamma_5 x at the finest
      vector_PRECISION_copy( h->mlmc_b1, h->mlmc_testing, start, end, finest_l );

      // Restrict it from finest to NEXT level
      for (int i = 0; i < l_index+1; ++i) {
        level_struct *fine   = levels[i];
        level_struct *coarse = levels[i + 1];

        compute_core_start_end(0, fine->inner_vector_size, &start_tmp, &end_tmp, fine, threading);
        compute_core_start_end(0, coarse->inner_vector_size, &start_tmp_c, &end_tmp_c, coarse, threading);
        // restrict from fine to coarse
        apply_R_PRECISION(h->mlmc_b2, h->mlmc_b1, fine, threading);
        vector_PRECISION_copy(h->mlmc_b1, h->mlmc_b2, start_tmp_c, end_tmp_c, coarse);

        if (g.my_rank == 0)
          printf("TERM 2: Restricting from depth %d to %d, \tfunction called at depth %d \n", fine->depth, coarse->depth, l->depth);
      }

      // Copy restricted vector to rhs at the next level
      compute_core_start_end( 0, l->next_level->inner_vector_size, &start_tmp, &end_tmp, l->next_level, threading );
      vector_PRECISION_copy(l->next_level->p_PRECISION.b, h->mlmc_b1, start_tmp, end_tmp, l->next_level);
    } else {
    // apply_R_PRECISION( l->next_level->p_PRECISION.b, l->powerit_PRECISION.vecs[type_appl], l, threading );
    }

    // the input of this solve is l->next_level->p_PRECISION.x, the output l->next_level->p_PRECISION.b
    apply_solver_PRECISION( l->next_level, threading );

    // Prolongate Solution from NEXT level to finest
    vector_PRECISION_copy(h->mlmc_b1, l->next_level->p_PRECISION.x, start_tmp, end_tmp, l->next_level);
    for (int i = l_index+1; i > 0; --i) {
      level_struct *coarse = levels[i];
      // finer level (towards finest)
      level_struct *fine   = levels[i - 1];

      compute_core_start_end(0, fine->inner_vector_size, &start_tmp, &end_tmp, fine, threading);
      // prolongate from coarse to fine
      apply_P_PRECISION(h->mlmc_b2, h->mlmc_b1, fine, threading);
      // copy result to fine level
      vector_PRECISION_copy(h->mlmc_b1, h->mlmc_b2, start_tmp, end_tmp, fine);

      if(g.my_rank == 0) {
        printf("TERM 2: Prolongating from depth %d to %d, \tfunction called at depth %d\n\n", coarse->depth, fine->depth, l->depth);
      }
    }

    // Copy result in second term vector
    compute_core_start_end(0, finest_l->inner_vector_size, &start_tmp, &end_tmp, finest_l, threading);
    //copy for dot product
    vector_PRECISION_copy(second_term, h->mlmc_b1, start_tmp, end_tmp, finest_l);
  }
  // subtract the results
  {
    int start, end;
    //gmres_PRECISION_struct* p = get_p_struct_PRECISION( finest_l);
    compute_core_start_end( 0, finest_l->inner_vector_size, &start, &end, finest_l, threading );
    vector_PRECISION_minus( h->mlmc_b1, first_term, second_term, start, end, finest_l);
  }

  // CONNECTED STEP #4. apply \Pi_{t'}
  // CONNECTED STEP #5. apply G5

  {
    int start, end;
    //gmres_PRECISION_struct* p = get_p_struct_PRECISION( l );
    compute_core_start_end( 0, finest_l->inner_vector_size, &start, &end, finest_l, threading );

    // apply \Pi_{t'}
    //vector_PRECISION_ghg(  p->x, 0, l->inner_vector_size, l );
    int bufft = g.time_slice;
    // TODO : check : is this assuming periodic in time ?
    g.time_slice = g.time_slice_inner_connected;
    //g.time_slice = g.time_slice%g.global_lattice[0][0];
    //vector_PRECISION_ghg( h->rademacher_vector, 0, l->inner_vector_size, l );
    vector_PRECISION_ghg( h->mlmc_b1, 0, finest_l->inner_vector_size, finest_l );
    g.time_slice = bufft;

    // Apply gamma_5 AT FINEST level and save in buffer
    gamma5_PRECISION( h->mlmc_testing, h->mlmc_b1, finest_l, threading );
  }

  // CONNECTED STEP #6. apply level difference operator

  // LEVEL DIFFERENCE : FIRST TERM : result stored in first_term
  {
    int start, end;
    gmres_PRECISION_struct* p = get_p_struct_PRECISION( l );
    compute_core_start_end( 0, finest_l->inner_vector_size, &start, &end, finest_l, threading );

    // copy buffer to be restricted below
    if ( type_appl==-1 ) {
      vector_PRECISION_copy( h->mlmc_b1, h->mlmc_testing, start, end, finest_l );
    } else {
    // vector_PRECISION_copy( p->b, l->powerit_PRECISION.vecs[type_appl], start, end, l );
    }

    // Restrict rhs from finest to current level
    int start_tmp, end_tmp, start_tmp_c, end_tmp_c;

    for (int i = 0; i < l_index; ++i) {
      level_struct *fine   = levels[i];
      level_struct *coarse = levels[i + 1];

      compute_core_start_end(0, fine->inner_vector_size, &start_tmp, &end_tmp, fine, threading);
      compute_core_start_end(0, coarse->inner_vector_size, &start_tmp_c, &end_tmp_c, coarse, threading);
      // restrict from fine to coarse
      apply_R_PRECISION(h->mlmc_b2, h->mlmc_b1, fine, threading);
      vector_PRECISION_copy(h->mlmc_b1, h->mlmc_b2, start_tmp_c, end_tmp_c, coarse);

      if (g.my_rank == 0)
        printf("TERM 1: Restricting from depth %d to %d, \tfunction called at depth %d \n", fine->depth, coarse->depth, l->depth);
    }

    // Solve at current level
    compute_core_start_end( 0, l->inner_vector_size, &start, &end, l, threading );
    vector_PRECISION_copy( p->b, h->mlmc_b1, start, end, l);
    apply_solver_PRECISION( l, threading );

    // Prolongate solution to finest level

    vector_PRECISION_copy(h->mlmc_b1, p->x, start, end, l);
    for (int i = l_index; i > 0; --i) {
      level_struct *coarse = levels[i];
      level_struct *fine   = levels[i - 1];  /* finer level (towards finest) */

      compute_core_start_end(0, fine->inner_vector_size, &start_tmp, &end_tmp, fine, threading);
      apply_P_PRECISION(h->mlmc_b2, h->mlmc_b1, fine, threading);  /* prolongate from coarse to fine */
      vector_PRECISION_copy(h->mlmc_b1, h->mlmc_b2, start_tmp, end_tmp, fine);  /* copy result to fine level */

      if(g.my_rank == 0) {
        printf("TERM 1: Prolongating from depth %d to %d, \tfunction called at depth %d\n\n", coarse->depth, fine->depth, l->depth);
      }

    }

    // Copy result in first term vector at finest
    compute_core_start_end(0, finest_l->inner_vector_size, &start_tmp, &end_tmp, finest_l, threading);
    // copy for dot product
    vector_PRECISION_copy(first_term, h->mlmc_b1, start_tmp, end_tmp, finest_l);
  }
  // LEVEL DIFFERENCE : SECOND TERM : result stored in h->mlmc_b2
  // 1. Restrict (iteratively from finest)
  // 2. invert
  // 3. Prolongate (iteratively to finest)
  {
    int start, end, start_tmp, end_tmp, start_tmp_c, end_tmp_c;
    compute_core_start_end( 0, finest_l->inner_vector_size, &start, &end, finest_l, threading );

    if ( type_appl==-1 ) {
      // Copy gamma_5 x at the finest
      vector_PRECISION_copy( h->mlmc_b1, h->mlmc_testing, start, end, finest_l );

      // Restrict it from finest to NEXT level
      for (int i = 0; i < l_index+1; ++i) {
        level_struct *fine   = levels[i];
        level_struct *coarse = levels[i + 1];

        compute_core_start_end(0, fine->inner_vector_size, &start_tmp, &end_tmp, fine, threading);
        compute_core_start_end(0, coarse->inner_vector_size, &start_tmp_c, &end_tmp_c, coarse, threading);
        // restrict from fine to coarse
        apply_R_PRECISION(h->mlmc_b2, h->mlmc_b1, fine, threading);
        vector_PRECISION_copy(h->mlmc_b1, h->mlmc_b2, start_tmp_c, end_tmp_c, coarse);

        if (g.my_rank == 0)
          printf("TERM 2: Restricting from depth %d to %d, \tfunction called at depth %d \n", fine->depth, coarse->depth, l->depth);
      }

      // Copy restricted vector to rhs at the next level
      compute_core_start_end( 0, l->next_level->inner_vector_size, &start_tmp, &end_tmp, l->next_level, threading );
      vector_PRECISION_copy(l->next_level->p_PRECISION.b, h->mlmc_b1, start_tmp, end_tmp, l->next_level);
    } else {
    // apply_R_PRECISION( l->next_level->p_PRECISION.b, l->powerit_PRECISION.vecs[type_appl], l, threading );
    }

    // the input of this solve is l->next_level->p_PRECISION.x, the output l->next_level->p_PRECISION.b
    apply_solver_PRECISION( l->next_level, threading );

    // Prolongate Solution from NEXT level to finest

    vector_PRECISION_copy(h->mlmc_b1, l->next_level->p_PRECISION.x, start_tmp, end_tmp, l->next_level);
    for (int i = l_index+1; i > 0; --i) {
      level_struct *coarse = levels[i];
      // finer level (towards finest)
      level_struct *fine   = levels[i - 1];

      compute_core_start_end(0, fine->inner_vector_size, &start_tmp, &end_tmp, fine, threading);
      // prolongate from coarse to fine
      apply_P_PRECISION(h->mlmc_b2, h->mlmc_b1, fine, threading);
      // copy result to fine level
      vector_PRECISION_copy(h->mlmc_b1, h->mlmc_b2, start_tmp, end_tmp, fine);

      if(g.my_rank == 0) {
        printf("TERM 2: Prolongating from depth %d to %d, \tfunction called at depth %d\n\n", coarse->depth, fine->depth, l->depth);
      }
    }

    // Copy result in second term vector
    compute_core_start_end(0, finest_l->inner_vector_size, &start_tmp, &end_tmp, finest_l, threading);
    //copy for dot product
    vector_PRECISION_copy(second_term, h->mlmc_b1, start_tmp, end_tmp, finest_l);
  }
  // subtract the results
  {
    int start, end;
    //gmres_PRECISION_struct* p = get_p_struct_PRECISION( finest_l);
    compute_core_start_end( 0, finest_l->inner_vector_size, &start, &end, finest_l, threading );
    vector_PRECISION_minus( h->mlmc_b1, first_term, second_term, start, end, finest_l);
  }

  // do the dot product AT FINEST
  {
    gmres_PRECISION_struct* p = get_p_struct_PRECISION( finest_l);

    if ( type_appl==-1 ) {
      //if(g.trace_deflation_type[l->depth] != 0){
        //hutchinson_deflate_vector_PRECISION(h->mlmc_b1, l, threading);
      //}
      aux = trace_inner_product_PRECISION( h->rademacher_vector, h->mlmc_b1, 0, finest_l->inner_vector_size, finest_l, threading );
    } else {
      //aux = global_inner_product_PRECISION( l->powerit_PRECISION.vecs[type_appl], h->mlmc_b1, p->v_start, p->v_end, l, threading );
    }
  }

  PUBLIC_FREE( first_term, complex_PRECISION, finest_l->inner_vector_size );
  PUBLIC_FREE( second_term, complex_PRECISION, finest_l->inner_vector_size );

  return aux;
}




// complex_PRECISION g5_3D_hutchinson_mlmc_coarsest_PRECISION( int type_appl, level_struct *l, hutchinson_PRECISION_struct* h, struct Thread *threading ){
//   /* Store from fine to coarsest lvl in an array*/
//   level_struct *finest_l = h->finest_level;
//   level_struct *levels[g.num_levels];
//   int lvl_nr = 0, l_index = -1;  (void)l_index;

//   if(g.my_rank == 0)
//     printf("\n\n------------------function at COARSEST%d ------------------\n\n", l->depth);fflush(0);

//   for (level_struct *l_tmp = finest_l; l_tmp != NULL; l_tmp = l_tmp->next_level) {
//     levels[lvl_nr] = l_tmp;
//     if(g.my_rank == 0) printf("Storing depth %d \t Function called at depth %d\n", l_tmp->depth, l->depth);
//     if (l_tmp == l) l_index = lvl_nr;
//     lvl_nr++;
//   }

//  // Preallocate coarse term vector at the finest
//   vector_PRECISION *only_term;
//   PUBLIC_MALLOC( only_term, complex_PRECISION, finest_l->inner_vector_size );

//   // FIRST TERM : result stored in only_term
//   {
//     int start, end;
//     //gmres_PRECISION_struct* p = get_p_struct_PRECISION( finest_l );
//     compute_core_start_end( 0, finest_l->inner_vector_size, &start, &end, finest_l, threading );

//     /* Apply gamma_5 AT FINEST level and save in buffer */
//     gamma5_PRECISION( h->mlmc_testing, h->rademacher_vector, finest_l, threading );

//     /* copy buffer to be restricted below */
//     if ( type_appl==-1 ) {
//       vector_PRECISION_copy( h->mlmc_b1, h->mlmc_testing, start, end, finest_l );
//     } else {
//      // vector_PRECISION_copy( p->b, l->powerit_PRECISION.vecs[type_appl], start, end, l );
//     }

//     /* Restrict rhs from finest to current level */
//     int start_tmp, end_tmp;

//     for (int i = 0; i < l_index; ++i) {
//       level_struct *fine   = levels[i];
//       level_struct *coarse = levels[i + 1];

//       compute_core_start_end(0, fine->inner_vector_size, &start_tmp, &end_tmp, fine, threading);
//       apply_R_PRECISION(h->mlmc_b2, h->mlmc_b1, fine, threading);  // restrict from fine to coarse
//       vector_PRECISION_copy(h->mlmc_b1, h->mlmc_b2, start_tmp, end_tmp, coarse);

//       if (g.my_rank == 0)
//         printf("Coarsest 1: Restricting from depth %d to %d, \tfunction called at depth %d \n", fine->depth, coarse->depth, l->depth);
//     }
//     /* Solve at current level */
//     compute_core_start_end( 0, l->inner_vector_size, &start, &end, l, threading );
//     gmres_PRECISION_struct* p = get_p_struct_PRECISION( l );
//     vector_PRECISION_copy( p->b, h->mlmc_b1, start, end, l);

//     apply_solver_PRECISION( l, threading );

//     /* Prolongate solution to finest level */
//     vector_PRECISION_copy(h->mlmc_b1, p->x, start, end, l);

//     for (int i = l_index; i > 0; --i) {
//       level_struct *coarse = levels[i];
//       level_struct *fine   = levels[i - 1];  /* finer level (towards finest) */

//       compute_core_start_end(0, fine->inner_vector_size, &start_tmp, &end_tmp, fine, threading);
//       apply_P_PRECISION(h->mlmc_b2, h->mlmc_b1, fine, threading);  /* prolongate from coarse to fine */
//       vector_PRECISION_copy(h->mlmc_b1, h->mlmc_b2, start_tmp, end_tmp, fine);  /* copy result to fine level */

//       if(g.my_rank == 0)
//           printf("Coarsest: Prolongating from depth %d to %d, \tfunction called at depth %d\n\n", coarse->depth, fine->depth, l->depth);

//     }

//     // Copy result in the term vector
//     compute_core_start_end(0, finest_l->inner_vector_size, &start_tmp, &end_tmp, finest_l, threading);
//     vector_PRECISION_copy(only_term, h->mlmc_b1, start_tmp, end_tmp, finest_l);  /* copy for dot product */

//   }

//   //  perform dot product with rademacher AT FINEST
//   {
//     int start, end;
//     complex_PRECISION aux;
//     gmres_PRECISION_struct* p = get_p_struct_PRECISION( finest_l );

//     if ( type_appl==-1 ) {
//       //if(g.trace_deflation_type[l->depth] != 0){
//         //hutchinson_deflate_vector_PRECISION(h->mlmc_b1, l, threading);
//       //}
//       aux = global_inner_product_PRECISION( h->rademacher_vector, only_term, p->v_start, p->v_end, finest_l, threading );
//     } else {
//       //aux = global_inner_product_PRECISION( l->powerit_PRECISION.vecs[type_appl], h->mlmc_b1, p->v_start, p->v_end, l, threading );
//     }

//     return aux;
//   }

// }


complex_PRECISION gamma_3D_mlmc_hutchinson_driver_PRECISION( level_struct *l, struct Thread *threading ){
  int i;
  complex_PRECISION trace = 0.0;
  struct sample estimate;
  hutchinson_PRECISION_struct* h = &(l->h_PRECISION);
  level_struct* lx;

  // for all but coarsest level
  lx = l;
  for( i=0; i<g.num_levels-1; i++ ){
    // set the pointer to the mlmc difference operator
    h->hutch_compute_one_sample = gamma_3D_mlmc_difference_PRECISION;
    // colors and dilution are looped inside the blind routine, with one shared z0 per sample
    estimate = probing_hutchinson_blind_PRECISION( lx, h, 0, i, threading );
    trace += estimate.acc_trace / estimate.sample_size;
    { char label[64]; snprintf( label, 64, "level %d - level %d", i+1, i+2 ); report_piece_PRECISION( label, estimate, 1.0, threading ); }
    
    //if deflation vectors are available
    //if(g.trace_deflation_type[lx->depth] != 0){
    //  trace += hutchinson_deflated_direct_term_PRECISION( lx, h, threading );
    //}
    lx = lx->next_level;
  }

  // coarsest level
  // set the pointer to the coarsest-level Hutchinson estimator
  h->hutch_compute_one_sample = gamma_3D_mlmc_coarsest_PRECISION;
  
  // colors and dilution are looped inside the blind routine, with one shared z0 per sample
    estimate = probing_hutchinson_blind_PRECISION( lx, h, 0, i, threading );
    trace += estimate.acc_trace / estimate.sample_size;
    { char label[64]; snprintf( label, 64, "coarsest (level %d)", g.num_levels ); report_piece_PRECISION( label, estimate, 1.0, threading ); }
  //if deflation vectors are available
  //if(g.trace_deflation_type[lx->depth] != 0){
  //  trace += hutchinson_deflated_direct_term_PRECISION( lx, h, threading );
  //}

  return trace;
}

complex_PRECISION g5_3D_connected_mlmc_driver_PRECISION( level_struct *l, struct Thread *threading ){
  int i, j;
  complex_PRECISION trace = 0.0;
  complex_PRECISION global_connected_trace = 0.0;
  struct sample estimate;
  hutchinson_PRECISION_struct* h = &(l->h_PRECISION);

  h->lx_i = h->finest_level;
  h->lx_j = h->finest_level;

  for(i = 0; i < g.num_levels; i++){
    h->lx_j = h->finest_level;
    for(j = 0; j < g.num_levels; j++){
    
      if(g.my_rank == 0) printf("\n Computing trace for G_{%d, %d} \n", i,j);
      
      trace = 0.0;
      h->hutch_compute_one_sample = connected_outer_PRECISION;
      
      // Call to blind from t+t’ with t’ = 0, ..., T-1
      for ( g.time_slice_inner_connected=0; g.time_slice_inner_connected<g.global_lattice[0][0]; g.time_slice_inner_connected++ ) {
        // colors and dilution are looped inside the blind routine, with one shared z0 per sample
    estimate = probing_hutchinson_blind_PRECISION( h->finest_level, h, 0, i, threading );
    trace += estimate.acc_trace / estimate.sample_size;
    { char label[64]; snprintf( label, 64, "G_{%d,%d} (average over t')", h->lx_i->depth, h->lx_j->depth ); report_piece_PRECISION( label, estimate, 1.0/g.global_lattice[0][0], threading ); }
      }
      
      trace = trace/g.time_slice_inner_connected;
      global_connected_trace += trace;
      //int nlevs = g.num_levels;
      //int idx = h->lx_i->depth * nlevs + h->lx_j->depth;  
      if (g.my_rank == 0) printf("\n Trace of G_{%d, %d}(t=%d) = %f +%f\n", h->lx_i->depth, h->lx_j->depth, g.time_slice, CSPLIT(trace));
      
      if (j < g.num_levels - 1)    
        h->lx_j = h->lx_j->next_level;
    }
      if (i < g.num_levels - 1)
        h->lx_i = h->lx_i->next_level;
  }

  return global_connected_trace;
  
}


/* complex_PRECISION g5_3D_connected_mlmc_driver_PRECISION( level_struct *l, struct Thread *threading ){
  int i;
  complex_PRECISION trace = 0.0;
  struct sample estimate;
  hutchinson_PRECISION_struct* h = &(l->h_PRECISION);
  level_struct* lx;

  // for all but coarsest level
  lx = l;
  //for( i=0; i<1; i++ ) {
  for( i=0; i<g.num_levels-1; i++ ){
    // set the pointer to the mlmc difference operator
    h->hutch_compute_one_sample = g5_3D_connected_mlmc_difference_PRECISION;  
    
    for ( g.time_slice_inner_connected=0; g.time_slice_inner_connected<g.global_lattice[0][0]; g.time_slice_inner_connected++ ) {
      if (g.probing) {
        for (g.coloring_count = 1; g.coloring_count < g.num_colors[i] + 1; g.coloring_count++) {
          if(g.my_rank==0) {
            printf("\nEstimating trace at color %d\n", g.coloring_count);
          }

          estimate = hutchinson_blind_PRECISION(lx, h, 0, threading);
          trace += estimate.acc_trace / estimate.sample_size;
        }
      } else {
        estimate = hutchinson_blind_PRECISION(lx, h, 0, threading);
        trace += estimate.acc_trace / estimate.sample_size;
      }
    }

    //if deflation vectors are available
    //if(g.trace_deflation_type[lx->depth] != 0){
    //  trace += hutchinson_deflated_direct_term_PRECISION( lx, h, threading );
    //}
    lx = lx->next_level;
  }

//  // coarsest level
//  // set the pointer to the coarsest-level Hutchinson estimator
//  h->hutch_compute_one_sample = gamma_3D_mlmc_coarsest_PRECISION;

//  if (g.probing) {
//    for (g.coloring_count = 1; g.coloring_count < g.num_colors[i] + 1; g.coloring_count++) {
//      if(g.my_rank==0) {
//        printf("\nEstimating trace at color %d\n", g.coloring_count);
//      }

//      estimate = hutchinson_blind_PRECISION(lx, h, 0, threading);
//      trace += estimate.acc_trace / estimate.sample_size;
//    }
//  } else {
//      estimate = hutchinson_blind_PRECISION(lx, h, 0, threading);
//      trace += estimate.acc_trace / estimate.sample_size;
//  }
//  //if deflation vectors are available
//  //if(g.trace_deflation_type[lx->depth] != 0){
//  //  trace += hutchinson_deflated_direct_term_PRECISION( lx, h, threading );
//  //}

  return trace;
}
*/

complex_PRECISION g5_3D_connected_split_orthogonal_PRECISION( int type_appl, level_struct *l, hutchinson_PRECISION_struct* h, struct Thread *threading ){
  // store from fine to coarsest lvl in an array
  level_struct *finest_l = h->finest_level;
  level_struct *levels[g.num_levels];
  (void)levels;
  int lvl_nr = 0, l_index = -1;  (void)l_index;

  complex_PRECISION aux = 0;

  if(g.my_rank == 0) {
    printf("\n\n------------------function at level %d ------------------\n\n", l->depth);fflush(0);
  }

  for (level_struct *l_tmp = finest_l; l_tmp != NULL; l_tmp = l_tmp->next_level) {
    levels[lvl_nr] = l_tmp;
    if(g.my_rank == 0) printf("Storing depth %d \t Function called at depth %d\n", l_tmp->depth, l->depth);
    if (l_tmp == l) l_index = lvl_nr;
    lvl_nr++;
  }

  // CONNECTED STEP #1. apply \Pi_{t+t'}
  // CONNECTED STEP #2. apply G5

  {
    int start, end;
    //gmres_PRECISION_struct* p = get_p_struct_PRECISION( l );
    compute_core_start_end( 0, finest_l->inner_vector_size, &start, &end, finest_l, threading );

    // Project rademacher vector into t+t'
    int bufft = g.time_slice;
    // TODO : check : is this assuming periodic in time ?
    g.time_slice = g.time_slice + g.time_slice_inner_connected;
    g.time_slice = g.time_slice%g.global_lattice[0][0];
    vector_PRECISION_ghg( h->rademacher_vector, 0, finest_l->inner_vector_size, finest_l );
    g.time_slice = bufft;
    //vector_PRECISION_ghg(  h->rademacher_vector, 0, finest_l->inner_vector_size, finest_l );

    // Apply gamma_5 AT FINEST level and save in buffer
    gamma5_PRECISION( h->mlmc_testing, h->rademacher_vector, finest_l, threading );
  }

  // CONNECTED STEP #3. apply level orthogonal difference operator
  {
    int start, end;
    gmres_PRECISION_struct* p = get_p_struct_PRECISION( l );
    compute_core_start_end( 0, finest_l->inner_vector_size, &start, &end, finest_l, threading );

    // copy buffer to be restricted below
    if ( type_appl==-1 ) {
      vector_PRECISION_copy( h->mlmc_b1, h->mlmc_testing, start, end, finest_l );
    } else {
    // vector_PRECISION_copy( p->b, l->powerit_PRECISION.vecs[type_appl], start, end, l );
    }

    // Restrict rhs from finest to current level
    int start_tmp, end_tmp, start_tmp_c, end_tmp_c;

    for (int i = 0; i < l_index; ++i) {
      level_struct *fine   = levels[i];
      level_struct *coarse = levels[i + 1];

      compute_core_start_end(0, fine->inner_vector_size, &start_tmp, &end_tmp, fine, threading);
      compute_core_start_end(0, coarse->inner_vector_size, &start_tmp_c, &end_tmp_c, coarse, threading);
      // restrict from fine to coarse
      apply_R_PRECISION(h->mlmc_b2, h->mlmc_b1, fine, threading);
      vector_PRECISION_copy(h->mlmc_b1, h->mlmc_b2, start_tmp_c, end_tmp_c, coarse);

      if (g.my_rank == 0)
        printf("TERM 1: Restricting from depth %d to %d, \tfunction called at depth %d \n", fine->depth, coarse->depth, l->depth);
    }

    // Orthogonal difference (I - P P^H )(D_1^{-1} gamma_5 G_1^H) ) x
    {
      compute_core_start_end( 0, l->inner_vector_size, &start, &end, l, threading );
      vector_PRECISION_copy( p->b, h->mlmc_b1, start, end, l);
      apply_solver_PRECISION( l, threading );
      
      // P^H (D_1^{-1} gamma_5 G_1^H) x
      apply_R_PRECISION( h->mlmc_b1, p->x, l, threading );
      // P P^H (D_1^{-1} gamma_5 G_1^H) x
      apply_P_PRECISION( h->mlmc_b2, h->mlmc_b1, l, threading );
      // (I - P P^H )(D_1^{-1} gamma_5 G_1^H) ) x
      vector_PRECISION_minus( h->mlmc_b1,  p->x, h->mlmc_b2, start, end, l );
    }
    // Prolongate solution to finest level
    //vector_PRECISION_copy(h->mlmc_b1, p->x, start, end, l);
    for (int i = l_index; i > 0; --i) {
      level_struct *coarse = levels[i];
      level_struct *fine   = levels[i - 1];  /* finer level (towards finest) */

      compute_core_start_end(0, fine->inner_vector_size, &start_tmp, &end_tmp, fine, threading);
      apply_P_PRECISION(h->mlmc_b2, h->mlmc_b1, fine, threading);  /* prolongate from coarse to fine */
      vector_PRECISION_copy(h->mlmc_b1, h->mlmc_b2, start_tmp, end_tmp, fine);  /* copy result to fine level */

      if(g.my_rank == 0) {
        printf("TERM 1: Prolongating from depth %d to %d, \tfunction called at depth %d\n\n", coarse->depth, fine->depth, l->depth);
      }
    }
  }

  // CONNECTED STEP #4. apply \Pi_{t'}
  // CONNECTED STEP #5. apply G5
  {
    int start, end;
    //gmres_PRECISION_struct* p = get_p_struct_PRECISION( l );
    compute_core_start_end( 0, finest_l->inner_vector_size, &start, &end, finest_l, threading );

    // apply \Pi_{t'}
    //vector_PRECISION_ghg(  p->x, 0, l->inner_vector_size, l );
    int bufft = g.time_slice;
    // TODO : check : is this assuming periodic in time ?
    g.time_slice = g.time_slice_inner_connected;
    //g.time_slice = g.time_slice%g.global_lattice[0][0];
    //vector_PRECISION_ghg( h->rademacher_vector, 0, l->inner_vector_size, l );
    vector_PRECISION_ghg( h->mlmc_b1, 0, finest_l->inner_vector_size, finest_l );
    g.time_slice = bufft;

    // Apply gamma_5 AT FINEST level and save in buffer
    gamma5_PRECISION( h->mlmc_testing, h->mlmc_b1, finest_l, threading );
  }

  // CONNECTED STEP #6. apply level orthogonal difference operator
  {
    int start, end;
    gmres_PRECISION_struct* p = get_p_struct_PRECISION( l );
    compute_core_start_end( 0, finest_l->inner_vector_size, &start, &end, finest_l, threading );

    // copy buffer to be restricted below
    if ( type_appl==-1 ) {
      vector_PRECISION_copy( h->mlmc_b1, h->mlmc_testing, start, end, finest_l );
    } else {
    // vector_PRECISION_copy( p->b, l->powerit_PRECISION.vecs[type_appl], start, end, l );
    }

    // Restrict rhs from finest to current level
    int start_tmp, end_tmp, start_tmp_c, end_tmp_c;

    for (int i = 0; i < l_index; ++i) {
      level_struct *fine   = levels[i];
      level_struct *coarse = levels[i + 1];

      compute_core_start_end(0, fine->inner_vector_size, &start_tmp, &end_tmp, fine, threading);
      compute_core_start_end(0, coarse->inner_vector_size, &start_tmp_c, &end_tmp_c, coarse, threading);
      // restrict from fine to coarse
      apply_R_PRECISION(h->mlmc_b2, h->mlmc_b1, fine, threading);
      vector_PRECISION_copy(h->mlmc_b1, h->mlmc_b2, start_tmp_c, end_tmp_c, coarse);

      if (g.my_rank == 0)
        printf("TERM 1: Restricting from depth %d to %d, \tfunction called at depth %d \n", fine->depth, coarse->depth, l->depth);
    }

    // Orthogonal difference (I - P P^H )(D_1^{-1} gamma_5 G_1^H) ) x
    {
      compute_core_start_end( 0, l->inner_vector_size, &start, &end, l, threading );
      vector_PRECISION_copy( p->b, h->mlmc_b1, start, end, l);
      apply_solver_PRECISION( l, threading );
      
      // P^H (D_1^{-1} gamma_5 G_1^H) x
      apply_R_PRECISION( h->mlmc_b1, p->x, l, threading );
      // P P^H (D_1^{-1} gamma_5 G_1^H) x
      apply_P_PRECISION( h->mlmc_b2, h->mlmc_b1, l, threading );
      // (I - P P^H )(D_1^{-1} gamma_5 G_1^H) ) x
      vector_PRECISION_minus( h->mlmc_b1,  p->x, h->mlmc_b2, start, end, l );
    }
    // Prolongate solution to finest level
    //vector_PRECISION_copy(h->mlmc_b1, p->x, start, end, l);
    for (int i = l_index; i > 0; --i) {
      level_struct *coarse = levels[i];
      level_struct *fine   = levels[i - 1];  /* finer level (towards finest) */

      compute_core_start_end(0, fine->inner_vector_size, &start_tmp, &end_tmp, fine, threading);
      apply_P_PRECISION(h->mlmc_b2, h->mlmc_b1, fine, threading);  /* prolongate from coarse to fine */
      vector_PRECISION_copy(h->mlmc_b1, h->mlmc_b2, start_tmp, end_tmp, fine);  /* copy result to fine level */

      if(g.my_rank == 0) {
        printf("TERM 1: Prolongating from depth %d to %d, \tfunction called at depth %d\n\n", coarse->depth, fine->depth, l->depth);
      }
    }
  }
  
  
  // do the dot product AT FINEST
  {
    gmres_PRECISION_struct* p = get_p_struct_PRECISION( finest_l);

    if ( type_appl==-1 ) {
      //if(g.trace_deflation_type[l->depth] != 0){
        //hutchinson_deflate_vector_PRECISION(h->mlmc_b1, l, threading);
      //}
      aux = trace_inner_product_PRECISION( h->rademacher_vector, h->mlmc_b1, 0, finest_l->inner_vector_size, finest_l, threading );
    } else {
      //aux = global_inner_product_PRECISION( l->powerit_PRECISION.vecs[type_appl], h->mlmc_b1, p->v_start, p->v_end, l, threading );
    }
  }

  return aux;
}



complex_PRECISION g5_3D_connected_split_intermediate_PRECISION( int type_appl, level_struct *l, hutchinson_PRECISION_struct* h, struct Thread *threading ){
  // store from fine to coarsest lvl in an array
  level_struct *finest_l = h->finest_level;
  level_struct *second_l = h->finest_level->next_level;
  level_struct *levels[g.num_levels];
  (void)levels;
  int lvl_nr = 0, l_index = -1;  (void)l_index;

  complex_PRECISION aux = 0.0;

  if(g.my_rank == 0) {
    printf("\n\n------------------function at level %d ------------------\n\n", l->depth);fflush(0);
  }

  for (level_struct *l_tmp = finest_l; l_tmp != NULL; l_tmp = l_tmp->next_level) {
    levels[lvl_nr] = l_tmp;
    if(g.my_rank == 0) printf("Storing depth %d \t Function called at depth %d\n", l_tmp->depth, l->depth);
    if (l_tmp == l) l_index = lvl_nr;
    lvl_nr++;
  }

  // CONNECTED STEP #1. apply \Pi_{t'} at SECOND LEVEL
  {
    int start, end;
    //gmres_PRECISION_struct* p = get_p_struct_PRECISION( l );
    compute_core_start_end( 0, second_l->inner_vector_size, &start, &end, second_l, threading );

    int bufft = g.time_slice;
    // TODO : check : is this assuming periodic in time ?
    g.time_slice = g.time_slice_inner_connected;
    vector_PRECISION_ghg( h->rademacher_vector, 0, second_l->inner_vector_size, second_l );
    g.time_slice = bufft;
  }

  // CONNECTED STEP #2. apply level right difference operator (I -  PP^H) (D_1^{-1} P)  \Pi_{t'}  x
  {
    int start, end;
    compute_core_start_end( 0, second_l->inner_vector_size, &start, &end, second_l, threading );

    // copy buffer to be projected below
    if ( type_appl==-1 ) {
      vector_PRECISION_copy( h->mlmc_b1, h->rademacher_vector, start, end, second_l );
    } else {
    // vector_PRECISION_copy( p->b, l->powerit_PRECISION.vecs[type_appl], start, end, l );
    }

    // Right term solve (D_1^{-1} P)  \Pi_{t'}  x
    {
      compute_core_start_end( 0, l->inner_vector_size, &start, &end, l, threading );
    
      apply_P_PRECISION( h->mlmc_b2, h->mlmc_b1, l, threading );
      gmres_PRECISION_struct* p = get_p_struct_PRECISION( l );
      vector_PRECISION_copy( p->b, h->mlmc_b2, start, end, l);
      apply_solver_PRECISION( l, threading );           
    }
    
    // Intermediate difference (I -  PP^H) (D_1^{-1} P)  \Pi_{t'}  x
    {
      gmres_PRECISION_struct* p = get_p_struct_PRECISION( l );
      compute_core_start_end( 0, l->inner_vector_size, &start, &end, l, threading );
      apply_R_PRECISION( h->mlmc_b1, p->x, l, threading );
      apply_P_PRECISION( h->mlmc_b2, h->mlmc_b1, l, threading );
      vector_PRECISION_minus( h->mlmc_b1, p->x, h->mlmc_b2, start, end, l);

    }
    
  }
  // CONNECTED STEP #4. apply \Pi_{t+t'}     
  // CONNECTED STEP #5. apply Gamma_5 \Pi_{t+t') (I -  PP^H) (D_1^{-1} P)  \Pi_{t'}  x
  {
    int start, end;
    //gmres_PRECISION_struct* p = get_p_struct_PRECISION( l );
    compute_core_start_end( 0, finest_l->inner_vector_size, &start, &end, finest_l, threading );

    // Project rademacher vector into t+t'
    int bufft = g.time_slice;
    // TODO : check : is this assuming periodic in time ?
    g.time_slice = g.time_slice + g.time_slice_inner_connected;
    g.time_slice = g.time_slice%g.global_lattice[0][0];
    vector_PRECISION_ghg( h->mlmc_b1, 0, finest_l->inner_vector_size, finest_l );
    g.time_slice = bufft;

    // Apply gamma_5 AT FINEST level and save in buffer
    gamma5_PRECISION( h->mlmc_b1, h->mlmc_b1, finest_l, threading );
  }
  
  // CONNECTED STEP #6. apply left difference (P^H D^{-1} - D_c^{-1} P^H)
  {
    //second term: D_c^{-1} P^H
      {
      int start, end;
      compute_core_start_end( 0, finest_l->inner_vector_size, &start, &end, finest_l, threading );
      
      gmres_PRECISION_struct* p = get_p_struct_PRECISION( second_l );
      
      apply_R_PRECISION(p->b, h->mlmc_b1, l, threading );

      apply_solver_PRECISION( second_l, threading );       
            
      compute_core_start_end( 0, second_l->inner_vector_size, &start, &end, second_l, threading );
      vector_PRECISION_copy( h->mlmc_b2, p->x, start, end, second_l );
    }
    
    //First term P^H D^{-1}
    {
      int start, end;
      compute_core_start_end( 0, finest_l->inner_vector_size, &start, &end, finest_l, threading );
      gmres_PRECISION_struct* p = get_p_struct_PRECISION( finest_l);
      vector_PRECISION_copy( p->b, h->mlmc_b1, start, end, finest_l );

      apply_solver_PRECISION( finest_l, threading );      

      apply_R_PRECISION(h->mlmc_b1, p->x, l, threading );
      
    }
    
    // Diffrence 
    {
      int start, end;
      compute_core_start_end( 0, second_l->inner_vector_size, &start, &end, second_l, threading );
      vector_PRECISION_minus( h->mlmc_b1, h->mlmc_b1, h->mlmc_b2, start, end, second_l);
    }
    
    // Dot product
    { 
      int start, end;
      compute_core_start_end( 0, second_l->inner_vector_size, &start, &end, second_l, threading );
      aux = trace_inner_product_PRECISION( h->rademacher_vector, h->mlmc_b1, start, end, second_l, threading );
    }
  }
  
  return aux;
}






void connected_mlmc_PRECISION_non_difference( vector_PRECISION out, vector_PRECISION in, level_struct *l,  hutchinson_PRECISION_struct* h, struct Thread *threading ){

  // store from fine to coarsest lvl in an array
  level_struct *finest_l = h->finest_level;
  level_struct *levels[g.num_levels];
  int lvl_nr = 0, l_index = -1;  (void)l_index;

  if(g.my_rank == 0) {
    printf("\n\n------------------function at level %d ------------------\n\n", l->depth);fflush(0);
  }

  for (level_struct *l_tmp = finest_l; l_tmp != NULL; l_tmp = l_tmp->next_level) {
    levels[lvl_nr] = l_tmp;
    if(g.my_rank == 0) printf("Storing depth %d \t Function called at depth %d\n", l_tmp->depth, l->depth);
    if (l_tmp == l) l_index = lvl_nr;
    lvl_nr++;
  }
  
  // Apply Gamma_5 at the finest
  {    
    gamma5_PRECISION( h->mlmc_testing, in, finest_l, threading );
  }
  
  // Solve at current level : FIRST TERM : result stored in first_term
  {
    int start, end;
    gmres_PRECISION_struct* p = get_p_struct_PRECISION( l );
    compute_core_start_end( 0, finest_l->inner_vector_size, &start, &end, finest_l, threading );

    // copy buffer to be restricted below
    //if ( type_appl==-1 ) {
      vector_PRECISION_copy( h->mlmc_b1, h->mlmc_testing, start, end, finest_l );
    //} else {
    // // vector_PRECISION_copy( p->b, l->powerit_PRECISION.vecs[type_appl], start, end, l );
    // }

    // Restrict rhs from finest to current level
    int start_tmp, end_tmp, start_tmp_c, end_tmp_c;

    for (int i = 0; i < l_index; ++i) {
      level_struct *fine   = levels[i];
      level_struct *coarse = levels[i + 1];

      compute_core_start_end(0, fine->inner_vector_size, &start_tmp, &end_tmp, fine, threading);
      compute_core_start_end(0, coarse->inner_vector_size, &start_tmp_c, &end_tmp_c, coarse, threading);
      // restrict from fine to coarse
      apply_R_PRECISION(h->mlmc_b2, h->mlmc_b1, fine, threading);
      vector_PRECISION_copy(h->mlmc_b1, h->mlmc_b2, start_tmp_c, end_tmp_c, coarse);

      if (g.my_rank == 0)
        printf("TERM 1: Restricting from depth %d to %d, \tfunction called at depth %d \n", fine->depth, coarse->depth, l->depth);
    }

    // Solve at current level
    compute_core_start_end( 0, l->inner_vector_size, &start, &end, l, threading );
    vector_PRECISION_copy( p->b, h->mlmc_b1, start, end, l);
    apply_solver_PRECISION( l, threading );

    // Prolongate solution to finest level
    vector_PRECISION_copy(h->mlmc_b1, p->x, start, end, l);
    for (int i = l_index; i > 0; --i) {
      level_struct *coarse = levels[i];
      level_struct *fine   = levels[i - 1];  /* finer level (towards finest) */

      compute_core_start_end(0, fine->inner_vector_size, &start_tmp, &end_tmp, fine, threading);
      apply_P_PRECISION(h->mlmc_b2, h->mlmc_b1, fine, threading);  /* prolongate from coarse to fine */
      vector_PRECISION_copy(h->mlmc_b1, h->mlmc_b2, start_tmp, end_tmp, fine);  /* copy result to fine level */

      if(g.my_rank == 0) {
        printf("TERM 1: Prolongating from depth %d to %d, \tfunction called at depth %d\n\n", coarse->depth, fine->depth, l->depth);
      }
    }

    // Copy result in first term vector at finest
    compute_core_start_end(0, finest_l->inner_vector_size, &start_tmp, &end_tmp, finest_l, threading);
    // copy for dot product
    vector_PRECISION_copy(out, h->mlmc_b1, start_tmp, end_tmp, finest_l);
  }
  

}







void connected_mlmc_PRECISION_difference( vector_PRECISION out, vector_PRECISION in, level_struct *l,  hutchinson_PRECISION_struct* h, struct Thread *threading ){
  // store from fine to coarsest lvl in an array
  level_struct *finest_l = h->finest_level;
  level_struct *levels[g.num_levels];
  (void)levels;
  int lvl_nr = 0, l_index = -1;  (void)l_index;

  if(g.my_rank == 0) {
    printf("\n\n------------------function at level %d ------------------\n\n", l->depth);fflush(0);
  }

  for (level_struct *l_tmp = finest_l; l_tmp != NULL; l_tmp = l_tmp->next_level) {
    levels[lvl_nr] = l_tmp;
    if(g.my_rank == 0) printf("Storing depth %d \t Function called at depth %d\n", l_tmp->depth, l->depth);
    if (l_tmp == l) l_index = lvl_nr;
    lvl_nr++;
  }
  
  // Preallocate first and second term vectors at the finest
  vector_PRECISION first_term = NULL;
  PUBLIC_MALLOC( first_term, complex_PRECISION, finest_l->inner_vector_size );
  vector_PRECISION second_term = NULL;
  PUBLIC_MALLOC( second_term, complex_PRECISION, finest_l->inner_vector_size );
    
  // Apply Gamma_5 at the finest
  {    
    gamma5_PRECISION( h->mlmc_testing, in, finest_l, threading );
  }
  
  // LEVEL DIFFERENCE : FIRST TERM : result stored in first_term
  {
    int start, end;
    gmres_PRECISION_struct* p = get_p_struct_PRECISION( l );
    compute_core_start_end( 0, finest_l->inner_vector_size, &start, &end, finest_l, threading );

    // copy buffer to be restricted below
    //if ( type_appl==-1 ) {
      vector_PRECISION_copy( h->mlmc_b1, h->mlmc_testing, start, end, finest_l );
    //} else {
    // // vector_PRECISION_copy( p->b, l->powerit_PRECISION.vecs[type_appl], start, end, l );
    // }

    // Restrict rhs from finest to current level
    int start_tmp, end_tmp, start_tmp_c, end_tmp_c;

    for (int i = 0; i < l_index; ++i) {
      level_struct *fine   = levels[i];
      level_struct *coarse = levels[i + 1];

      compute_core_start_end(0, fine->inner_vector_size, &start_tmp, &end_tmp, fine, threading);
      compute_core_start_end(0, coarse->inner_vector_size, &start_tmp_c, &end_tmp_c, coarse, threading);
      // restrict from fine to coarse
      apply_R_PRECISION(h->mlmc_b2, h->mlmc_b1, fine, threading);
      vector_PRECISION_copy(h->mlmc_b1, h->mlmc_b2, start_tmp_c, end_tmp_c, coarse);

      if (g.my_rank == 0)
        printf("TERM 1: Restricting from depth %d to %d, \tfunction called at depth %d \n", fine->depth, coarse->depth, l->depth);
    }

    // Solve at current level
    compute_core_start_end( 0, l->inner_vector_size, &start, &end, l, threading );
    vector_PRECISION_copy( p->b, h->mlmc_b1, start, end, l);
    apply_solver_PRECISION( l, threading );

    // Prolongate solution to finest level
    vector_PRECISION_copy(h->mlmc_b1, p->x, start, end, l);
    for (int i = l_index; i > 0; --i) {
      level_struct *coarse = levels[i];
      level_struct *fine   = levels[i - 1];  /* finer level (towards finest) */

      compute_core_start_end(0, fine->inner_vector_size, &start_tmp, &end_tmp, fine, threading);
      apply_P_PRECISION(h->mlmc_b2, h->mlmc_b1, fine, threading);  /* prolongate from coarse to fine */
      vector_PRECISION_copy(h->mlmc_b1, h->mlmc_b2, start_tmp, end_tmp, fine);  /* copy result to fine level */

      if(g.my_rank == 0) {
        printf("TERM 1: Prolongating from depth %d to %d, \tfunction called at depth %d\n\n", coarse->depth, fine->depth, l->depth);
      }
    }

    // Copy result in first term vector at finest
    compute_core_start_end(0, finest_l->inner_vector_size, &start_tmp, &end_tmp, finest_l, threading);
    // copy for dot product
    vector_PRECISION_copy(first_term, h->mlmc_b1, start_tmp, end_tmp, finest_l);
  }
  
  // LEVEL DIFFERENCE : SECOND TERM : result stored in h->mlmc_b2
  // 1. Restrict (iteratively from finest)
  // 2. invert
  // 3. Prolongate (iteratively to finest)
  {
    int start, end, start_tmp, end_tmp, start_tmp_c, end_tmp_c;
    compute_core_start_end( 0, finest_l->inner_vector_size, &start, &end, finest_l, threading );

    //if ( type_appl==-1 ) {
      // Copy gamma_5 x at the finest
      vector_PRECISION_copy( h->mlmc_b1, h->mlmc_testing, start, end, finest_l );

      // Restrict it from finest to NEXT level
      for (int i = 0; i < l_index+1; ++i) {
        level_struct *fine   = levels[i];
        level_struct *coarse = levels[i + 1];

        compute_core_start_end(0, fine->inner_vector_size, &start_tmp, &end_tmp, fine, threading);
        compute_core_start_end(0, coarse->inner_vector_size, &start_tmp_c, &end_tmp_c, coarse, threading);
        // restrict from fine to coarse
        apply_R_PRECISION(h->mlmc_b2, h->mlmc_b1, fine, threading);
        vector_PRECISION_copy(h->mlmc_b1, h->mlmc_b2, start_tmp_c, end_tmp_c, coarse);

        if (g.my_rank == 0)
          printf("TERM 2: Restricting from depth %d to %d, \tfunction called at depth %d \n", fine->depth, coarse->depth, l->depth);
      }

      // Copy restricted vector to rhs at the next level
      compute_core_start_end( 0, l->next_level->inner_vector_size, &start_tmp, &end_tmp, l->next_level, threading );
      vector_PRECISION_copy(l->next_level->p_PRECISION.b, h->mlmc_b1, start_tmp, end_tmp, l->next_level);
  // } else {
    // apply_R_PRECISION( l->next_level->p_PRECISION.b, l->powerit_PRECISION.vecs[type_appl], l, threading );
  // }

    // the input of this solve is l->next_level->p_PRECISION.x, the output l->next_level->p_PRECISION.b
    apply_solver_PRECISION( l->next_level, threading );

    // Prolongate Solution from NEXT level to finest
    vector_PRECISION_copy(h->mlmc_b1, l->next_level->p_PRECISION.x, start_tmp, end_tmp, l->next_level);
    for (int i = l_index+1; i > 0; --i) {
      level_struct *coarse = levels[i];
      // finer level (towards finest)
      level_struct *fine   = levels[i - 1];

      compute_core_start_end(0, fine->inner_vector_size, &start_tmp, &end_tmp, fine, threading);
      // prolongate from coarse to fine
      apply_P_PRECISION(h->mlmc_b2, h->mlmc_b1, fine, threading);
      // copy result to fine level
      vector_PRECISION_copy(h->mlmc_b1, h->mlmc_b2, start_tmp, end_tmp, fine);

      if(g.my_rank == 0) {
        printf("TERM 2: Prolongating from depth %d to %d, \tfunction called at depth %d\n\n", coarse->depth, fine->depth, l->depth);
      }
    }

    // Copy result in second term vector
    compute_core_start_end(0, finest_l->inner_vector_size, &start_tmp, &end_tmp, finest_l, threading);
    //copy for dot product
    vector_PRECISION_copy(second_term, h->mlmc_b1, start_tmp, end_tmp, finest_l);
  }
  
  // subtract the results
  {
    int start, end;
    //gmres_PRECISION_struct* p = get_p_struct_PRECISION( finest_l);
    compute_core_start_end( 0, finest_l->inner_vector_size, &start, &end, finest_l, threading );
    vector_PRECISION_minus( out, first_term, second_term, start, end, finest_l);
  }
  
  PUBLIC_FREE( first_term, complex_PRECISION, finest_l->inner_vector_size );
  PUBLIC_FREE( second_term, complex_PRECISION, finest_l->inner_vector_size );  

}



//TODO: better naming (?)
// builds C_ij by calling mlmc_difference and mlmc_non_difference
complex_PRECISION connected_outer_PRECISION( int type_appl, level_struct *l, hutchinson_PRECISION_struct* h, struct Thread *threading ){

  complex_PRECISION trace = 0.0;
  
  level_struct *finest_l = h->finest_level;
  level_struct* lx_i;
  level_struct* lx_j;
  lx_i = h->lx_i;
  lx_j = h->lx_j;
  
  int i = lx_i->depth;
  int j = lx_j->depth;
  
  // Allocate buffer at finest
  vector_PRECISION op1 = NULL;
  PUBLIC_MALLOC( op1, complex_PRECISION, finest_l->inner_vector_size );
  
  // for all but coarsest level
  //lx = l;
  
  // CONNECTED STEP #1. apply \Pi_{t+t'}   (\Pi_{t+t'} x)
  {
    // Project rademacher vector into t+t'
    int bufft = g.time_slice;
    // TODO : check : is this assuming periodic in time ?
    g.time_slice = g.time_slice + g.time_slice_inner_connected;
    g.time_slice = g.time_slice%g.global_lattice[0][0];
    vector_PRECISION_ghg( h->rademacher_vector, 0, finest_l->inner_vector_size, finest_l );
    g.time_slice = bufft;
  }

  // CONNECTED STEP #2. apply Operator S_j    (S_j \Pi_{t+t'} x)
  {
    if( j == g.num_levels - 1){   
      connected_mlmc_PRECISION_non_difference( op1, h->rademacher_vector, lx_j,  h, threading );
    }
    else{
      connected_mlmc_PRECISION_difference( op1, h->rademacher_vector,  lx_j,  h, threading );
    }
  }
  
  // CONNECTED STEP #3. apply \Pi{t’}    (\Pi_{t'} S_j \Pi_{t+t'} x) 
  {
    int bufft = g.time_slice;
    // TODO : check : is this assuming periodic in time ?
    g.time_slice = g.time_slice_inner_connected;
    vector_PRECISION_ghg( op1, 0, finest_l->inner_vector_size, finest_l );
    g.time_slice = bufft;
  }
  
  // CONNECTED STEP #4. apply Difference operator S_i    (S_i \Pi_{t'} S_j \Pi_{t+t'} x)
  {
    if( i == g.num_levels - 1){   
      connected_mlmc_PRECISION_non_difference( op1, op1,  lx_i,  h, threading );
    }
    else{
      connected_mlmc_PRECISION_difference( op1, op1,  lx_i,  h, threading );
    }
  }

  
  // CONNECTED STEP #5 do the dot product AT FINEST       (x^H \Pi_{t+t'} S_i \Pi_{t'} S_j \Pi_{t+t'} x))
  {
    gmres_PRECISION_struct* p = get_p_struct_PRECISION( finest_l);

    //if ( type_appl==-1 ) {
    //  if(g.trace_deflation_type[l->depth] != 0){
    //    hutchinson_deflate_vector_PRECISION(h->mlmc_b1, l, threading);
    //  }
    trace = trace_inner_product_PRECISION( h->rademacher_vector, op1, 0, finest_l->inner_vector_size, finest_l, threading );
    //  } else {
    //  aux = global_inner_product_PRECISION( l->powerit_PRECISION.vecs[type_appl], h->mlmc_b1, p->v_start, p->v_end, l, threading );
    // }
  }
  
  PUBLIC_FREE( op1, complex_PRECISION, finest_l->inner_vector_size );
  
  return trace;

}

void connected_split_PRECISION_orthogonal( vector_PRECISION out, vector_PRECISION in, level_struct *l,  hutchinson_PRECISION_struct* h, struct Thread *threading ){
  // store from fine to coarsest lvl in an array
  level_struct *finest_l = h->finest_level;
  level_struct *levels[g.num_levels];
  int lvl_nr = 0, l_index = -1;  (void)l_index;

  if(g.my_rank == 0) {
    printf("\n\n------------------function at level %d ------------------\n\n", l->depth);fflush(0);
  }

  for (level_struct *l_tmp = finest_l; l_tmp != NULL; l_tmp = l_tmp->next_level) {
    levels[lvl_nr] = l_tmp;
    if(g.my_rank == 0) printf("Storing depth %d \t Function called at depth %d\n", l_tmp->depth, l->depth);
    if (l_tmp == l) l_index = lvl_nr;
    lvl_nr++;
  }
  
  // Preallocate first and second term vectors at the finest
  vector_PRECISION first_term = NULL;
  PUBLIC_MALLOC( first_term, complex_PRECISION, finest_l->inner_vector_size );
  vector_PRECISION second_term = NULL;
  PUBLIC_MALLOC( second_term, complex_PRECISION, finest_l->inner_vector_size );
    
  // Apply Gamma_5 at the finest
  {    
    gamma5_PRECISION( h->mlmc_testing, in, finest_l, threading );
  }
  
  // LEVEL DIFFERENCE : FIRST TERM : result stored in first_term
  {
    int start, end;
    gmres_PRECISION_struct* p = get_p_struct_PRECISION( l );

    compute_core_start_end( 0, finest_l->inner_vector_size, &start, &end, finest_l, threading );

    // copy buffer to be restricted below
    //if ( type_appl==-1 ) {
      vector_PRECISION_copy( h->mlmc_b1, h->mlmc_testing, start, end, finest_l );
    //} else {
    // // vector_PRECISION_copy( p->b, l->powerit_PRECISION.vecs[type_appl], start, end, l );
    // }

    // Restrict rhs from finest to current level
    int start_tmp, end_tmp, start_tmp_c, end_tmp_c;

    for (int i = 0; i < l_index ; ++i) {
      level_struct *fine   = levels[i];
      level_struct *coarse = levels[i + 1];

      compute_core_start_end(0, fine->inner_vector_size, &start_tmp, &end_tmp, fine, threading);
      compute_core_start_end(0, coarse->inner_vector_size, &start_tmp_c, &end_tmp_c, coarse, threading);
      // restrict from fine to coarse
      apply_R_PRECISION(h->mlmc_b2, h->mlmc_b1, fine, threading);
      vector_PRECISION_copy(h->mlmc_b1, h->mlmc_b2, start_tmp_c, end_tmp_c, coarse);

      if (g.my_rank == 0)
        printf("TERM 1: Restricting from depth %d to %d, \tfunction called at depth %d \n", fine->depth, coarse->depth, l->depth);
    }

    // Compute first term D_l^{-1}\hat{P}
    compute_core_start_end( 0, l->inner_vector_size, &start, &end, l, threading );
    vector_PRECISION_copy( p->b, h->mlmc_b1, start, end, l);
    apply_solver_PRECISION( l, threading );
    vector_PRECISION_copy( first_term, p->x, start, end, l);

    // Compute second term P_lP_l^H D_l^{-1}\hat{P}
    vector_PRECISION_copy( h->mlmc_b1, first_term, start, end, l);
    //Apply P_lP_l^H at current level
    apply_R_PRECISION(h->mlmc_b2, h->mlmc_b1, l, threading);
    apply_P_PRECISION(second_term, h->mlmc_b2, l, threading);
    
    // compute (D_l^{-1} - P_lP_l^H D_l^{-1})\hat{P} at current level
    vector_PRECISION_minus( h->mlmc_b1, first_term, second_term, start, end, l);

    
    // Prolongate vector to finest level
    for (int i = l_index; i > 0; --i) {
      level_struct *coarse = levels[i];
      level_struct *fine   = levels[i - 1];  /* finer level (towards finest) */

      compute_core_start_end(0, fine->inner_vector_size, &start_tmp, &end_tmp, fine, threading);
      apply_P_PRECISION(h->mlmc_b2, h->mlmc_b1, fine, threading);  /* prolongate from coarse to fine */
      vector_PRECISION_copy(h->mlmc_b1, h->mlmc_b2, start_tmp, end_tmp, fine);  /* copy result to fine level */

      if(g.my_rank == 0) {
        printf("TERM 1: Prolongating from depth %d to %d, \tfunction called at depth %d\n\n", coarse->depth, fine->depth, l->depth);
      }
    }
    
    int start_f, end_f;
    compute_core_start_end(0, finest_l->inner_vector_size, &start_f, &end_f, finest_l, threading);
    vector_PRECISION_copy(out, h->mlmc_b1, start_f, end_f, finest_l);
  }
   
  PUBLIC_FREE( first_term, complex_PRECISION, finest_l->inner_vector_size );
  PUBLIC_FREE( second_term, complex_PRECISION, finest_l->inner_vector_size );  

}



void connected_split_PRECISION_intermediate( vector_PRECISION out, vector_PRECISION in, level_struct *l,  hutchinson_PRECISION_struct* h, struct Thread *threading ){
  // store from fine to coarsest lvl in an array
  level_struct *finest_l = h->finest_level;
  level_struct *levels[g.num_levels];
  int lvl_nr = 0, l_index = -1;  (void)l_index;

  if(g.my_rank == 0) {
    printf("\n\n------------------function at level %d ------------------\n\n", l->depth);fflush(0);
  }

  for (level_struct *l_tmp = finest_l; l_tmp != NULL; l_tmp = l_tmp->next_level) {
    levels[lvl_nr] = l_tmp;
    if(g.my_rank == 0) printf("Storing depth %d \t Function called at depth %d\n", l_tmp->depth, l->depth);
    if (l_tmp == l) l_index = lvl_nr;
    lvl_nr++;
  }
  
  // Preallocate first and second term vectors at the finest
  vector_PRECISION first_term = NULL;
  PUBLIC_MALLOC( first_term, complex_PRECISION, finest_l->inner_vector_size );
  vector_PRECISION second_term = NULL;
  PUBLIC_MALLOC( second_term, complex_PRECISION, finest_l->inner_vector_size );
    
  // Apply Gamma_5 at the finest
  {    
    gamma5_PRECISION( h->mlmc_testing, in, finest_l, threading );
  }
  
  // LEVEL DIFFERENCE : FIRST TERM : result stored in first_term
  {
    int start, end;
    gmres_PRECISION_struct* p = get_p_struct_PRECISION( l );
    gmres_PRECISION_struct* p_next = get_p_struct_PRECISION( l->next_level );

    compute_core_start_end( 0, finest_l->inner_vector_size, &start, &end, finest_l, threading );

    // copy buffer to be restricted below
    //if ( type_appl==-1 ) {
      vector_PRECISION_copy( h->mlmc_b1, h->mlmc_testing, start, end, finest_l );
    //} else {
    // // vector_PRECISION_copy( p->b, l->powerit_PRECISION.vecs[type_appl], start, end, l );
    // }

    // Restrict rhs from finest to current level
    int start_tmp, end_tmp, start_tmp_c, end_tmp_c;

    for (int i = 0; i < l_index ; ++i) {
      level_struct *fine   = levels[i];
      level_struct *coarse = levels[i + 1];

      compute_core_start_end(0, fine->inner_vector_size, &start_tmp, &end_tmp, fine, threading);
      compute_core_start_end(0, coarse->inner_vector_size, &start_tmp_c, &end_tmp_c, coarse, threading);
      // restrict from fine to coarse
      apply_R_PRECISION(h->mlmc_b2, h->mlmc_b1, fine, threading);
      vector_PRECISION_copy(h->mlmc_b1, h->mlmc_b2, start_tmp_c, end_tmp_c, coarse);

      if (g.my_rank == 0)
        printf("TERM 1: Restricting from depth %d to %d, \tfunction called at depth %d \n", fine->depth, coarse->depth, l->depth);
    }

    // Solve at current level
    compute_core_start_end( 0, l->inner_vector_size, &start, &end, l, threading );
    vector_PRECISION_copy( p->b, h->mlmc_b1, start, end, l);
    apply_solver_PRECISION( l, threading );

    //Apply P_lP_l^H at current level
    apply_R_PRECISION(h->mlmc_b2, p->x, l, threading);
    apply_P_PRECISION(first_term, h->mlmc_b2, l, threading);
    
    // compute second term P_lD_{l+1}^{-1}P_l^H
    apply_R_PRECISION(p_next->b, h->mlmc_b1, l, threading);
    apply_solver_PRECISION( l->next_level, threading );
    apply_P_PRECISION(second_term, p_next->x, l, threading);
    
    // compute P_lP_l^H -  P_lD_{l+1}^{-1}P_l^H at current level
    vector_PRECISION_minus( h->mlmc_b1, first_term, second_term, start, end, l);

    
    // Prolongate vector to finest level
    for (int i = l_index; i > 0; --i) {
      level_struct *coarse = levels[i];
      level_struct *fine   = levels[i - 1];  /* finer level (towards finest) */

      compute_core_start_end(0, fine->inner_vector_size, &start_tmp, &end_tmp, fine, threading);
      apply_P_PRECISION(h->mlmc_b2, h->mlmc_b1, fine, threading);  /* prolongate from coarse to fine */
      vector_PRECISION_copy(h->mlmc_b1, h->mlmc_b2, start_tmp, end_tmp, fine);  /* copy result to fine level */

      if(g.my_rank == 0) {
        printf("TERM 1: Prolongating from depth %d to %d, \tfunction called at depth %d\n\n", coarse->depth, fine->depth, l->depth);
      }
    }

    int start_f, end_f;
    compute_core_start_end(0, finest_l->inner_vector_size,
                       &start_f, &end_f, finest_l, threading);

    vector_PRECISION_copy(out, h->mlmc_b1, start_f, end_f, finest_l);

  }
   
  PUBLIC_FREE( first_term, complex_PRECISION, finest_l->inner_vector_size );
  PUBLIC_FREE( second_term, complex_PRECISION, finest_l->inner_vector_size );  

}



//TODO: better naming (?)
// builds C_ij for split by calling split_intermediate, split_orthogonal and mlmc_non_difference
complex_PRECISION connected_outer_split_PRECISION( int type_appl, level_struct *l, hutchinson_PRECISION_struct* h, struct Thread *threading ){

  complex_PRECISION trace = 0.0;
  
  level_struct *finest_l = h->finest_level;
  level_struct* lx_i;
  level_struct* lx_j;
  lx_i = h->lx_i;
  lx_j = h->lx_j;

  int l_op = h->l_op; //0 is orthogonal and 1 is full rank
  int r_op = h->r_op;
  int i = lx_i->depth; //When i or j = g.num_levels-1 we deal with the coarsest level operator
  int j = lx_j->depth;
  
  // Allocate buffer at finest
  vector_PRECISION op1 = NULL;
  PUBLIC_MALLOC( op1, complex_PRECISION, finest_l->inner_vector_size );
  
  // for all but coarsest level
  //lx = l;
  
  // CONNECTED SPLIT STEP #1. apply \Pi_{t+t'}   (\Pi_{t+t'} x)
  {
    // Project rademacher vector into t+t'
    int bufft = g.time_slice;
    // TODO : check : is this assuming periodic in time ?
    g.time_slice = g.time_slice + g.time_slice_inner_connected;
    g.time_slice = g.time_slice%g.global_lattice[0][0];
    vector_PRECISION_ghg( h->rademacher_vector, 0, finest_l->inner_vector_size, finest_l );
    g.time_slice = bufft;
  }

  // CONNECTED SPLIT STEP #2. apply Operator S_j    (S_j \Pi_{t+t'} x)
  if(g.my_rank == 0) printf("\napplying right operator, r_op = %d, j = %d", r_op, j);
  {
    if(r_op == 0 && j != g.num_levels-1){
      if(g.my_rank == 0) printf("\napplying Right ORTHOGONAL\n");
      connected_split_PRECISION_orthogonal(op1, h->rademacher_vector, lx_j,  h, threading );
    }
    else if(r_op == 1 && j != g.num_levels-1){
      if(g.my_rank == 0) printf("\napplying Right FULL RANK\n");
      connected_split_PRECISION_intermediate(op1, h->rademacher_vector, lx_j,  h, threading );
    }
    else if(j == g.num_levels-1){
      if(g.my_rank == 0) printf("\napplying Right COARSEST\n");
      connected_mlmc_PRECISION_non_difference( op1, h->rademacher_vector, lx_j,  h, threading ); 
    }
  }
  if(g.my_rank == 0) printf("\nright operator applied\n");
  
  // CONNECTED SPLIT STEP #3. apply \Pi{t’}    (\Pi_{t'} S_j \Pi_{t+t'} x) 
  {
    int bufft = g.time_slice;
    // TODO : check : is this assuming periodic in time ?
    g.time_slice = g.time_slice_inner_connected;
    vector_PRECISION_ghg( op1, 0, finest_l->inner_vector_size, finest_l );
    g.time_slice = bufft;
  }
  
  // CONNECTED SPLIT STEP #4. apply Difference operator S_i    (S_i \Pi_{t'} S_j \Pi_{t+t'} x)
  if(g.my_rank == 0) printf("\napplying left operator, l_op = %d, i = %d", l_op, i);
  {
    if(l_op == 0 && i != g.num_levels-1){
      if(g.my_rank == 0) printf("\napplying Left ORTHOGONAL\n");
      connected_split_PRECISION_orthogonal(op1, op1, lx_i,  h, threading );
    }
    else if(l_op == 1 && i != g.num_levels-1){
        if(g.my_rank == 0) printf("\napplying Left FULL RANK\n");
        connected_split_PRECISION_intermediate(op1, op1, lx_i,  h, threading );
    }
    else if(i == g.num_levels-1){
      if(g.my_rank == 0) printf("\napplying Left COARSEST\n");
      connected_mlmc_PRECISION_non_difference( op1, op1, lx_i,  h, threading );
    }
  }
  if(g.my_rank == 0) printf("\nleft operator applied\n");

  
  // CONNECTED SPLIT STEP #5 do the dot product AT FINEST       (x^H \Pi_{t+t'} S_i \Pi_{t'} S_j \Pi_{t+t'} x))
  {
    gmres_PRECISION_struct* p = get_p_struct_PRECISION( finest_l);

    //if ( type_appl==-1 ) {
    //  if(g.trace_deflation_type[l->depth] != 0){
    //    hutchinson_deflate_vector_PRECISION(h->mlmc_b1, l, threading);
    //  }
    trace = trace_inner_product_PRECISION( h->rademacher_vector, op1, 0, finest_l->inner_vector_size, finest_l, threading );
    //  } else {
    //  aux = global_inner_product_PRECISION( l->powerit_PRECISION.vecs[type_appl], h->mlmc_b1, p->v_start, p->v_end, l, threading );
    // }
  }
  
  PUBLIC_FREE( op1, complex_PRECISION, finest_l->inner_vector_size );
  
  return trace;

}

//This function returns 0 if the corresponding CL terms has not been computed yet and 1 otherwise
int CL_check_PRECISION(int l_op, int r_op, int i, int j){
  int CL_status = 0;

  if(l_op == 0 && r_op == 0)
    CL_status = 0;

  if(l_op == 0 && r_op == 1){
    if( i == g.num_levels-1 && j != g.num_levels-1)
      CL_status = 0;
    else
      CL_status = 1;
  }

  if(l_op == 1 && r_op == 0){
    if(i != g.num_levels-1)
      CL_status = 0;
    else
      CL_status = 1;
  }

  if(l_op == 1 && r_op == 1)
    CL_status = 1;

  return CL_status;
}



complex_PRECISION g5_3D_connected_split_driver_PRECISION( level_struct *l, struct Thread *threading ){
  int i, j;
  complex_PRECISION trace = 0.0;
  complex_PRECISION global_connected_trace = 0.0;
  struct sample estimate;
  hutchinson_PRECISION_struct* h = &(l->h_PRECISION);

  h->l_op = 0; //0:ORT 1:FR
  h->r_op = 0; //we have CL when i or j == g.num_levels-1
  
  //To avoid creating the same operator more than once we need the CL_check function for the CL terms
  //In a n lvl method the number of operators for the connected split is (2n-1)*(2n-1)
  //E.g. with the 4 nested for loops in a n=3 lvl method we would have 36 iterations but only 25 operators
  //The extra iterations will create repeated CL terms because when i or j == g.num_levels-1 the index of the corresponding G term
  //is always CL and not ORT or FR
  //In a n lvl method the total number of ORT-CL FR-CL CL-ORT CL-FR CL-CL terms is 4n-3 

  for(h->l_op = 0; h->l_op < 2; h->l_op++){
    for(h->r_op = 0; h->r_op < 2; h->r_op++){

      h->lx_i = h->finest_level;
      h->lx_j = h->finest_level;

      for(i = 0; i < g.num_levels; i++){
        h->lx_j = h->finest_level;
	for(j = 0; j < g.num_levels; j++){


          const char *lop_string, *rop_string;
          if(h->l_op == 0 && i != g.num_levels-1)
            lop_string = "ORT";
          else if(h->l_op == 1 && i != g.num_levels-1)
            lop_string = "FR";
          else if(i == g.num_levels-1)
            lop_string = "CL";

          if(h->r_op == 0 && j != g.num_levels-1)
              rop_string = "ORT";
            else if(h->r_op == 1 && j != g.num_levels-1)
              rop_string = "FR";
            else if(j == g.num_levels-1)
              rop_string = "CL";

          if (i == g.num_levels-1 || j == g.num_levels-1) {
            if (CL_check_PRECISION(h->l_op, h->r_op, i, j)) {
              if (j < g.num_levels - 1) {
                h->lx_j = h->lx_j->next_level;
              }
              continue;
            }
          }

          if(g.my_rank == 0) printf("\n Computing trace for G_{%d, %d}^{%s-%s}(t=%d) \n", h->lx_i->depth, h->lx_j->depth, lop_string, rop_string, g.time_slice);

          trace = 0.0;
          h->hutch_compute_one_sample = connected_outer_split_PRECISION;

          // Call to blind from t+t’ with t’ = 0, ..., T-1
          for ( g.time_slice_inner_connected=0; g.time_slice_inner_connected < g.global_lattice[0][0]; g.time_slice_inner_connected++ ) {

            // colors and dilution are looped inside the blind routine, with one shared z0 per sample
    estimate = probing_hutchinson_blind_PRECISION( h->finest_level, h, 0, i, threading );
    trace += estimate.acc_trace / estimate.sample_size;
    { char label[96]; snprintf( label, 96, "G_{%d,%d}^{%s-%s} (average over t')", h->lx_i->depth, h->lx_j->depth, lop_string, rop_string ); report_piece_PRECISION( label, estimate, 1.0/g.global_lattice[0][0], threading ); }
           }

        trace = trace/g.time_slice_inner_connected;
        global_connected_trace += trace;
        //int nlevs = g.num_levels;
        //int idx = h->lx_i->depth * nlevs + h->lx_j->depth;
        if (g.my_rank == 0) printf("\n Trace of G_{%d, %d}^{%s, %s}(t=%d) = %f +%f\n", h->lx_i->depth, h->lx_j->depth, lop_string, rop_string, g.time_slice, CSPLIT(trace));



        if (j < g.num_levels - 1)
          h->lx_j = h->lx_j->next_level;
      }
        if (i < g.num_levels - 1)
          h->lx_i = h->lx_i->next_level;
    }
  }
}

  return global_connected_trace;
  
}




/// if type_appl==-1 : Hutchinson-like
// else : direct term, where type_appl is the index of the deflation vector to apply
//        the operator on
complex_PRECISION hutchinson_plain_PRECISION( int type_appl, level_struct *l, hutchinson_PRECISION_struct* h, struct Thread *threading ){
  complex_PRECISION aux = 0.0;
  {
    int start, end;
    gmres_PRECISION_struct* p = get_p_struct_PRECISION( l );
    compute_core_start_end( 0, l->inner_vector_size, &start, &end, l, threading );

    if ( type_appl==-1 ) {
      vector_PRECISION_copy( p->b, h->rademacher_vector, start, end, l );
   // } else {
   //   vector_PRECISION_copy( p->b, l->powerit_PRECISION.vecs[type_appl], start, end, l );
   //}
    }
  }

  {
    apply_solver_PRECISION( l, threading );
  }

  // subtract the results and perform dot product
  {
    int start, end;
    gmres_PRECISION_struct* p = get_p_struct_PRECISION( l );
    compute_core_start_end( 0, l->inner_vector_size, &start, &end, l, threading );

    if ( type_appl==-1 ) {
   //   if(g.trace_deflation_type[l->depth] != 0){
   //     hutchinson_deflate_vector_PRECISION(p->x, l, threading); 
   //   }
      aux = trace_inner_product_PRECISION( h->rademacher_vector, p->x, 0, l->inner_vector_size, l, threading );
   //} else {
   //   vector_PRECISION_copy(l->powerit_PRECISION.vecs_buff1, p->x, start, end, l);  
   //   aux = global_inner_product_PRECISION( l->powerit_PRECISION.vecs[type_appl], l->powerit_PRECISION.vecs_buff1, p->v_start, p->v_end, l, threading );
    }

    return aux;  
  }
}

complex_PRECISION g5_hutchinson_plain_PRECISION( int type_appl, level_struct *l, hutchinson_PRECISION_struct* h, struct Thread *threading ){
  complex_PRECISION aux = 0.0;
  {
    int start, end;
    gmres_PRECISION_struct* p = get_p_struct_PRECISION( l );
    compute_core_start_end( 0, l->inner_vector_size, &start, &end, l, threading );

    // Apply Gamma out-of-place: b = gamma * z.
    // NEVER apply spin-off-diagonal gammas in place (out==in overwrites
    // spin components before their partners are read).
    if ( type_appl==-1 ) {
      if(g.gamma_idx == 5)      gamma5_PRECISION( p->b, h->rademacher_vector, l, threading );
      else if(g.gamma_idx == 0) gamma0_PRECISION( p->b, h->rademacher_vector, l, threading );
      else if(g.gamma_idx == 1) gamma1_PRECISION( p->b, h->rademacher_vector, l, threading );
      else if(g.gamma_idx == 2) gamma2_PRECISION( p->b, h->rademacher_vector, l, threading );
      else if(g.gamma_idx == 3) gamma3_PRECISION( p->b, h->rademacher_vector, l, threading );
      else vector_PRECISION_copy( p->b, h->rademacher_vector, start, end, l );  // identity
    }
  }

  {
    apply_solver_PRECISION( l, threading );
  }
  
    
  // subtract the results and perform dot product
  {
    int start, end;
    gmres_PRECISION_struct* p = get_p_struct_PRECISION( l );
    compute_core_start_end( 0, l->inner_vector_size, &start, &end, l, threading );

    if ( type_appl==-1 ) {
      aux = trace_inner_product_PRECISION( h->rademacher_vector, p->x, 0, l->inner_vector_size, l, threading );
    }

    return aux;  
  }
}

complex_PRECISION g5_hutchinson_driver_PRECISION( level_struct *l, struct Thread *threading ){
  complex_PRECISION trace = 0.0;
  struct sample estimate;
  hutchinson_PRECISION_struct* h = &(l->h_PRECISION);
  level_struct* lx = l;

  // set the pointer to the finest-level Hutchinson estimator
  h->hutch_compute_one_sample = g5_hutchinson_plain_PRECISION;

  // colors and dilution are looped inside the blind routine, with one shared z0 per sample
    estimate = probing_hutchinson_blind_PRECISION( lx, h, 0, 0, threading );
    trace += estimate.acc_trace / estimate.sample_size;
    report_piece_PRECISION( "plain Hutchinson (level 1)", estimate, 1.0, threading );

  return trace;
}



// the term tr( A_{l}^{-1} - P A_{l+1}^{-1} R )
complex_PRECISION hutchinson_mlmc_difference_PRECISION( int type_appl, level_struct *l, hutchinson_PRECISION_struct* h, struct Thread *threading ){
  // FIRST TERM : result stored in p->x
  // apply A_{l}^{-1}
  complex_PRECISION aux = 0.0;
  {
    int start, end;
    gmres_PRECISION_struct* p = get_p_struct_PRECISION( l );
    compute_core_start_end( 0, l->inner_vector_size, &start, &end, l, threading );

    if ( type_appl==-1 ) {
      vector_PRECISION_copy( p->b, h->rademacher_vector, start, end, l );
    //} else {
    //  vector_PRECISION_copy( p->b, l->powerit_PRECISION.vecs[type_appl], start, end, l );
    }

    // solution of this solve is in l->p_PRECISION.x
    apply_solver_PRECISION( l, threading );
  }

  // SECOND TERM : result stored in h->mlmc_b2
  // 1. Restrict
  // 2. invert
  // 3. Prolongate
  {

    if ( type_appl==-1 ) {
      apply_R_PRECISION( l->next_level->p_PRECISION.b, h->rademacher_vector, l, threading );
    //} else {
    //  apply_R_PRECISION( l->next_level->p_PRECISION.b, l->powerit_PRECISION.vecs[type_appl], l, threading );
    }

    // the input of this solve is l->next_level->p_PRECISION.x, the output l->next_level->p_PRECISION.b
    apply_solver_PRECISION( l->next_level, threading );
    apply_P_PRECISION( h->mlmc_b2, l->next_level->p_PRECISION.x, l, threading );
  }

  // subtract the results and perform dot product
  {
    int start, end;
    gmres_PRECISION_struct* p = get_p_struct_PRECISION( l);
    compute_core_start_end( 0, l->inner_vector_size, &start, &end, l, threading );
    vector_PRECISION_minus( h->mlmc_b1, p->x, h->mlmc_b2, start, end, l);

    if ( type_appl==-1 ) {
      //if(g.trace_deflation_type[l->depth] != 0){
      //  hutchinson_deflate_vector_PRECISION(h->mlmc_b1, l, threading); 
      //}
      aux = trace_inner_product_PRECISION( h->rademacher_vector, h->mlmc_b1, 0, l->inner_vector_size, l, threading );
    //} else {
    //  aux = global_inner_product_PRECISION( l->powerit_PRECISION.vecs[type_appl], h->mlmc_b1, p->v_start, p->v_end, l, threading );
    }

    return aux;
  }
}


complex_PRECISION mlmc_hutchinson_driver_PRECISION( level_struct *l, struct Thread *threading ){
  int i;
  complex_PRECISION trace = 0.0;
  struct sample estimate;
  hutchinson_PRECISION_struct* h = &(l->h_PRECISION);
  level_struct* lx;

  // for all but coarsest level
  lx = l;
  for( i=0; i<g.num_levels-1; i++ ){ 
    // set the pointer to the mlmc difference operator
    h->hutch_compute_one_sample = hutchinson_mlmc_difference_PRECISION;
    // colors and dilution are looped inside the blind routine, with one shared z0 per sample
    estimate = probing_hutchinson_blind_PRECISION( lx, h, 0, i, threading );
    trace += estimate.acc_trace / estimate.sample_size;
    { char label[64]; snprintf( label, 64, "level %d - level %d", i+1, i+2 ); report_piece_PRECISION( label, estimate, 1.0, threading ); }
    // if deflation vectors are available
    //if(g.trace_deflation_type[lx->depth] != 0){
    //  trace += hutchinson_deflated_direct_term_PRECISION( lx, h, threading );
    //}
    lx = lx->next_level;
  }

  // coarsest level
  // set the pointer to the coarsest-level Hutchinson estimator
  h->hutch_compute_one_sample = hutchinson_plain_PRECISION;
  // colors and dilution are looped inside the blind routine, with one shared z0 per sample
    estimate = probing_hutchinson_blind_PRECISION( lx, h, 0, i, threading );
    trace += estimate.acc_trace / estimate.sample_size;
    { char label[64]; snprintf( label, 64, "coarsest (level %d)", g.num_levels ); report_piece_PRECISION( label, estimate, 1.0, threading ); }
  // if deflation vectors are available
  //if(g.trace_deflation_type[lx->depth] != 0){
  //  trace += hutchinson_deflated_direct_term_PRECISION( lx, h, threading );
  //}

  return trace;
}

// the term tr( A_{l}^{-1}(I - P_{l} P_{l}^{H})  )
complex_PRECISION hutchinson_split_orthogonal_PRECISION( int type_appl, level_struct *l, hutchinson_PRECISION_struct* h, struct Thread *threading ){
  // 1. project
  // 2. invert

  // FIRST TERM
  {
    int start, end;
    gmres_PRECISION_struct* p = get_p_struct_PRECISION( l );

    if ( type_appl==-1 ) {
      apply_R_PRECISION( h->mlmc_b2, h->rademacher_vector, l, threading );
    //} else {
    //  apply_R_PRECISION( h->mlmc_b2, l->powerit_PRECISION.vecs[type_appl], l, threading );
    }

    apply_P_PRECISION( h->mlmc_b1, h->mlmc_b2, l, threading );
    compute_core_start_end( 0, l->inner_vector_size, &start, &end, l, threading );

    if ( type_appl==-1 ) {
      vector_PRECISION_minus( p->b, h->rademacher_vector, h->mlmc_b1, start, end, l );
    //} else {
    //  vector_PRECISION_minus( p->b, l->powerit_PRECISION.vecs[type_appl], h->mlmc_b1, start, end, l );
    }
  }
    
  // SECOND "factor"
  {
    apply_solver_PRECISION( l, threading );
  }

  // perform dot product
  {
    int start, end;
    complex_PRECISION aux = 0.0;
    gmres_PRECISION_struct* p = get_p_struct_PRECISION( l );
    compute_core_start_end( 0, l->inner_vector_size, &start, &end, l, threading );
      
    
    apply_R_PRECISION( h->mlmc_b2, p->x, l, threading );
    apply_P_PRECISION( h->mlmc_b1, h->mlmc_b2, l, threading );

    vector_PRECISION_minus( h->mlmc_b1, p->x, h->mlmc_b1, start, end, l );
    
    if ( type_appl==-1 ) {
    //  if(g.trace_deflation_type[l->depth] != 0){
    //    hutchinson_deflate_vector_PRECISION(h->mlmc_b1, l, threading); 
    //  }
      aux = trace_inner_product_PRECISION( h->rademacher_vector, h->mlmc_b1, 0, l->inner_vector_size, l, threading );   
    //} else {  
    //  aux = global_inner_product_PRECISION( l->powerit_PRECISION.vecs[type_appl], h->mlmc_b1, p->v_start, p->v_end, l, threading );
    }
    return aux; 
  }
}

// the term tr( R A_{l}^{-1} P - A_{l+1}^{-1} )
complex_PRECISION hutchinson_split_intermediate_PRECISION( int type_appl, level_struct *l, hutchinson_PRECISION_struct* h, struct Thread *threading ){

  complex_PRECISION aux = 0.0;
  // FIRST TERM : result stored in h->mlmc_b1
  // 1. prolongate
  // 2. invert
  // 3. restrict
  {
    gmres_PRECISION_struct* p = get_p_struct_PRECISION( l );

    if ( type_appl==-1 ) {
      apply_P_PRECISION( p->b, h->rademacher_vector, l, threading );
   // } else {
   //   apply_P_PRECISION( p->b, l->powerit_PRECISION.vecs[type_appl], l, threading );
    }
    // the input of this solve is p->x, the output p->b
    apply_solver_PRECISION( l, threading );
    apply_R_PRECISION( h->mlmc_b1, p->x, l, threading );
  }

  // SECOND TERM : result stored in h->mlmc_b2

  // apply A_{l+1}^{-1}
  {
    int start, end;
    gmres_PRECISION_struct* p = get_p_struct_PRECISION( l->next_level );
    compute_core_start_end( 0, l->next_level->inner_vector_size, &start, &end, l->next_level, threading );

    if ( type_appl==-1 ) {
      vector_PRECISION_copy( p->b, h->rademacher_vector, start, end, l->next_level );
    //} else {
    //  vector_PRECISION_copy( p->b, l->powerit_PRECISION.vecs[type_appl], start, end, l->next_level );
    }

    // solution of this solve is in l->next_level->p_PRECISION.x
    apply_solver_PRECISION( l->next_level, threading );
    vector_PRECISION_copy( h->mlmc_b2, l->next_level->p_PRECISION.x, start, end, l->next_level );
  }

  // subtract the results and perform dot product
  {
    int start, end;
    gmres_PRECISION_struct* p = get_p_struct_PRECISION( l->next_level );
    compute_core_start_end( 0, l->next_level->inner_vector_size, &start, &end, l->next_level, threading );
    vector_PRECISION_minus( h->mlmc_b1, h->mlmc_b1, h->mlmc_b2, start, end, l->next_level ); 

    //if(g.trace_deflation_type[l->depth] != 0){
    //  hutchinson_deflate_vector_PRECISION(h->mlmc_b1, l, threading); 
    //}

    if ( type_appl==-1 ) {
      aux = trace_inner_product_PRECISION( h->rademacher_vector, h->mlmc_b1, 0, l->next_level->inner_vector_size, l->next_level, threading );         
    //} else {
    //  aux = global_inner_product_PRECISION( l->powerit_PRECISION.vecs[type_appl], h->mlmc_b1, p->v_start, p->v_end, l->next_level, threading );         
    }
      
    return aux;
  }
}

complex_PRECISION split_mlmc_hutchinson_driver_PRECISION( level_struct *l, struct Thread *threading ){
  int i;
  complex_PRECISION trace = 0.0;
  struct sample estimate;
  hutchinson_PRECISION_struct* h = &(l->h_PRECISION);
  level_struct* lx=0;

  // for all but coarsest level
  lx = l;
  for( i=0; i<g.num_levels-1 ;i++ ){  
    // set the pointer to the split full rank operator
    h->hutch_compute_one_sample = hutchinson_split_intermediate_PRECISION;
    // colors and dilution are looped inside the blind routine, with one shared z0 per sample
    estimate = probing_hutchinson_blind_PRECISION( lx, h, 1, i+1, threading );
    trace += estimate.acc_trace / estimate.sample_size;
    { char label[64]; snprintf( label, 64, "split full rank (level %d)", i+1 ); report_piece_PRECISION( label, estimate, 1.0, threading ); }
    START_MASTER(threading)
    if(g.my_rank == 0){
      g.variances[i]=0.0;
    }
    END_MASTER(threading)

    // if deflation vectors are available
    //if( g.trace_deflation_type[lx->depth] != 0 ){
    //  if( g.trace_deflation_type[lx->depth]==3 || g.trace_deflation_type[lx->depth]==5 ){
    //    trace += hutchinson_deflated_direct_term_PRECISION(lx, h, threading);
    //  }
    //}
    lx = lx->next_level;    
  }

  // for all but coarsest level
  lx = l;
  for( i=0; i<g.num_levels-1;i++ ){      
    // set the pointer to the split orthogonal operator
    h->hutch_compute_one_sample = hutchinson_split_orthogonal_PRECISION;
    // colors and dilution are looped inside the blind routine, with one shared z0 per sample
    estimate = probing_hutchinson_blind_PRECISION( lx, h, 0, i, threading );
    trace += estimate.acc_trace / estimate.sample_size;
    { char label[64]; snprintf( label, 64, "split orthogonal (level %d)", i+1 ); report_piece_PRECISION( label, estimate, 1.0, threading ); }

    // if deflation vectors are available
    //if( g.trace_deflation_type[lx->depth] != 0 ){
    //  if( g.trace_deflation_type[lx->depth]==4 || g.trace_deflation_type[lx->depth]==5 ){
    //    trace += hutchinson_deflated_direct_term_PRECISION(lx, h, threading);
    //  }
    //}
    lx = lx->next_level; 
  }

  // coarsest level
  // set the pointer to the coarsest-level Hutchinson estimator
  h->hutch_compute_one_sample = hutchinson_plain_PRECISION;
  // colors and dilution are looped inside the blind routine, with one shared z0 per sample
    estimate = probing_hutchinson_blind_PRECISION( lx, h, 0, i, threading );
    trace += estimate.acc_trace / estimate.sample_size;
    { char label[64]; snprintf( label, 64, "coarsest (level %d)", g.num_levels ); report_piece_PRECISION( label, estimate, 1.0, threading ); }

  return trace;
}

