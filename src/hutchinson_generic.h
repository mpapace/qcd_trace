#ifndef HUTCHINSON_PRECISION_HEADER
  #define HUTCHINSON_PRECISION_HEADER

  struct Thread;

  complex_PRECISION trace_inner_product_PRECISION( vector_PRECISION a, vector_PRECISION b, int start, int end, level_struct *l, struct Thread *threading );
  void apply_P_PRECISION( vector_PRECISION out, vector_PRECISION in, level_struct* l, struct Thread *threading );
  void apply_R_PRECISION( vector_PRECISION out, vector_PRECISION in, level_struct* l, struct Thread *threading );
  void apply_gamma_idx_PRECISION( vector_PRECISION out, vector_PRECISION in, level_struct *l, struct Thread *threading );
  void apply_M_timeslice_PRECISION( vector_PRECISION out, vector_PRECISION in, level_struct *l, hutchinson_PRECISION_struct* h, struct Thread *threading );
  int apply_solver_PRECISION( level_struct* l, struct Thread *threading );
   
  gmres_PRECISION_struct* get_p_struct_PRECISION( level_struct* l );

  complex_PRECISION hutchinson_plain_PRECISION( int type_appl, level_struct *l, hutchinson_PRECISION_struct* h, struct Thread *threading );

  struct sample probing_hutchinson_blind_PRECISION( level_struct *l, hutchinson_PRECISION_struct* h, int type,
                                                    int dil_level, struct Thread *threading );
  level_struct* trace_vector_level_PRECISION( level_struct *l, hutchinson_PRECISION_struct* h, int type );
  void probe_create_PRECISION( level_struct *vl, hutchinson_PRECISION_struct* h, int dil_level, struct Thread *threading );
  void rademacher_create_PRECISION( level_struct *vl, hutchinson_PRECISION_struct* h, struct Thread *threading );

  void hutchinson_diver_PRECISION_init( level_struct *l, struct Thread *threading );
  void hutchinson_diver_PRECISION_alloc( level_struct *l, struct Thread *threading );
  void hutchinson_diver_PRECISION_free( level_struct *l, struct Thread *threading );
  complex_PRECISION hutchinson_driver_PRECISION( level_struct *l, struct Thread *threading );
  complex_PRECISION g5_hutchinson_driver_PRECISION( level_struct *l, struct Thread *threading );
  complex_PRECISION g5_3D_connected_hutchinson_driver_PRECISION( level_struct *l, struct Thread *threading );
  complex_PRECISION g5_3D_connected_hutchinson_plain_PRECISION( int type_appl, level_struct *l, hutchinson_PRECISION_struct* h, struct Thread *threading );
  complex_PRECISION g5_3D_connected_mlmc_driver_PRECISION( level_struct *l, struct Thread *threading );
  complex_PRECISION g5_3D_connected_split_driver_PRECISION( level_struct *l, struct Thread *threading );


  complex_PRECISION gamma_3D_hutchinson_driver_PRECISION( level_struct *l, struct Thread *threading );
  complex_PRECISION gamma_3D_hutchinson_plain_PRECISION( int type_appl, level_struct *l, hutchinson_PRECISION_struct* h, struct Thread *threading );
  complex_PRECISION gamma_3D_mlmc_difference_PRECISION( int type_appl, level_struct *l, hutchinson_PRECISION_struct* h, struct Thread *threading );
  complex_PRECISION gamma_3D_mlmc_hutchinson_driver_PRECISION( level_struct *l, struct Thread *threading );
  complex_PRECISION gamma_3D_mlmc_coarsest_PRECISION( int type_appl, level_struct *l, hutchinson_PRECISION_struct* h, struct Thread *threading );
  
  complex_PRECISION connected_outer_PRECISION( int type_appl, level_struct *l, hutchinson_PRECISION_struct* h, struct Thread *threading );

  complex_PRECISION mlmc_hutchinson_driver_PRECISION( level_struct *l, struct Thread *threading );
  complex_PRECISION split_mlmc_hutchinson_driver_PRECISION( level_struct *l, struct Thread *threading );


  
  //-------------------- Jose’s code below -----------------------------------------
  complex_PRECISION mlmc_hutchinson_g5_driver_PRECISION( level_struct *l, struct Thread *threading );
  complex_PRECISION hutchinson_mlmc_difference_PRECISION( int type_appl, level_struct *l, hutchinson_PRECISION_struct* h, struct Thread *threading );
  struct sample hutchinson_blind_g5_PRECISION( level_struct *l, int depth, hutchinson_PRECISION_struct* h, int type, struct Thread *threading );


#endif
