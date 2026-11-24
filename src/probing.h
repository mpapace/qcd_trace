void allocate_variances();
void mlmc_connected_print_variances();
void set_probing_variances_to_zero();
void setup_local_colors( level_struct *l );
void level_site_order( level_struct *lv, int *lex_to_vec );
void setup_site_time( level_struct *l );
int level_time_slice( int depth );
int local_colors_size( int level );
void dilution_check(int level);
#ifndef PROBING_HEADER
  #define PROBING_HEADER
    void graph_coloring( level_struct *l );
    //void set_probing_variances_to_zero();
#endif
