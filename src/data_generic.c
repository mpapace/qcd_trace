/*
 * Copyright (C) 2016, Matthias Rottmann, Artur Strebel, Gustavo Ramirez, Simon Heybrock, Simone Bacchio, Bjoern Leder, Issaku Kanamori, Tilmann Matthaei, Ke-Long Zhang.
 * 
 * This file is part of the DDalphaAMG solver library.
 * 
 * The DDalphaAMG solver library is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * The DDalphaAMG solver library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * 
 * You should have received a copy of the GNU General Public License
 * along with the DDalphaAMG solver library. If not, see http://www.gnu.org/licenses/.
 * 
 */

#include "main.h"
#include<unistd.h>


// vector storage for PRECISION precision
void vector_PRECISION_define( vector_PRECISION phi, complex_PRECISION value, int start, int end, level_struct *l ) {
  
  int thread = omp_get_thread_num();
  if(thread == 0 && start != end)
  PROF_PRECISION_START( _SET );
  if ( phi != NULL ) {
    int i;
    for ( i=start; i<end; i++ )
      phi[i] = value;
  } else {
    error0("Error in \"vector_PRECISION_define\": pointer is null\n");
  }
  if(thread == 0 && start != end)
  PROF_PRECISION_STOP( _SET, 1 );
}


void vector_PRECISION_define_random( vector_PRECISION phi, int start, int end, level_struct *l ) {
  
  int thread = omp_get_thread_num();
  if(thread == 0 && start != end)
  PROF_PRECISION_START( _SET );
  if ( phi != NULL ) {
    int i;
    for ( i=start; i<end; i++ )
      phi[i] = (PRECISION)(((double)rand()/(double)RAND_MAX))-0.5 + ( (PRECISION)((double)rand()/(double)RAND_MAX)-0.5)*_Complex_I;
  } else {
    error0("Error in \"vector_PRECISION_define_random\": pointer is null\n");
  }
  if(thread == 0 && start != end)
  PROF_PRECISION_STOP( _SET, 1 );
}

int compute_dilution_idx_PRECISION(level_struct *l, int global_idx, int dof){
  int dilution_idx = 1;
  int dilution_group = dof/g.dilution[l->depth];

  if(g.dilution[l->depth] != 3 && l->depth == 0)
    dilution_idx = (global_idx/dilution_group)%g.dilution[l->depth] + 1;

  if(g.dilution[l->depth] == 3 && l->depth == 0)
    dilution_idx = global_idx%3 + 1;

  return dilution_idx;

}

int get_dilution_value_PRECISION(int idx, int level){

    int val = 0;

    if(g.dilution[level] == 1) val = 1;

    if(g.dilution[level] == 2){
      if(g.dilution_count == 1 && idx <= 5) val = 1;
      if(g.dilution_count == 2 && idx > 5) val = 1;
    }

    if(g.dilution[level] == 3){
      if(g.dilution_count == 1 && (idx == 0 || idx == 3 || idx == 6 || idx == 9)) val = 1;
      if(g.dilution_count == 2 && (idx == 1 || idx == 4 || idx == 7 || idx == 10)) val = 1;
      if(g.dilution_count == 3 && (idx == 2 || idx == 5 || idx == 8 || idx == 11)) val = 1;
    }

    if(g.dilution[level] == 4){
      if(g.dilution_count == 1 && idx <= 2) val = 1;
      if(g.dilution_count == 2 && idx > 2 && idx <= 5) val = 1;
      if(g.dilution_count == 3 && idx > 5 && idx <= 8) val = 1;
      if(g.dilution_count == 4 && idx > 8) val = 1;
    }

    if(g.dilution[level] == 12){
      if(g.dilution_count == idx+1) val=1;
    }

    return val;
}

// Coarse-grid dilution: a coarse site carries nvar = 2*n_tv dof, chirality outer
// and test vector inner (index = chi*n_tv + v, see restrict_PRECISION).
// D = 1: none; D = 2*m with m | n_tv: chirality x m contiguous blocks of test
// vectors (D = 2 is pure chirality, D = nvar is full dof dilution).
// Piece g.dilution_count (1-based) is kept.
static int coarse_dilution_value_PRECISION( int k, int nvar, int D ) {
  if( D == 1 ) return 1;
  int n_tv = nvar/2, m = D/2, block = n_tv/m;
  int piece = (k/n_tv)*m + (k%n_tv)/block;
  return ( piece == g.dilution_count-1 ) ? 1 : 0;
}

// Fill e_tab[0..nvar-1] with the dilution pattern of the current g.dilution_count
// for a vector on level l, using the dilution scheme of level dil_level:
// spin/color schemes on the finest grid, chirality/test-vector schemes on coarse grids.
static void dilution_table_PRECISION( int *e_tab, level_struct *l, int dil_level ) {
  int nvar = l->num_lattice_site_var, D = g.dilution[dil_level];
  if( l->depth == 0 ){
    if( D != 1 && D != 2 && D != 3 && D != 4 && D != 12 )
      error0("dilution %d (level %d) is not a finest-grid dilution (1, 2, 3, 4, 12)\n", D, dil_level);
    for( int k=0; k<nvar; k++ )
      e_tab[k] = get_dilution_value_PRECISION( k, dil_level );
  } else {
    if( D != 1 && ( D%2 != 0 || (nvar/2)%(D/2) != 0 ) )
      error0("dilution %d (level %d) on a grid with %d dof per site: use 1 or 2*m with m dividing %d\n",
             D, dil_level, nvar, nvar/2);
    for( int k=0; k<nvar; k++ )
      e_tab[k] = coarse_dilution_value_PRECISION( k, nvar, D );
  }
}

// Plain (unprobed) diluted vector: the dilution pattern on every site.
// phi is a vector on level l; dil_level selects the dilution scheme.
void vector_PRECISION_dilution( vector_PRECISION phi, int start, int end, level_struct *l, int dil_level ) {
  int thread = omp_get_thread_num();
  if(thread == 0 && start != end)
    PROF_PRECISION_START( _SET );

  if ( phi != NULL ) {
    int nvar = l->num_lattice_site_var;
    int e_tab[nvar];
    dilution_table_PRECISION( e_tab, l, dil_level );

    for ( int i=start; i<end; i++ )
      phi[i] = (PRECISION)( e_tab[i%nvar] );

  } else {
    error0("Error in \"vector_PRECISION_dilution\": pointer is null\n");
  }

  if(thread == 0 && start != end)
    PROF_PRECISION_STOP( _SET, 1 );
}

// Multiplier-based probing vector: indicator of color g.coloring_count+1 on the
// grid of level l, times the dilution pattern of level dil_level.
void vector_PRECISION_probing( vector_PRECISION phi, int start, int end, level_struct *l, int dil_level ) {
  int thread = omp_get_thread_num();
  if(thread == 0 && start != end)
    PROF_PRECISION_START( _SET );

  if ( phi != NULL ) {
    int nvar = l->num_lattice_site_var;
    int e_tab[nvar];
    dilution_table_PRECISION( e_tab, l, dil_level );

    int  mycol = g.coloring_count + 1;          // colors are stored 1-based
    int *cols  = g.local_colors[l->depth];

    for ( int i=start; i<end; i++ )
      phi[i] = (PRECISION)( (cols[i/nvar] == mycol) ? e_tab[i%nvar] : 0 );

  } else {
    error0("Error in \"vector_PRECISION_probing\": pointer is null\n");
  }

  if(thread == 0 && start != end)
    PROF_PRECISION_STOP( _SET, 1 );
}

// Hierarchical-probing (Hadamard) vector number g.coloring_count, times the
// dilution pattern of level dil_level. Only available on the finest grid.
void vector_PRECISION_hadamard( vector_PRECISION phi, int start, int end, level_struct *l, int dil_level ) {
  int thread = omp_get_thread_num();
  if(thread == 0 && start != end)
    PROF_PRECISION_START( _SET );

  if( l->depth != 0 )
    error0("hierarchical probing is only implemented on the finest grid (requested at depth %d)\n", l->depth);

  if ( phi != NULL ) {
    int e_tab[12];
    dilution_table_PRECISION( e_tab, l, dil_level );

    const int       tb      = (g.probing_dimension == 4) ? g_hp_total_bits_4d : g_hp_total_bits_3d;
    const long long m       = (long long)g.coloring_count;
    const long long r       = m & ((1LL << tb) - 1);   // Hadamard column index (low bits)
    const long long q       = m >> tb;                 // base-3 Fourier digit: 0, 1, or 2
    const uint32_t *loc_tab = (g.probing_dimension == 4) ? g_hp_loc_4d : g_hp_loc_3d;
    const uint8_t  *c3_tab  = (g.probing_dimension == 4) ? g_hp_c3_4d  : g_hp_c3_3d;

    static const PRECISION re3[3] = { 1.0, -0.5, -0.5 };
    static const PRECISION im3[3] = { 0.0,  0.8660254037844386, -0.8660254037844386 };

    for ( int i=start; i<end; i++ ){
      int site = i / 12;
      int hij  = (__builtin_popcount(loc_tab[site] & (uint32_t)r) & 1) ? -1 : +1;

      complex_PRECISION val = (complex_PRECISION)(hij * e_tab[i%12]);
      if( q != 0 ){
        int ph = (int)((c3_tab[site] * q) % 3);
        val *= (re3[ph] + I*im3[ph]);
      }
      phi[i] = val;
    }
  } else {
    error0("Error in \"vector_PRECISION_hadamard\": pointer is null\n");
  }

  if(thread == 0 && start != end)
    PROF_PRECISION_STOP( _SET, 1 );
}


void vector_PRECISION_define_random_rademacher( vector_PRECISION phi, int start, int end, level_struct *l ) {

  int thread = omp_get_thread_num();
  if(thread == 0 && start != end)
  PROF_PRECISION_START( _SET );
  if ( phi != NULL ) {
    int i;
    int j = 0;

    //int size = g.global_lattice[l->depth][0]*g.global_lattice[l->depth][1]*g.global_lattice[l->depth][2]*g.global_lattice[l->depth][3];
    int dof = l->num_lattice_site_var;

    for ( i=start; i<end; i++ ){
      if(   (PRECISION)((double)rand()<(double)RAND_MAX/2.0)   ) phi[i]=  (double) (-1);
      else phi[i]= (PRECISION)(1);
    }
  } else {
    error0("Error in \"vector_PRECISION_define_random\": pointer is null\n");
  }

  if(thread == 0 && start != end)
  PROF_PRECISION_STOP( _SET, 1 );

  //print_rademacher_PRECISION( phi, l);

}


// TODO: implement separate 3D funciton, this is just for testing
/*  Takes a 4D vector and sets it’s components to 0
 * at all components except in g.time_slice
 */
void vector_PRECISION_ghg(vector_PRECISION phi,
                                               int start, int end,
                                               level_struct *l)
{
  int thread = omp_get_thread_num();
  if(thread == 0 && start != end)
    PROF_PRECISION_START( _SET );
   if(thread == 0){
    //for (int i = start; i < end; ++i) phi[i] = 0.0;
    if ( phi != NULL ) {
      int depth = l->depth;
      int r_t = g.my_coords[T];
      //int r_z = g.my_coords[Z];
      //int r_y = g.my_coords[Y];
      //int r_x = g.my_coords[X];
      //printf("rank= %d \t %d %d %d %d \n",g.my_rank, r_t, r_z, r_y, r_x);

      /* global_splitting: Num Procs in each dim
       * Nt_loc: Local Lattice sites per direcction*/
      int Nt_loc = g.global_lattice[depth][T]/l->global_splitting[T];
      int Nz_loc = g.global_lattice[depth][Z]/l->global_splitting[Z];
      int Ny_loc = g.global_lattice[depth][Y]/l->global_splitting[Y];
      int Nx_loc = g.global_lattice[depth][X]/l->global_splitting[X];

      // Initial global coordinates based on cartesian coordinates
      int t0 = r_t * Nt_loc;
      //int x0 = r_x * Nx_loc;
      //int y0 = r_y * Ny_loc;
      //int z0 = r_z * Nz_loc;

      //printf("RANK = %d \t Nt= %d,  Nz = %d, Ny= %d, Nx = %d\n T = %d, Z = %d, Y = %d, X = %d\n", g.my_rank, Nt_loc, Nz_loc, Ny_loc, Nx_loc, T,Z,Y,X);
      for (int lt = 0; lt< Nt_loc; ++lt) {
        int t = t0 + lt;               // global t
        bool ontimeslice = (t == g.time_slice);
        for (int lz = 0; lz < Nz_loc; ++lz) {
         // int z = z0 + lz;
          for (int ly = 0; ly < Ny_loc; ++ly) {
          //  int y = y0 + ly;           // global y
            for (int lx = 0; lx < Nx_loc; ++lx) {
            //  int x = x0 + lx;             // global x

              int loc_site  = ((lt * Nz_loc + lz) * Ny_loc + ly) * Nx_loc + lx;

              for( int site_el = 0; site_el < l->num_lattice_site_var; site_el++){
              //for (int d=0; d<4; ++d) {
              //  for (int c=0; c<3; ++c){
              //      int i = loc_site*12    /* site base  */
              //                + d * 3      /* stride is 3 colors  */
              //                + c;         /* color component   */
                int i = loc_site * l->num_lattice_site_var + site_el;
                if (!ontimeslice){
                  phi[i] = 0.0;
                }
              }
            }
          }
        }
      }
    }else {
      error0("Error in \"vector_PRECISION_define_random\": pointer is null\n");
    }
  }
  if(thread == 0 && start != end)
    PROF_PRECISION_STOP( _SET, 1 );
}

// Zero phi outside timeslice t of the grid of level l (t in that level's global
// coordinates). Valid on every level and site ordering (uses g.site_t).
void vector_PRECISION_timeslice( vector_PRECISION phi, int t, level_struct *l, struct Thread *threading ) {
  START_MASTER(threading)
  int nvar = l->num_lattice_site_var, *st = g.site_t[l->depth];
  for( int s=0; s<l->num_inner_lattice_sites; s++ )
    if( st[s] != t )
      for( int k=0; k<nvar; k++ ) phi[s*nvar+k] = 0.0;
  END_MASTER(threading)
  SYNC_MASTER_TO_ALL(threading)
}

void vector_PRECISION_define_spin_color( vector_PRECISION phi, int start, int end, level_struct *l, struct Thread *threading) {
  
  int thread = omp_get_thread_num();
  if(thread == 0 && start != end)
    PROF_PRECISION_START( _SET );
   if(thread == 0){
    if ( phi != NULL ) {
      int depth = l->depth;
      int r_t = g.my_coords[T];
      int r_z = g.my_coords[Z];
      int r_y = g.my_coords[Y];
      int r_x = g.my_coords[X];
      //printf("rank= %d \t %d %d %d %d \n",g.my_rank, r_t, r_z, r_y, r_x);
      
      /* global_splitting: Num Procs in each dim
       * Nt_loc: Local Lattice sites per direcction*/      
      int Nt_loc = g.global_lattice[depth][T]/l->global_splitting[T]; 
      int Nz_loc = g.global_lattice[depth][Z]/l->global_splitting[Z];
      int Ny_loc = g.global_lattice[depth][Y]/l->global_splitting[Y];
      int Nx_loc = g.global_lattice[depth][X]/l->global_splitting[X];
      
      // Initial global coordinates based on cartesian coordinates
      int t0 = r_t * Nt_loc;
      int x0 = r_x * Nx_loc;
      int y0 = r_y * Ny_loc;
      int z0 = r_z * Nz_loc;    
      
      //printf("RANK = %d \t Nt= %d,  Nz = %d, Ny= %d, Nx = %d\n T = %d, Z = %d, Y = %d, X = %d\n", g.my_rank, Nt_loc, Nz_loc, Ny_loc, Nx_loc, T,Z,Y,X);
      for (int lt=0; lt<Nt_loc; ++lt) {
        int t = t0 + lt;               // global t
        for (int lz=0; lz<Nz_loc; ++lz) {
          int z = z0 + lz;
            for (int ly=0; ly<Ny_loc; ++ly) {
              int y = y0 + ly;           // global y
                for (int lx=0; lx<Nx_loc; ++lx) {
                  int x = x0 + lx;             // global x

                  int loc_site  = ((lt * Nz_loc + lz) * Ny_loc + ly) * Nx_loc + lx;

                  for (int d=0; d<4; ++d) {
                    for (int c=0; c<3; ++c){
                        PRECISION re =  t
                                      + 3.14 * z
                                      - 2.72 * y
                                      + 0.58 * x
                                      - 1.41 * c
                                      + 1.20 * d;

                        int i = loc_site*12    /* site base  */
                                  + d * 3      /* stride is 3 colors  */
                                  + c;         /* color component   */
                        PRECISION non_re = 1.0 / (1.2345 + re); 
                        phi[i] = re + I* non_re;;
                }
              }
            }
          }
        }
      }


    }else {
      error0("Error in \"vector_PRECISION_define_random\": pointer is null\n");
    }
  }
  if(thread == 0 && start != end)
    PROF_PRECISION_STOP( _SET, 1 );
  //exit(0);
}


void print_rademacher_PRECISION( vector_PRECISION phi, level_struct *l){

  int depth = l->depth;
  int dof   = l->num_lattice_site_var;

  int Nt_loc = g.global_lattice[depth][T] / l->global_splitting[T];
  int Nz_loc = g.global_lattice[depth][Z] / l->global_splitting[Z];
  int Ny_loc = g.global_lattice[depth][Y] / l->global_splitting[Y];
  int Nx_loc = g.global_lattice[depth][X] / l->global_splitting[X];

  int r_t = g.my_coords[T];
  int r_z = g.my_coords[Z];
  int r_y = g.my_coords[Y];
  int r_x = g.my_coords[X];

  //global offset of this rank’s block
  int t0 = r_t * Nt_loc;
  int z0 = r_z * Nz_loc;
  int y0 = r_y * Ny_loc;
  int x0 = r_x * Nx_loc;

  // global lattice sizes
  int T_g = g.global_lattice[depth][0];
  int Z_g = g.global_lattice[depth][1];
  int Y_g = g.global_lattice[depth][2];
  int X_g = g.global_lattice[depth][3];


  MPI_Barrier(g.comm_cart);

  // go through global sites in order
  for (int gsite = 0; gsite < T_g * Z_g * Y_g * X_g; ++gsite) {
      // get (t,z,y,x) from lexicographic index:
      int t = gsite / (Z_g * Y_g * X_g);
      int z = (gsite / (Y_g * X_g)) % Z_g;
      int y = (gsite /  X_g)          % Y_g;
      int x =  gsite %  X_g;

      // check if this rank owns the site
      if ( t >= t0 && t < t0 + Nt_loc &&
          z >= z0 && z < z0 + Nz_loc &&
          y >= y0 && y < y0 + Ny_loc &&
          x >= x0 && x < x0 + Nx_loc )
      {
          // local indices
          int lt = t - t0;
          int lz = z - z0;
          int ly = y - y0;
          int lx = x - x0;

          int local_site = ((lt * Nz_loc + lz) * Ny_loc + ly) * Nx_loc + lx;
          int color      = g.local_colors[depth][local_site];

          // print first component of phi for this site
          printf("x = %+6.1f , color = %6d , site = %3d , t = %3d , rank = %3d\n",
                creal(phi[local_site * dof]),
                color,
                gsite,
                t,
                g.my_rank);

          fflush(stdout);
      }

      // "ordered" output
      MPI_Barrier(g.comm_cart);
  }

}

/*void vector_PRECISION_define_spin_color( vector_PRECISION phi, int start, int end, level_struct *l, struct Thread *threading) {
  
  int thread = omp_get_thread_num();
  if(thread == 0 && start != end)
    PROF_PRECISION_START( _SET );
   if(thread == 0){
    if ( phi != NULL ) {
      int i_global, i_local, site_index, owner;
      int dof_per_site = l->num_lattice_site_var; // variables per site
      
      int depth = l->depth;
      int global_sites = g.global_lattice[depth][0] *
      g.global_lattice[depth][1] *
      g.global_lattice[depth][2] *
      g.global_lattice[depth][3];
      
      //  degrees of freedom in the full global vector
      int total_dof = dof_per_site * global_sites;
      
      // Number of dofs owned by each process //TODO: is always even?
      int chunk = total_dof / g.num_processes;
      
      for (int t = 0; t < g.global_lattice[0][0]; t++) {
        for (int x = 0; x < g.global_lattice[0][3]; x++) {
          for (int y = 0; y < g.global_lattice[0][2]; y++) {
            for (int z = 0; z < g.global_lattice[0][1]; z++) {
              
              // Get the global site index at (t_slice, z, y, x)
              site_index = lex_index(t, x, y, z, g.global_lattice[0]);
              if(g.my_rank==2)printf("[%d, %d, %d, %d]\t site_index \t\t %d\n",
                          t,z,y,x, site_index);
              // Loop over the dof of that site
              int counter = 0; int color = 0; int spin = 0;
              for (i_global = dof_per_site * site_index;
                  i_global < dof_per_site * site_index + dof_per_site;
                  i_global++) {
                  
                  if(color == 3){ 
                    color = 0;
                    spin ++;
                  }
                  // Compute owner of the entry and local index 
                  owner = i_global / chunk;
                  i_local = i_global % chunk;
                  
                  // Only the owner assigns the value
                  if(owner == g.my_rank){
//START_LOCKED_MASTER(threading);
                     PRECISION re = t + x*3.14 - y*2.72 + z*0.58 - color*1.41 + spin*1.20;
                     PRECISION non_re = 1.0 / (1.2345 + re);
                     phi[i_local ] = re + I* non_re;
                    // debug print
                    //printf("[%d, %d, %d, %d]\t c = %d s= %d \t\t %d\n",
                    //     t,z,y,x,color,spin,counter);
//END_LOCKED_MASTER(threading);
//SYNC_MASTER_TO_ALL(threading);
                  }
                  
                  color++;
                  counter++;
                
              }
            }
          }
        }
      }
    }else {
      error0("Error in \"vector_PRECISION_define_random\": pointer is null\n");
    }
  }
  if(thread == 0 && start != end)
    PROF_PRECISION_STOP( _SET, 1 );
  //exit(0);
}
*/

/*for (int i_global = 0; i_global < N; i_global++) {
  owner = i_global / chunk;
  int i_local = i_global % chunk;
  if(g.my_rank==0) printf("component %d, by processor %d (local index %d)\n", i_global, owner, i_local);
}
  */  
/*void vector_PRECISION_define_random_rademacher( vector_PRECISION phi, int start, int end, level_struct *l ) {
  
  int thread = omp_get_thread_num();
  if(thread == 0 && start != end)
  PROF_PRECISION_START( _SET );
  if ( phi != NULL ) {
    int i;
    for ( i=start; i<end; i++ )
      if(   (PRECISION)((double)rand()<(double)RAND_MAX/2.0)   ) phi[i]=  (double) (-1);
      else phi[i]= (PRECISION)(1);
  } else {
    error0("Error in \"vector_PRECISION_define_random\": pointer is null\n");
  }
  if(thread == 0 && start != end)
  PROF_PRECISION_STOP( _SET, 1 );
  
}
*/
