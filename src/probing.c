#include "main.h"
#include "data_layout.h"
#include "stdbool.h"

//TODO: We need a function that frees the memory allocated by the calls of probing, there are memory leaks

int max(int v[], int size) {
    int max = v[0];

    for (int i = 1; i < size; i++) {
        if (v[i] > max) {
            max = v[i];
        }
    }

    return max;
}

int min(int v[], int size) {
    int min = v[0];

    for (int i = 1; i < size; i++) {
        if (v[i] < min) {
            min = v[i];
        }
    }

    return min;
}

int pow_int(int base, int exp) {
    int res = 1;
    for (int i = 0; i < exp; ++i) {
        res *= base;
    }
    return res;
}

int pow2_valuation(unsigned int x) {
    return __builtin_ctz(x);
}

void vector_copy(int *dest, int *src, int size) {
    for (int i = 0; i < size; i++) {
        dest[i] = src[i];
    }
}

bool contains(int *array, int size, int value) {
    for (int i = 0; i < size; i++) {
        if (array[i] == value)
            return true;
    }
    return false;
}

//TODO: move all the variance related functions here into a new file
void mlmc_connected_print_variances(){
  if(g.my_rank==0){
    for(int i=0; i<g.num_levels; i++){
        for(int j=0; j<g.num_levels; j++){
            int nlevs = g.num_levels;
            int idx = i*nlevs + j;
            printf("\n Variance of G_{%d,%d}(t=%d) = %f ", i,j,g.time_slice,g.variances[idx]);
        }
    }
  }
}

// Variances must be set to zero for each time-slice trace estimation
void set_probing_variances_to_zero(){
  if(g.my_rank == 0){
    if(g.trace_op_type==7){
    	for(int level = 0; level < g.num_levels*g.num_levels; level++){
            g.variances[level] = 0.0;
    	}

  }else{
	for(int level = 0; level < g.num_levels; level++){
            g.variances[level] = 0.0;
        }

    }
  }
}

//PRINT COLORS PER MPI RANK AND PER TIMESLICE
void print_colors(){

  char filename[100];
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int num_processes;
  MPI_Comm_size(MPI_COMM_WORLD, &num_processes);
  sprintf(filename, "print_files/print_colors/colors_timeslice_%d_rank_%d.txt", g.time_slice, rank);
  FILE *file = fopen(filename, "w");

  fprintf(file, "\nTimeslice %d, process rank %d, t = %d, z = %d, y = %d, x = %d\n", g.time_slice, rank, g.my_coords[0], g.my_coords[1], g.my_coords[2], g.my_coords[3]);

  for(int level = 0; level < g.num_levels; level++){
    fprintf(file, "\nColors at level %d\n[", level);
    int size = g.global_lattice[level][0]*g.global_lattice[level][1]*g.global_lattice[level][2]*g.global_lattice[level][3];
    int local_size = size/num_processes;
    for(int i = 0; i < local_size; i++){
      fprintf(file, " %d", g.local_colors[level][i]);
    }
    fprintf(file, "]\n");
  }
  fclose(file); 
}

void allocate_variances(){
//If we are doing mlmc with connected operator we have g.num_levels^2 operators
   if(g.my_rank==0){
      if(g.trace_op_type==7)
          MALLOC(g.variances, double, g.num_levels*g.num_levels);
      else
          MALLOC(g.variances, double, g.num_levels);

      set_probing_variances_to_zero();
   }
}


void get_coloring_dimension(){

  if(g.trace_op_type == 1 || g.trace_op_type == 2 || g.trace_op_type == 6 || g.trace_op_type == 7 || g.trace_op_type == 8)
    g.probing_dimension = 3;

  if(g.trace_op_type == 9 || g.trace_op_type == 10 || g.trace_op_type == 11 || g.trace_op_type == 12)
    g.probing_dimension = 4;
}

// number of lattice sites stored per process at a given depth
// (coarse levels may merge processes: use the level's own local lattice)
int local_colors_size( int level ){
    return g.local_lattice[level][0] * g.local_lattice[level][1]
         * g.local_lattice[level][2] * g.local_lattice[level][3];
}

void setup_local_colors_direct(){
    // calloc: entries must start as NULL for the MALLOC macro
    g.local_colors = (int**)calloc(g.num_levels, sizeof(int*));

    for(int level = 0; level < g.num_levels; level++){
        // local lattice of this level (not global/num_processes: coarse levels may merge processes)
        int Nt_loc = g.local_lattice[level][0];
        int Nz_loc = g.local_lattice[level][1];
        int Ny_loc = g.local_lattice[level][2];
        int Nx_loc = g.local_lattice[level][3];
        const int local_size = local_colors_size( level );

        MALLOC(g.local_colors[level], int, local_size);

        // global coordinate offsets of this rank's block (only used at depth 0,
        // where every process owns exactly one block of the process grid)
        int t0 = g.my_coords[0] * Nt_loc;
        int z0 = g.my_coords[1] * Nz_loc;
        int y0 = g.my_coords[2] * Ny_loc;
        int x0 = g.my_coords[3] * Nx_loc;

        int idx = 0;
        for (int lt = 0; lt < Nt_loc; lt++)
            for (int lz = 0; lz < Nz_loc; lz++)
                for (int ly = 0; ly < Ny_loc; ly++)
                    for (int lx = 0; lx < Nx_loc; lx++, idx++)
                    {
                        if( level == 0 ){
                            int t = t0 + lt;
                            int z = z0 + lz;
                            int y = y0 + ly;
                            int x = x0 + lx;
                            int col = ( g.sigma[level][0]*t + g.sigma[level][1]*z
                                      + g.sigma[level][2]*y + g.sigma[level][3]*x ) % g.num_colors[level];
                            // manual two-form override for d=7 would replace the line above:
                            //   96^3: int col = 3*(( z + 2*y + 11*x ) % 48) + ( y % 3 );
                            //   64^3: int col = 2*(( z + 3*y + 13*x ) % 64) + (( y + x ) % 2);
                            g.local_colors[level][idx] = col + 1;   // 1-based
                        } else {
                            g.local_colors[level][idx] = 1;         // coarse levels: single color
                        }
                    }
    }
}

// Multiplier-based coloring tables, selected by the geometry of each level.
// 4D tables are keyed by T x Z x Y x X, 3D tables by Z x Y x X.
typedef void (*sigma_table_fn)( int dist, int *sigma, int *nc );
typedef struct { int dim; int ext[4]; sigma_table_fn fn; } sigma_table_entry;

static const sigma_table_entry sigma_tables[] = {
  // fine lattices
  { 4, {   8,  8,  8,  8 }, get_sigma_4D_8x8x8x8         },
  { 4, {  64, 32, 32, 32 }, get_sigma_4D_64x32x32x32     },
  { 4, { 128, 64, 64, 64 }, get_sigma_4D_128x64x64x64    },
  { 4, { 192, 96, 96, 96 }, get_sigma_4D_192x96x96x96    },
  { 3, {   0,  8,  8,  8 }, get_sigma_3D_8x8x8           },
  { 3, {   0, 32, 32, 32 }, get_sigma_3D_32x32x32        },
  { 3, {   0, 64, 64, 64 }, get_sigma_3D_64x64x64        },
  { 3, {   0, 96, 96, 96 }, get_sigma_3D_96x96x96        },
  // coarse levels (placeholders, see sigma.c)
  { 4, {  16,  8,  8,  8 }, get_sigma_4D_16x8x8x8        },
  { 4, {   4,  2,  2,  2 }, get_sigma_4D_4x2x2x2         },
  { 4, {  32, 16, 16, 16 }, get_sigma_4D_32x16x16x16     },
  { 4, {   8,  4,  4,  4 }, get_sigma_4D_8x4x4x4         },
  { 4, {  48, 24, 24, 24 }, get_sigma_4D_48x24x24x24     },
  { 4, {  12,  6,  6,  6 }, get_sigma_4D_12x6x6x6        },
  { 4, {   2,  2,  2,  2 }, get_sigma_4D_2x2x2x2         },
  { 4, {   4,  4,  4,  4 }, get_sigma_4D_4x4x4x4         },
  { 3, {   0,  2,  2,  2 }, get_sigma_3D_2x2x2           },
  { 3, {   0,  4,  4,  4 }, get_sigma_3D_4x4x4           },
  { 3, {   0,  6,  6,  6 }, get_sigma_3D_6x6x6           },
  { 3, {   0, 16, 16, 16 }, get_sigma_3D_16x16x16        },
  { 3, {   0, 24, 24, 24 }, get_sigma_3D_24x24x24        },
};

// Fill g.sigma[level] and g.num_colors[level] for the coloring distance of that level
void get_sigma( int level ){
  int *n = g.global_lattice[level];
  int dist = g.coloring_dist[level];
  int found = 0;

  for( int mu=0; mu<4; mu++ ) g.sigma[level][mu] = 0;
  g.num_colors[level] = 0;

  for( int e=0; e<(int)(sizeof(sigma_tables)/sizeof(sigma_tables[0])); e++ ){
    const sigma_table_entry *t = &sigma_tables[e];
    if( t->dim != g.probing_dimension ) continue;
    int match = 1;
    for( int mu = (t->dim==4) ? 0 : 1; mu<4; mu++ )
      if( t->ext[mu] != n[mu] ) match = 0;
    if( match ){
      t->fn( dist, g.sigma[level], &(g.num_colors[level]) );
      found = 1;
      break;
    }
  }

  if( !found )
    error0("No sigma table for the %dD coloring of the %dx%dx%dx%d lattice at depth %d\n",
           g.probing_dimension, n[0], n[1], n[2], n[3], level);
  if( g.num_colors[level] == 0 )
    error0("No sigma entry for coloring distance %d of the %dx%dx%dx%d lattice at depth %d\n",
           dist, n[0], n[1], n[2], n[3], level);

  // torus validity: sigma_mu * n_mu must vanish mod nc, otherwise the coloring
  // silently violates its distance across the periodic boundary
  for( int mu = (g.probing_dimension==4) ? 0 : 1; mu<4; mu++ )
    if( ( g.sigma[level][mu] * n[mu] ) % g.num_colors[level] != 0 )
      error0("multiplier coloring INVALID on the torus at depth %d: sigma[%d]=%d, n=%d, nc=%d\n",
             level, mu, g.sigma[level][mu], n[mu], g.num_colors[level]);
}

void dilution_check(int level){

  if(g.dilution[level] != 1 && g.dilution[level] != 2 && g.dilution[level] != 3 && g.dilution[level] != 4 && g.dilution[level] != 12){
    printf("\nError: choose a correct dilution value (1, 2, 3, 4, 12) at level %d", level);
    exit(1);
  }

  if(g.dilution[level] == 1)
    if(g.my_rank==0) printf("\nNo dilution at level %d   ", level);

  if(g.dilution[level] == 2)
    if(g.my_rank==0) printf("\nPartial spin dilution at level %d   ", level);

  if(g.dilution[level] == 3)
    if(g.my_rank==0) printf("\nColor dilution at level %d   ", level);

  if(g.dilution[level] == 4)
    if(g.my_rank==0) printf("\nComplete spin dilution at level %d   ", level);

  if(g.dilution[level] == 12)
    if(g.my_rank==0) printf("\nSpin-Color dilution at level %d   ", level);
}

void coloring_scheme(){
  // sigma / num_colors selection: pure table lookup on ALL ranks, one coloring per level
  for( int level = 0; level < g.num_levels; level++ )
    get_sigma( level );

  // coarse levels are not probed yet: their probing vectors need the coarse site
  // ordering, which comes with the next step
  for( int level = 1; level < g.num_levels; level++ )
    g.num_colors[level] = 1;

  if(g.my_rank==0){
    printf("\nProbing = %d - Classical probing\n", g.probing);
    printf("Grids to be colored = %d\n", g.colored_grids);
    printf("Coloring dimension = %d\n", g.probing_dimension);
    for( int level = 0; level < g.num_levels; level++ )
      printf("Depth %d: coloring distance = %d, sigma = %d %d %d %d\n", level, g.coloring_dist[level],
             g.sigma[level][0], g.sigma[level][1], g.sigma[level][2], g.sigma[level][3]);
  }

  for( int level = 0; level < g.num_levels; level++ )
    dilution_check(level);

  setup_local_colors_direct();

  if(g.my_rank==0)
    for( int level = 0; level < g.num_levels; level++ )
      printf("\n Colors at depth %d : \t %d ", level, g.num_colors[level]);

  MPI_Barrier(MPI_COMM_WORLD);
}

int* find_indices(int *array, int size, int value, int *count){
    // first pass: count how many matches
    *count = 0;
    for(int i = 0; i < size; i++){
        if(array[i] == value)
            (*count)++;
    }
    // allocate result
    int *indices = (int*)malloc(*count * sizeof(int));
    if(indices == NULL)
        return NULL;
    // second pass: fill indices
    int j = 0;
    for(int i = 0; i < size; i++){
        if(array[i] == value)
            indices[j++] = i;
    }
    return indices;
}

void graph_coloring(){

  get_coloring_dimension();
  if(g.probing == 1) coloring_scheme();

 if(g.probing == 2){
   if(g.my_rank==0) printf("\nProbing = %d - Hierarchical probing\n", g.probing);
   if(g.my_rank==0) printf("Applied to num levels = %d\n", g.colored_grids);
   if(g.my_rank==0) printf("Probing dimension = %d\n", g.probing_dimension);
   for(int level = 0; level<g.num_levels; level++){
     g.global_k[level][0] = pow2_valuation(g.global_lattice[level][0]);
     g.global_k[level][1] = pow2_valuation(g.global_lattice[level][1]);
     g.global_k[level][2] = pow2_valuation(g.global_lattice[level][2]);
     g.global_k[level][3] = pow2_valuation(g.global_lattice[level][3]);
     dilution_check(level);
     g.num_colors[level] = g.n_had[level];

     if(g.my_rank==0){
       printf("Number of Hadamard vectors at level %d = %d\n", level, g.num_colors[level]);
       printf("Global k_t at level %d = %d\n", level, g.global_k[level][0]);
       printf("Global k_z at level %d = %d\n", level, g.global_k[level][1]);
       printf("Global k_y at level %d = %d\n", level, g.global_k[level][2]);
       printf("Global k_x at level %d = %d\n", level, g.global_k[level][3]);
       printf("Hierarchy at level %d = %d\n", level, g.k[level]);
       printf("Anisotropy at level %d = %d\n", level, g.anisotropic[level]);
     }
   }
//MOVE THIS CHECK INSIDE THE PER-LEVEL LOOP WHEN DOING MGMLMC
   int tb = (g.probing_dimension == 4)
             ? g.global_k[0][0] + 3*g.global_k[0][1]
             : 3*g.global_k[0][1];
   long long limit = 3LL * (1LL << tb);
   if( (long long)g.num_colors[0] > limit )
     error0("HP: num_colors %d exceeds 3*2^%d = %lld (one three-coloring level past the bit levels); undefined beyond\n",
             g.num_colors[0], tb, limit);

 }
 
 if(g.probing == 0){
   // plain (diluted) Hutchinson: one "color" per level, no probing vectors
   if(g.my_rank==0) printf("\nProbing = %d - Plain estimator\n", g.probing);
   for(int level = 0; level<g.num_levels; level++){
     dilution_check(level);
     g.num_colors[level] = 1;
   }
 }
}
