#include "main.h"

void get_sigma_4D_8x8x8x8( int dist, int *sigma, int *nc ){

  if(dist == 1){
    sigma[0] = 1;
    sigma[1] = 1;
    sigma[2] = 1;
    sigma[3] = 1;

    *nc = 2;
  }

  if(dist == 2){
    sigma[0] = 1;
    sigma[1] = 2;
    sigma[2] = 3;
    sigma[3] = 4;

    *nc = 10;
  }
}

void get_sigma_3D_8x8x8( int dist, int *sigma, int *nc ){

  if(dist == 1){
    sigma[0] = 1;
    sigma[1] = 1;
    sigma[2] = 1;
    sigma[3] = 1;

    *nc = 2;
  }

  if(dist == 2){
    sigma[0] = 1;
    sigma[1] = 2;
    sigma[2] = 3;
    sigma[3] = 4;

    *nc = 8;
  }
}

void get_sigma_4D_64x32x32x32( int dist, int *sigma, int *nc ){

  if(dist == 1){
    sigma[0] = 1;
    sigma[1] = 1;
    sigma[2] = 1;
    sigma[3] = 1;

    *nc = 2;
  }

  if(dist == 2){
    sigma[0] = 1;
    sigma[1] = 2;
    sigma[2] = 3;
    sigma[3] = 4;

    *nc = 10;
  }

  if(dist == 3){
    sigma[0] = 1;
    sigma[1] = 5;
    sigma[2] = 55;
    sigma[3] = 61;

    *nc = 16;
  }

  if(dist == 4){
    sigma[0] = 1;
    sigma[1] = 8;
    sigma[2] = 12;
    sigma[3] = 18;

    *nc = 64;
  }

    if(dist == 5){
    sigma[0] = 38;
    sigma[1] = 1;
    sigma[2] = 12;
    sigma[3] = 16;

    *nc = 128;
  }

  if(dist == 6){
    sigma[0] = 3;
    sigma[1] = 20;
    sigma[2] = 48;
    sigma[3] = 50;

    *nc = 320;
  }

  if(dist == 7){
    sigma[0] = 40;
    sigma[1] = 32;
    sigma[2] = 33;
    sigma[3] = 61;

    *nc = 512;
  }

}

void get_sigma_3D_32x32x32( int dist, int *sigma, int *nc ){

  if(dist == 1){
    sigma[0] = 0;
    sigma[1] = 1;
    sigma[2] = 1;
    sigma[3] = 1;

    *nc = 2;
  }

  if(dist == 2){
    sigma[0] = 0;
    sigma[1] = 1;
    sigma[2] = 2;
    sigma[3] = 3;

    *nc = 8;
  }

  if(dist == 3){
    sigma[0] = 0;
    sigma[1] = 1;
    sigma[2] = 3;
    sigma[3] = 5;

    *nc = 16;
  }

  if(dist == 4){
    sigma[0] = 0;
    sigma[1] = 1;
    sigma[2] = 6;
    sigma[3] = 9;

    *nc = 32;
  }

    if(dist == 5){
    sigma[0] = 0;
    sigma[1] = 1;
    sigma[2] = 11;
    sigma[3] = 27;

    *nc = 88;
  }

  if(dist == 6){
    sigma[0] = 0;
    sigma[1] = 1;
    sigma[2] = 8;
    sigma[3] = 44;

    *nc = 128;
  }

  if(dist == 7){
    sigma[0] = 0;
    sigma[1] = 1;
    sigma[2] = 9;
    sigma[3] = 33;

    *nc = 176;
  }

  if(dist == 8){
    sigma[0] = 0;
    sigma[1] = 7;
    sigma[2] = 48;
    sigma[3] = 51;

    *nc = 272;
  }

  if(dist == 9){
    sigma[0] = 0;
    sigma[1] = 1;
    sigma[2] = 33;
    sigma[3] = 45;

    *nc = 352;
  }

}

void get_sigma_4D_128x64x64x64( int dist, int *sigma, int *nc ){

  if(dist == 1){
    sigma[0] = 1;
    sigma[1] = 1;
    sigma[2] = 1;
    sigma[3] = 1;

    *nc = 2;
  }

  if(dist == 2){
    sigma[0] = 1;
    sigma[1] = 2;
    sigma[2] = 3;
    sigma[3] = 4;

    *nc = 16;
  }

  if(dist == 3){
    sigma[0] = 1;
    sigma[1] = 3;
    sigma[2] = 5;
    sigma[3] = 7;

    *nc = 16;
  }

  if(dist == 4){
    sigma[0] = 1;
    sigma[1] = 4;
    sigma[2] = 10;
    sigma[3] = 17;

    *nc = 64;
  }

    if(dist == 5){
    sigma[0] = 1;
    sigma[1] = 6;
    sigma[2] = 16;
    sigma[3] = 52;

    *nc = 128;
  }

  if(dist == 6){
    sigma[0] = 6;
    sigma[1] = 1;
    sigma[2] = 48;
    sigma[3] = 120;

    *nc = 256;
  }

  if(dist == 7){
    sigma[0] = 1;
    sigma[1] = 7;
    sigma[2] = 33;
    sigma[3] = 89;

    *nc = 448;
  }

  if(dist == 8){
    sigma[0] = 4;
    sigma[1] = 26;
    sigma[2] = 114;
    sigma[3] = 117;

    *nc = 832;
  }

}

void get_sigma_3D_64x64x64( int dist, int *sigma, int *nc ){

  if(dist == 1){
    sigma[0] = 0;
    sigma[1] = 1;
    sigma[2] = 3;
    sigma[3] = 5;

    *nc = 2;
  }

  if(dist == 2){
    sigma[0] = 0;
    sigma[1] = 1;
    sigma[2] = 2;
    sigma[3] = 3;

    *nc = 8;
  }

  if(dist == 3){
    sigma[0] = 0;
    sigma[1] = 1;
    sigma[2] = 3;
    sigma[3] = 5;

    *nc = 16;
  }

  if(dist == 4){
    sigma[0] = 0;
    sigma[1] = 1;
    sigma[2] = 4;
    sigma[3] = 10;

    *nc = 32;
  }

    if(dist == 5){
    sigma[0] = 0;
    sigma[1] = 1;
    sigma[2] = 5;
    sigma[3] = 13;

    *nc = 64;
  }

  if(dist == 6){
    sigma[0] = 0;
    sigma[1] = 1;
    sigma[2] = 6;
    sigma[3] = 26;

    *nc = 128;
  }

  if(dist == 7){
    sigma[0] = 0;
    sigma[1] = 1;
    sigma[2] = 7;
    sigma[3] = 57;

    *nc = 152;
  }

  if(dist == 8){
    sigma[0] = 0;
    sigma[1] = 1;
    sigma[2] = 36;
    sigma[3] = 52;

    *nc = 256;
  }

  if(dist == 9){
    sigma[0] = 0;
    sigma[1] = 4;
    sigma[2] = 30;
    sigma[3] = 35;

    *nc = 320;
  }

}


void get_sigma_4D_192x96x96x96( int dist, int *sigma, int *nc ){

  if(dist == 1){
    sigma[0] = 1;
    sigma[1] = 1;
    sigma[2] = 1;
    sigma[3] = 1;

    *nc = 2;
  }

  if(dist == 2){
    sigma[0] = 1;
    sigma[1] = 2;
    sigma[2] = 3;
    sigma[3] = 4;

    *nc = 12;
  }

  if(dist == 3){
    sigma[0] = 1;
    sigma[1] = 3;
    sigma[2] = 5;
    sigma[3] = 7;

    *nc = 16;
  }

  if(dist == 4){
    sigma[0] = 1;
    sigma[1] = 4;
    sigma[2] = 11;
    sigma[3] = 17;

    *nc = 48;
  }

    if(dist == 5){
    sigma[0] = 1;
    sigma[1] = 5;
    sigma[2] = 13;
    sigma[3] = 29;

    *nc = 96;
  }

  if(dist == 6){
    sigma[0] = 1;
    sigma[1] = 6;
    sigma[2] = 50;
    sigma[3] = 76;

    *nc = 192;
  }

}

void get_sigma_3D_96x96x96( int dist, int *sigma, int *nc ){

  if(dist == 1){
    sigma[0] = 0;
    sigma[1] = 1;
    sigma[2] = 1;
    sigma[3] = 1;

    *nc = 2;
  }

  if(dist == 2){
    sigma[0] = 0;
    sigma[1] = 5;
    sigma[2] = 1;
    sigma[3] = 6;

    *nc = 8;
  }

  if(dist == 3){
    sigma[0] = 0;
    sigma[1] = 7;
    sigma[2] = 3;
    sigma[3] = 1;

    *nc = 12;
  }

  if(dist == 4){
    sigma[0] = 0;
    sigma[1] = 28;
    sigma[2] = 31;
    sigma[3] = 19;

    *nc = 32;
  }

    if(dist == 5){
    sigma[0] = 0;
    sigma[1] = 41;
    sigma[2] = 23;
    sigma[3] = 13;

    *nc = 48;
  }

  if(dist == 6){
    sigma[0] = 0;
    sigma[1] = 21;
    sigma[2] = 65;
    sigma[3] = 30;

    *nc = 96;
  }

  if(dist == 7){
    sigma[0] = 0;
    sigma[1] = 42;
    sigma[2] = 49;
    sigma[3] = 117;

    *nc = 144;
  }

}

// ---------------------------------------------------------------------------
// Coarse-level geometries (4^4 aggregation of 64x32^3, 128x64^3, 192x96^3 and
// the 8^4 test lattices). PLACEHOLDERS: every coloring distance currently gets
// the distance-1 (Red-Black) coloring sigma = 1, nc = 2. Replace with real
// tables once computed, using the same format as the fine-lattice tables above.
// ---------------------------------------------------------------------------

void get_sigma_4D_16x8x8x8( int dist, int *sigma, int *nc ){

  if(dist >= 1){   // placeholder: distance-1 coloring for every distance
    sigma[0] = 1;
    sigma[1] = 1;
    sigma[2] = 1;
    sigma[3] = 1;

    *nc = 2;
  }
}

void get_sigma_4D_4x2x2x2( int dist, int *sigma, int *nc ){

  if(dist >= 1){   // placeholder: distance-1 coloring for every distance
    sigma[0] = 1;
    sigma[1] = 1;
    sigma[2] = 1;
    sigma[3] = 1;

    *nc = 2;
  }
}

void get_sigma_4D_32x16x16x16( int dist, int *sigma, int *nc ){

  if(dist >= 1){   // placeholder: distance-1 coloring for every distance
    sigma[0] = 1;
    sigma[1] = 1;
    sigma[2] = 1;
    sigma[3] = 1;

    *nc = 2;
  }
}

void get_sigma_4D_8x4x4x4( int dist, int *sigma, int *nc ){

  if(dist >= 1){   // placeholder: distance-1 coloring for every distance
    sigma[0] = 1;
    sigma[1] = 1;
    sigma[2] = 1;
    sigma[3] = 1;

    *nc = 2;
  }
}

void get_sigma_4D_48x24x24x24( int dist, int *sigma, int *nc ){

  if(dist >= 1){   // placeholder: distance-1 coloring for every distance
    sigma[0] = 1;
    sigma[1] = 1;
    sigma[2] = 1;
    sigma[3] = 1;

    *nc = 2;
  }
}

void get_sigma_4D_12x6x6x6( int dist, int *sigma, int *nc ){

  if(dist >= 1){   // placeholder: distance-1 coloring for every distance
    sigma[0] = 1;
    sigma[1] = 1;
    sigma[2] = 1;
    sigma[3] = 1;

    *nc = 2;
  }
}

void get_sigma_4D_2x2x2x2( int dist, int *sigma, int *nc ){

  if(dist >= 1){   // placeholder: distance-1 coloring for every distance
    sigma[0] = 1;
    sigma[1] = 1;
    sigma[2] = 1;
    sigma[3] = 1;

    *nc = 2;
  }
}

void get_sigma_4D_4x4x4x4( int dist, int *sigma, int *nc ){

  if(dist >= 1){   // placeholder: distance-1 coloring for every distance
    sigma[0] = 1;
    sigma[1] = 1;
    sigma[2] = 1;
    sigma[3] = 1;

    *nc = 2;
  }
}

void get_sigma_3D_2x2x2( int dist, int *sigma, int *nc ){

  if(dist >= 1){   // placeholder: distance-1 coloring for every distance
    sigma[0] = 0;
    sigma[1] = 1;
    sigma[2] = 1;
    sigma[3] = 1;

    *nc = 2;
  }
}

void get_sigma_3D_4x4x4( int dist, int *sigma, int *nc ){

  if(dist >= 1){   // placeholder: distance-1 coloring for every distance
    sigma[0] = 0;
    sigma[1] = 1;
    sigma[2] = 1;
    sigma[3] = 1;

    *nc = 2;
  }
}

void get_sigma_3D_6x6x6( int dist, int *sigma, int *nc ){

  if(dist >= 1){   // placeholder: distance-1 coloring for every distance
    sigma[0] = 0;
    sigma[1] = 1;
    sigma[2] = 1;
    sigma[3] = 1;

    *nc = 2;
  }
}

void get_sigma_3D_16x16x16( int dist, int *sigma, int *nc ){

  if(dist >= 1){   // placeholder: distance-1 coloring for every distance
    sigma[0] = 0;
    sigma[1] = 1;
    sigma[2] = 1;
    sigma[3] = 1;

    *nc = 2;
  }
}

void get_sigma_3D_24x24x24( int dist, int *sigma, int *nc ){

  if(dist >= 1){   // placeholder: distance-1 coloring for every distance
    sigma[0] = 0;
    sigma[1] = 1;
    sigma[2] = 1;
    sigma[3] = 1;

    *nc = 2;
  }
}
