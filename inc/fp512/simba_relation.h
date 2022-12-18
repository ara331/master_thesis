#ifndef _SIMBA_PARAMETERS_H_
#define _SIMBA_PARAMETERS_H_

// SIMBA-(NUMBER_OF_BATCHES)-MY
#define NUMBER_OF_BATCHES 3
//#define MY 8
#define MY 0
//relation lattice
static uint64_t Number_of_isogenies_2[N]={
    159, 174, 177, 173, 177, 193, 172, 197, 196, 175, 219, 200, 184, 196, 210, 210, 210, 233, 241, 216, 238, 205, 249, 229, 226, 241, 197, 226, 210, 200, 224, 234, 243, 259, 230, 235, 245, 223, 231, 235, 281, 200, 242, 241, 242, 204, 243, 205, 253, 261, 230, 280, 225, 206, 239, 244, 231, 255, 240, 235, 255, 248, 228, 227, 264, 234, 247, 208, 258, 233, 243, 231, 245, 263
};
static int8_t Relation_lattice[N][N] ={{-8, 0, -3, -1, 0, -1, 6, -1, 2, 2, -2, -1, 6, 0, 4, -1, -5, 1, -2, 0, -2, 2, 0, -2, -2, 0, 1, 4, 0, 2, 1, -3, 2, -7, 6, -1, 2, 0, 1, 0, 3, -4, 3, 3, 2, -2, 4, 4, -2, -2, 0, -3, 2, 0, 1, -2, 0, -2, -3, -1, 2, 1, 1, -1, -1, -2, 2, -5, 3, 2, 4, -3, -4, -1}, {-2, -4, -2, 3, -2, -2, 3, 1, -1, 3, 1, 2, -9, 2, 3, -4, 0, -1, -2, -1, 3, -2, -3, 2, 5, -2, 2, -7, -1, -3, -5, 2, -1, 5, 1, 3, -4, -1, 4, -5, -3, 1, 1, 2, 1, 3, 1, -4, -1, 1, -4, 2, -1, 2, 0, -1, 0, -1, -1, 2, -1, -1, -1, 0, 4, -1, -4, 0, -7, 1, -5, 2, 2, -4}, {1, 3, 0, -2, -4, 1, -7, 0, 3, -4, 0, 0, -4, 6, 4, -1, -2, 3, -1, -2, 2, -4, -2, 3, 3, 3, 0, -3, -1, 1, -1, -1, 3, 4, 0, -1, 4, -3, 1, -1, -1, 3, -4, 3, 0, -6, 4, -7, -2, 1, 2, 2, 2, -1, -3, 2, -3, 5, 5, 0, 0, -1, -2, 3, 6, 2, -2, -3, -4, -3, -4, 0, -2, 0}, {2, 1, -1, -3, 3, -1, -1, 3, -2, 1, 4, 2, 1, -1, -8, -3, 2, 2, 6, 0, -3, -1, 3, 2, -8, -3, 0, 5, -1, 4, -3, 2, 3, -4, 1, -3, 0, 2, -1, 0, 5, 0, 5, -3, 1, 3, -3, 3, 7, 0, 4, -4, -6, 0, 0, 0, 1, 0, 2, 1, 0, 0, 0, 1, -6, 1, 4, 0, 2, -1, -1, 2, -5, -6}, {-1, -3, 0, 2, 1, 5, 0, 4, -1, 1, 3, 3, 1, 1, 7, -7, -3, 2, 5, -1, -2, -1, -4, -2, -4, 1, 0, 2, 2, 2, 1, -2, 1, -5, 2, -6, 3, 2, 5, -3, -1, -3, 6, 1, 3, 0, -3, -3, 0, 2, 0, -4, 2, -1, 2, -2, 4, 0, -1, 0, 5, -2, -5, 1, -3, 4, -1, -6, 1, -2, 0, -5, 1, -2}, {6, -3, 2, -2, -2, 4, -6, -1, 5, -2, 0, 6, 1, -1, -1, 0, 5, -2, 5, -2, -1, -2, -1, -6, -2, 7, -4, -5, -1, 0, 9, 2, -1, -1, -2, 1, -3, 1, 0, 6, -2, 0, 0, 2, 1, 1, -2, -2, -1, 0, -2, 1, 5, -1, -1, 2, 3, 0, -2, 2, 0, -1, -7, -4, -4, 1, -3, 2, 8, -1, 3, -2, 9, 7}, {0, 3, 2, 3, 1, 0, 1, 1, -2, 0, -1, 3, 2, -4, 1, -3, -2, 2, -2, -4, 2, -4, 1, 7, 5, 2, 0, 5, 4, 3, -1, -1, -1, 2, -2, -2, -3, -2, -3, 1, 0, 2, -3, 1, -5, 0, -1, -2, 2, 1, -2, 5, -1, 0, -2, 1, -5, 2, 6, 0, -1, 0, 2, 3, 6, 10, -1, 0, -8, -2, -7, -1, -1, -1}, {6, 0, 0, 0, 2, 3, 0, -1, -3, -2, -1, -3, 9, -2, 7, 1, -1, 2, -2, 3, -4, 4, -4, -4, -1, 1, 1, 5, 1, 2, -1, -4, -4, -1, -2, -3, 4, 3, -2, 2, -6, 3, -1, -3, 6, 0, -7, 4, -6, -2, 1, 0, -4, -2, 2, 0, 0, -7, -2, -1, 4, -4, 4, 1, -3, 0, 1, -3, 4, -7, 2, -3, 3, 0}, {-4, -7, 1, 0, -4, -4, -2, 5, 3, 2, 3, 0, -2, 3, 0, -1, -4, 5, 4, 0, -1, -2, -4, 2, -1, -7, 0, -5, -3, -1, 0, -1, 7, 2, 4, 2, 5, 0, -2, -3, 6, -5, -2, 0, 1, -4, 0, -4, 0, 3, 5, -2, 0, -3, 3, 3, 0, 3, -2, -3, 1, 1, 1, 1, 0, -3, 1, -5, 1, 6, 3, 6, -7, -5}, {-6, -6, 4, 5, -4, -1, -4, 0, -1, 4, -2, 2, 0, 2, -1, 1, 1, 3, -4, 1, 2, 5, 0, 0, 5, -1, -2, -3, -3, 2, 3, -1, 2, 0, 5, 1, 3, -4, 3, 0, 6, -7, 0, -1, -5, -1, -3, 0, -4, -1, 4, 2, 4, 3, -1, 1, -2, 6, 1, -5, -1, -1, -2, -1, 2, 2, 0, 2, -2, 4, 0, 1, 1, -2}, {-4, 2, 4, 4, 5, 0, 1, 2, -7, 2, -2, -5, -7, -3, 2, 5, 3, -3, -6, 3, 4, 2, 3, 2, 7, 4, 1, 1, 0, 3, -1, 0, -4, 2, -1, -4, -3, 0, 3, 1, -2, -2, 1, -2, -5, 3, 1, 2, -3, -1, -6, 2, 4, 1, 2, 5, -1, 0, 3, 4, -5, -1, 1, 3, 4, 10, -3, 11, -6, 1, -4, -1, 0, 3}, {-6, -3, -3, 4, 0, -3, 4, 4, -2, 1, -2, 0, 1, -4, -2, -3, 4, 7, 1, 1, -3, 0, -3, -1, 0, -4, 1, 1, 1, 3, -4, -4, 2, -3, 4, 1, -5, 4, -6, 5, 0, -7, 6, -1, 1, 5, -2, 4, 3, -4, 1, 1, -9, -2, 0, 8, 0, 1, -4, -2, -3, -2, 2, 0, 0, -2, 4, -1, -3, 1, 0, 1, -2, -8}, {-4, 0, 2, 3, 0, 2, -3, -7, -1, 0, -6, 4, 0, 2, -3, 0, -1, 6, 0, 1, -3, 5, -2, -7, -2, -4, -3, 3, 5, -2, -2, 0, 2, -2, 4, 2, 2, 0, 2, 5, 2, -7, -1, -7, -3, -2, 2, -4, -2, -4, 2, -2, -2, -1, -3, 1, -1, 4, -1, -6, 1, 2, -2, -2, 1, -2, 1, -5, 0, 1, 1, -3, 0, -4}, {5, -6, 0, -5, 1, -1, 1, -1, 1, 3, 3, -2, 0, 4, -7, -3, 3, -5, 0, 2, 2, 4, 2, 1, -1, 3, -1, 4, -5, -7, 0, 1, -1, -3, -2, 6, 2, 0, 1, -3, 2, 6, 0, 1, 2, 2, -6, 6, 3, 3, 5, -2, -2, 3, -2, -6, 1, 4, -6, 1, -3, 1, -3, 0, -2, -4, 1, 5, -1, -1, -1, 5, 2, 2}, {1, -10, -4, -3, -2, 0, 3, 5, 7, 5, -1, 3, 6, -3, 5, -5, 2, 1, -2, -2, 1, 0, -2, 6, 0, 1, 3, -3, -8, -1, 5, 3, -1, 0, 8, 5, -4, 5, 2, 1, -2, 1, 0, 3, 0, 4, -5, 5, -2, 0, -1, 4, 1, 4, 0, -5, -5, -3, -6, 1, 2, -4, -1, -3, 2, -1, -5, -5, 0, 0, 2, 1, -1, -2}, {-4, 5, -3, 0, 1, -1, -4, -3, 2, -2, -2, 3, -3, 0, -5, 0, 4, 5, 1, -2, 5, -3, 5, 10, -1, -2, 1, -2, -2, 4, 0, 3, 4, 5, 9, 3, -4, -1, 0, 5, 5, -3, 1, -1, -4, 4, 4, 2, 8, 1, -1, 1, -1, 6, -2, 7, -2, 1, 4, -1, -2, 0, -2, -3, 2, 7, 1, 0, -4, 3, 4, 0, 1, -3}, {2, -1, 0, 5, 2, 0, 1, -4, -3, -3, 10, -1, -3, 3, -2, 1, 6, -3, -2, 3, 0, -1, 0, 6, 0, 2, 2, 1, 5, -1, -5, -3, -3, 5, -5, 0, 0, 1, 0, -7, -1, 5, 5, 2, 2, 3, -5, 3, 2, 2, 6, 2, -11, 5, 1, 2, 1, -2, 3, -2, -1, -6, 0, 3, -2, -5, 5, 0, -5, -3, -6, 0, 1, -6}, {-7, 1, 2, 3, 1, 4, 0, -2, 5, 0, -4, -4, 6, -2, 5, 0, -1, 7, -9, 1, 6, 4, 3, 6, 1, 3, 1, 3, 8, -2, 4, -3, 1, 4, 0, -2, 5, 3, -3, 6, -3, -2, -4, 5, -5, -2, 5, 2, 0, 2, 6, 2, 4, 2, -5, -3, -5, 3, 2, -4, -2, 0, 2, 1, 4, 2, 0, -8, -7, 1, 3, -3, 1, -1}, {-7, -5, 2, 7, 2, 1, 6, -11, -4, 0, -1, 6, -4, 3, 1, -2, 3, -7, -2, 1, 0, 4, -2, -2, -3, -3, 1, 2, 2, 0, -5, -2, 2, -9, 8, 3, -3, 0, 10, -4, -3, -3, 4, -4, 0, 0, -2, 3, 1, -2, -1, -6, -4, 2, 4, -2, 4, 2, -2, -5, 5, 1, -9, 0, -1, -6, 2, -5, -3, 1, 0, -4, 1, -9}, {2, 2, 0, 6, 1, -3, -6, -1, -11, -4, 1, 3, -5, -3, -8, 4, 6, 0, 4, 1, -2, 5, 1, -3, 2, -10, 2, -4, 4, 3, -8, -2, 2, -2, -4, -4, 6, -5, 0, -2, 0, -2, 0, -9, 1, 1, -4, -1, 6, 3, 1, 2, -3, -2, 5, -1, 2, 0, 0, -4, 1, 3, -1, 0, 2, 0, 6, 1, -1, -2, 1, 2, 3, -4}, {-3, 7, -3, 0, 1, 0, -3, 0, 0, 2, -9, -9, 6, -3, 11, 0, -8, 9, -5, 3, 0, 0, 1, 1, 5, 0, 1, 4, 1, -1, 2, -4, 2, -1, 2, -3, 6, 3, -1, 11, -3, -1, -8, -4, 1, -2, 5, 1, -7, -2, 3, -3, 0, 2, -7, 4, -4, -4, -2, 1, -2, -4, 8, -4, 3, 4, 2, 0, -5, -1, 7, -1, -1, 1}, {-2, -1, 3, -2, 0, -2, -1, -5, -3, 2, -4, 3, -2, 4, 9, 0, -5, -1, -3, 0, -2, -3, -5, -5, 3, 0, 2, 0, 3, 1, -6, -1, 3, 13, 0, -1, 1, -1, -7, -2, -2, 0, -5, 0, 8, 0, 2, -5, -5, 1, -6, 2, 2, -5, 4, 8, 1, -2, 1, 0, 2, -1, 8, 4, 1, 3, -5, 1, -1, -3, 3, -1, 1, 1}, {1, 8, 4, 0, 2, -1, -3, 0, 4, -3, 5, -3, 0, 2, -2, 5, -1, -3, 4, 1, -4, -7, 5, 7, -1, -2, -1, 5, 5, 5, 6, 3, 6, 7, -6, -3, -3, 6, -1, -1, 6, 4, -4, -3, -6, 2, 5, -4, 6, 3, 0, -2, -1, -4, -1, 7, 2, 4, 14, -1, -1, -3, 1, 5, 1, 4, 1, 3, -6, 0, 1, 0, -6, -2}, {0, -3, -4, -5, 2, -3, -1, 7, 8, 3, 0, -6, -11, 4, 4, -3, 1, 8, 5, 6, 0, 3, -1, -4, -1, 0, -3, -7, 0, -4, 1, 2, 3, 2, 1, 4, -1, 2, 5, 2, 3, -2, 1, 3, -2, -2, 4, -7, 1, 4, 3, -3, 0, -3, -6, -2, 1, 3, -3, 5, -6, -1, -8, 0, 6, -2, -6, 1, -3, -3, 1, 1, -1, -2}, {4, -1, 2, 4, 1, 5, -10, 2, -2, -6, 7, 2, 9, 1, -6, 0, 3, 1, 1, 0, -2, 1, 2, 3, 0, -3, -1, 5, 1, 5, 2, 2, 3, 5, -3, -7, 1, 1, 0, -6, 0, -1, -1, -3, -4, 3, -8, -3, 3, 2, -1, 4, 4, 1, 1, -4, -5, 5, 7, -6, 5, 4, 1, 6, 3, 4, 6, -3, 1, 1, -3, 0, -1, -2}, {-2, 7, -1, 1, -2, -2, 1, 4, 0, 3, -2, -1, -8, 0, -7, -2, -1, 4, 4, -4, 1, -2, -1, 1, 0, -7, 3, -1, 2, 0, -3, 2, -3, 5, -9, -3, -5, 2, -6, -4, 2, 1, 0, 0, -3, 3, 4, -4, 6, 1, -2, 3, -3, -8, 0, -3, 2, 3, 4, 1, -2, 7, 3, 5, 3, 1, -3, -4, -8, 8, -9, 11, -1, -7}, {1, -1, 7, 4, 2, -3, 7, 1, -5, 3, -6, 3, 2, -7, 0, -1, 0, 0, 5, -1, -2, -3, -5, -5, 2, -7, -1, -3, 3, 0, 1, 0, 2, -4, -2, 0, 2, -7, -2, 7, -1, -8, 0, -1, -3, 1, 0, -3, 3, 6, -2, -5, -1, -2, 5, 2, 5, -2, -3, 0, 4, -1, 0, -2, -3, 0, 0, 1, 2, -2, 7, -1, 0, 4}, {-4, -1, -5, -2, -1, -1, -1, 5, 0, 1, -5, -5, -4, -2, 10, -3, -1, 7, -1, 0, 1, -7, 1, -1, 4, 2, -1, -6, 1, -4, 1, 1, -2, 0, 3, 3, 0, -3, 1, 3, -6, -7, 1, 11, -4, -5, 7, -3, -5, 1, -5, -1, 5, 1, -5, 0, -5, -9, -3, 5, -4, 1, -1, -4, 3, 0, -2, -3, 0, -4, 2, -3, 1, 5}, {2, 10, 5, 3, -3, 0, 1, -7, 5, 0, -7, 5, 0, 2, -5, 3, 3, -1, -4, -3, -5, 0, -1, -2, 2, 4, 1, -1, 1, -2, 2, 4, -1, 8, -5, -4, -3, -1, 1, -2, 1, 0, -1, -2, -2, 7, 0, -5, -6, -4, -2, 5, -4, 1, -4, 2, 1, 4, 6, 3, -1, -4, 3, 0, -3, -5, 0, 0, -1, 2, -7, 4, -1, 0}, {3, 2, 4, 4, 0, 2, 1, -4, -1, -5, -4, -2, 1, -2, 1, 3, 0, 0, -7, -1, 2, 0, 8, 2, -4, -3, 0, 0, -2, -5, -2, 2, 2, 0, 2, 3, -2, -2, 10, -3, -8, -3, -3, -3, -3, -5, -3, -1, -3, -8, -4, -6, 0, 7, -2, -3, 1, -3, 1, -2, 3, -1, -4, -4, 0, -4, 3, 1, 1, 0, -2, -4, 2, 1}, {-2, 1, 0, -5, 0, 1, -10, 0, -8, -3, 8, -6, -3, 6, -10, -1, -1, 1, 2, 1, -2, -1, 5, -2, -1, -2, 0, 6, 2, 0, 0, -3, 0, 3, -1, -1, 8, -3, -1, -8, 6, 5, 0, -2, 0, -5, -2, 1, 2, 2, 3, -1, 2, 3, -4, -5, -6, 5, 4, -6, 0, 2, -1, 1, 6, 1, 7, 7, -1, -1, -6, 2, 3, 4}, {2, 0, 3, -1, 3, 6, -2, 3, 1, 4, 3, -1, -1, -2, 4, -1, 1, -5, -2, 0, -1, 1, 3, -3, -3, 1, 3, -5, 7, -4, 10, 3, -6, 11, -5, -1, -2, 2, -3, -10, 2, 2, 2, 3, 0, 3, -3, -3, -3, 0, -2, 0, 5, -1, 0, -8, 3, -2, 1, -4, -1, 0, 0, -2, -8, -7, -7, -5, 7, 2, -2, 0, 8, 9}, {3, 6, -2, 3, 2, 5, -4, -1, 2, -3, 0, 5, -5, 1, 1, 0, 0, -1, -6, 1, 8, 1, 3, 8, -7, 5, 7, -3, 1, 6, -9, -1, 1, 6, 3, -2, -2, 3, 0, 0, -8, 5, 5, -2, 4, 3, 1, 0, 9, -4, 0, 2, -4, 3, -1, 5, -2, -1, 3, 2, -4, 2, -2, 5, -2, 3, -5, -5, -6, 0, -7, 3, 2, -6}, {4, -4, 1, -1, -2, -4, -2, 5, 7, 2, 6, -5, 2, -3, 8, 6, -9, 5, 3, 4, -2, -3, 2, 3, -5, -2, 2, -3, 0, 6, -1, -8, 2, -2, 1, -3, 9, 9, 2, -3, -1, 0, 1, 8, 7, -7, 0, 2, 0, 0, 6, -4, -3, -5, 5, 0, 3, -11, 1, 3, 1, -1, 0, 1, -6, 3, -2, -6, 3, -2, 3, 1, -4, -8}, {2, 3, -2, -1, -4, 5, -6, 8, 0, -2, 0, -6, 2, 1, 2, -5, -1, 5, 0, 1, -2, 0, -2, 3, -5, 1, 2, 0, -2, 9, 3, 3, 5, 1, 4, -4, 5, 2, 6, 1, -6, 1, 1, -3, 0, 4, -1, 3, 1, -3, 4, -4, -2, -3, -4, 9, -3, -4, 9, -6, -1, -8, -2, 0, -7, 3, -2, -6, -3, -7, 0, 0, 1, -3}, {-7, 1, -3, 1, 3, 5, -2, 5, -6, -2, 7, -1, -5, -3, -4, -5, -1, -2, 2, -2, -4, -3, 1, 2, -2, -3, 1, 6, 14, 1, 4, 1, 1, -5, -6, -7, 2, 5, 3, 0, 4, -2, 2, -10, 1, 3, 6, -3, 1, -1, 1, 2, -2, -7, 0, 2, 0, 4, -1, -4, 0, -2, 0, 3, 6, 4, 2, -3, 0, 0, -8, -7, -6, 0}, {-6, -9, 4, 2, 0, -3, 6, -2, -3, 5, 1, 3, 1, 1, 9, -6, -1, 1, 1, 1, -2, 2, -4, -4, -7, -2, 3, 1, -3, -5, -3, -2, 2, -3, 5, -3, -4, -5, -1, -2, -2, -1, -4, -2, 3, -2, -8, 4, -3, -1, 1, -1, -3, 6, 3, -3, 5, 1, -6, -5, 8, -10, -3, -4, -2, -1, 1, 3, 1, 2, 8, -1, 4, -5}, {-3, 0, -5, 1, -3, 4, 4, -3, 5, 6, -8, -2, -2, 4, 12, -3, -10, 1, -2, -3, 0, 1, -6, -2, 1, 3, -3, -1, 2, -6, -1, 0, 0, 0, 7, 3, 5, 4, 7, 0, -4, -2, 0, -2, 2, -3, 1, -2, -6, 1, 2, -4, -2, -2, -4, -2, -2, -9, -2, 2, -1, 0, 4, -4, -3, 0, -5, 1, 1, -2, 4, -2, 4, 5}, {-1, 1, 4, 3, 2, 1, -5, 0, 0, -3, -1, 7, -4, -2, -4, 2, 2, 3, -1, -7, 5, -5, 1, 2, 3, 3, 2, -3, 5, 6, 2, -2, 3, 0, -1, -4, 7, 0, 0, 9, 4, 2, -1, -2, 5, -2, 6, -7, 3, -1, -1, -3, 5, 3, 0, 2, -1, -5, 5, 1, -10, -2, -4, -4, -5, 7, 1, 2, -3, -2, -2, -3, 6, 5}, {-2, 9, -6, 0, 5, 3, -1, -1, -6, 1, 5, -3, -9, 3, 2, -5, 0, -4, 5, -1, 3, -4, 1, -1, 1, 0, -1, -3, 5, 3, -4, 5, -5, 10, 2, 2, -6, -6, 1, -11, 3, 5, 4, -1, -1, -1, 2, -5, -2, 5, -7, -5, 1, 0, -1, 8, 2, 0, 5, -2, 5, 1, 0, -3, -1, 0, -3, 4, -5, 3, -1, -1, 3, 0}, {3, -7, 1, 2, 2, 2, -5, 8, 1, 0, -2, -7, -1, -1, 7, 4, 11, -4, 3, 6, 1, 4, -6, -4, 2, 1, -5, -2, -9, 6, 5, -2, -3, -9, -4, 4, -3, 5, 7, 3, -4, -4, 4, -5, -2, 2, -7, 4, -1, -2, 2, -3, 2, -4, 2, -1, 2, 0, -6, -2, -4, 1, -10, 6, 7, -6, -2, 7, 3, -4, 2, -2, 5, -1}, {1, 1, -4, -4, 1, 2, -7, -1, 1, 1, -2, -1, 6, 0, 1, 3, -3, 8, -1, 4, -2, 0, 2, -4, -9, 2, 4, -1, 0, 0, 1, 3, 3, 1, -6, -5, 1, 0, 2, -1, -3, -6, 1, -2, -1, 0, 2, 0, 1, 2, 7, 3, 4, 8, -1, -7, -4, -9, 1, 2, 5, -1, 3, 0, 2, -3, 3, -5, 0, 1, 7, -1, -4, 2}, {-1, 3, -2, -2, -3, 0, 1, -4, 2, -5, 1, -2, -5, -6, -9, 0, 2, 2, 0, 2, 2, 0, -2, 2, 0, -5, 0, -3, -1, -4, -2, 3, 1, -2, -3, 2, -9, 4, 6, -3, -9, 9, -9, -12, 1, 9, 2, 3, 2, -1, -3, 3, -6, -2, -5, -6, -3, -1, -2, 1, -1, 1, 2, -4, 5, -5, 3, -5, 0, 1, -1, 7, 1, -11}, {1, -8, -2, -5, 4, -2, -6, -1, 1, -1, 0, 2, -2, -3, 5, 2, 3, -3, -3, 8, 4, 1, 7, -3, 4, 8, 5, -6, -3, 1, 2, 0, 5, 5, 2, 4, -1, 2, 1, 4, 2, 8, -5, -5, 5, 0, -1, 1, -2, 6, 4, 3, 6, 5, -3, 1, -6, -2, -2, 5, -9, -2, 0, 0, 2, -3, -7, 3, -3, 1, 3, 2, -1, 3}, {-8, 3, -6, 8, -2, -4, 17, 1, -2, -2, -1, -7, 6, 0, -7, 1, 0, -3, -3, 3, 0, 5, 1, 0, 4, -7, 3, 1, -5, 6, -3, -2, -3, -5, 6, -5, -3, -2, 8, -3, -1, 2, 8, 1, -2, 7, 0, 3, 0, 0, -6, 4, 1, 2, -2, -2, -2, -1, 2, 1, 0, -1, 0, -3, -1, -3, 3, 5, -2, 4, -7, -3, -5, -2}, {4, 1, -6, 0, 3, -3, 5, 3, 0, -2, 7, -1, -2, 4, 1, -1, -4, 1, 0, 1, 7, 5, -1, -2, -2, -5, -5, 1, -2, -5, -5, -1, -3, -2, 2, 6, 3, -7, -1, 0, 1, 8, -2, 2, -1, -5, 2, -3, -1, -1, -2, -1, 5, 5, -1, -4, -4, 4, -3, 1, 1, 4, -2, 2, 5, -2, 4, 5, -2, 0, 4, 0, 2, -1}, {-7, -11, 9, 2, 7, -2, 1, 2, -3, 6, 0, 2, 3, -6, 0, -1, -1, -1, -3, 5, -2, -1, 0, -6, 1, 5, 2, 4, 11, 6, -3, -8, -2, -4, 3, 0, 4, 3, -6, 4, 5, -4, -1, -2, 3, -4, -2, 3, -5, 0, -2, -2, -1, -8, 6, -4, -2, -4, -2, -2, -6, -1, -2, 5, 2, 2, 0, 4, -2, 2, -3, 0, -1, -4}, {-4, -4, -2, -3, 1, 0, -4, 3, -4, -2, 2, 5, 0, 0, 3, -2, 2, 3, -1, 1, 3, 1, -3, 4, 2, 3, 0, 2, -5, 5, -3, -4, 2, 0, 3, -5, 4, 5, -12, 7, -2, 2, 0, -1, 4, 1, -4, 2, 3, 4, 2, 3, 2, 7, -5, 1, -8, 3, -2, 0, -4, 1, 3, -1, 2, 7, -1, 2, 0, -4, -2, -1, 1, -1}, {-8, 2, 2, 4, 4, 4, 6, -5, -6, 3, -7, -1, 3, -7, 3, 4, 1, 1, -3, 1, -1, -3, 0, -5, 6, -6, 1, 4, -1, 6, 7, -3, -5, -2, 4, -1, -9, 6, 3, 3, -3, 0, 3, -2, -4, 5, 8, 3, -4, -4, -2, -3, 4, 0, 3, 2, 4, 0, 3, -3, -1, 3, 2, -1, -5, 3, -5, -2, -2, 5, 3, -3, 5, -2}, {5, 10, 1, -7, 8, 1, -6, -4, -6, 2, -8, -5, -4, 1, 5, 5, 6, 4, 5, 0, 1, -1, 1, -4, 0, -2, 2, -2, 3, -2, -1, 4, -6, 5, -7, 0, 3, 1, -10, 5, -3, 3, 0, -8, -1, -1, 2, -3, -6, 4, 0, -4, -2, -1, -3, -1, 1, -7, 0, 1, 6, -4, 7, -6, -1, -5, -2, 3, 6, -4, 4, 1, -2, 6}, {5, -2, 1, 3, 1, 2, -5, 3, 0, 2, -1, -4, 4, 2, -2, 1, 3, 1, 1, 3, 5, 5, -2, -1, 2, -5, -3, 2, -10, -4, 4, -1, 1, -4, -2, 0, -1, -5, 1, -4, 5, 0, 2, -3, -12, -4, -9, -4, 2, 2, 4, -3, 1, 6, -1, -5, -6, 3, 2, 0, -2, 5, 4, 3, -4, -3, 1, -4, 2, 4, -3, 4, -4, 5}, {-2, 8, 5, 3, 4, -5, 8, -1, -2, 0, 4, -4, 1, -8, 0, 4, -3, 5, 4, 0, -1, -5, 0, 13, 0, -5, -1, 4, 9, 3, 5, -6, -2, -1, 3, 6, -4, 2, -6, 0, 9, 2, -2, -9, -3, 4, -5, 2, 3, -1, -2, -2, -11, -5, 1, 3, -1, -11, 6, -2, -4, -1, 4, 2, -5, 2, 5, 0, 6, -2, 7, 5, -4, -2}, {-6, -1, 8, 5, -1, -1, -2, 2, 2, -3, 1, -3, 5, -6, -4, 6, 6, 0, -4, -1, -2, 1, -3, 1, 3, 4, -1, -3, -3, 6, 6, -7, 2, 1, 2, -4, 1, -1, 0, 9, 0, -6, -8, -1, -7, -1, 1, -1, -4, 0, -2, 3, 1, 0, 2, 7, 2, 4, -2, -6, -1, -5, -8, -5, 2, 1, 1, 3, 2, 5, -1, 2, -4, 0}, {2, -6, -1, 7, 2, 3, -1, -2, -5, 0, 3, 0, -3, 1, -5, -2, 4, -1, -1, 5, 3, -2, 4, -3, -3, -3, 5, 0, 3, 0, -6, -4, 2, 1, 0, -1, 4, 7, 0, 1, -5, 4, 5, -4, 1, 0, -2, 2, 0, 1, 5, 0, -2, 5, 1, -1, -5, -4, -3, 0, 0, -6, 4, 0, 4, -5, 6, 6, -6, 4, -5, -1, -1, 2}, {1, 9, -3, 7, 0, -2, -4, -1, 0, -4, 5, -4, -2, 0, -5, 9, 5, -1, -1, -5, 2, -6, 5, 5, 4, -5, -5, 3, -2, 1, 1, 1, 0, -3, -1, -4, -7, 0, 5, -2, 5, -3, 4, 1, 1, 1, -2, -1, 4, -9, 0, -6, -7, 0, -1, 6, 2, 2, 2, 0, -5, 0, 1, -1, -5, -2, 14, 5, -1, 8, -3, 0, -5, -2}, {9, 6, 2, -2, -1, 0, -4, -4, 0, -3, -1, -2, -1, 2, 0, 7, 3, -1, -3, 5, 1, 4, -6, -1, -6, 1, -2, -2, -5, 4, 2, 3, -7, 4, 5, -2, -1, -2, 1, 0, 0, 3, -8, -6, -7, -7, -4, -4, -12, -1, 3, 5, -2, -1, -1, 4, -3, 2, 3, -2, 3, -4, 0, 3, -6, -2, -6, 2, 5, -3, -6, -5, -6, 0}, {1, 1, -3, 4, -4, -3, 1, -6, -1, 0, 4, 2, -2, 1, 0, -1, 4, 7, -2, -3, 0, 0, -3, 4, 6, -1, -1, -7, -2, -3, -2, 1, -6, 6, -3, 3, -1, -4, -2, 1, 9, 1, 0, 4, 1, 2, -2, -2, -7, 0, 3, -1, -5, 2, -5, 6, -1, -1, 4, -3, -8, -4, 4, -1, -6, -7, -6, -4, -2, 4, -6, 5, -2, -7}, {-7, 1, -6, 3, 3, -12, 4, 0, 5, 3, 2, 0, -2, -4, 2, 3, -3, -4, 2, 3, -7, -5, 5, -3, 5, 0, -5, 3, 8, -1, -5, -1, 3, 1, 5, 3, 1, 2, 2, -5, 4, -8, 3, 4, -4, -1, 12, -2, 0, 7, -5, 1, -2, -2, 0, 10, -4, -4, 0, 3, -6, 5, 3, -4, 2, 0, 0, 1, -5, 4, -1, 0, -9, 0}, {-1, -7, 4, 1, -1, -1, 2, 2, 1, 7, -2, 2, 4, -3, 1, -5, -3, 1, 7, -4, 0, -10, 1, -9, -3, -1, -1, 3, 6, -9, 2, -4, -2, -8, -6, 4, 3, -8, -6, -1, -3, -6, -3, 0, -1, -4, 3, -4, 2, 6, 3, -1, 0, -5, 2, -4, 2, 4, -4, -6, 3, -2, -1, 0, 1, 0, -2, 6, 1, 2, 3, -1, 3, 6}, {5, 6, 1, 1, 5, -3, 0, 5, -2, 3, -6, 2, 0, -7, -4, -3, 1, -1, -1, 0, 3, -4, 7, 10, -6, -3, 5, -1, -1, -1, -3, 6, -5, -2, -5, 2, -5, 1, -2, 7, -3, 3, 1, 0, -3, 0, 1, 2, 8, -3, -2, -5, -2, -1, -1, -2, 2, -3, 5, 3, -7, 1, -4, 0, -3, -2, -2, 1, -8, 4, 3, 6, -4, -5}, {0, 7, -3, -3, -4, -5, 7, 6, -4, 2, -2, -2, -2, -4, -8, -6, 2, -5, 5, -6, -6, -4, 4, 6, 4, 1, -2, 1, -2, 1, 10, 8, 1, -1, -2, 2, -4, -4, -2, 1, 7, -1, -1, 6, 1, 3, 0, 1, 3, 0, 0, -2, 1, 0, -1, 9, 1, -2, -3, 4, -9, 4, 5, -5, -4, 4, 5, 6, 1, 0, 5, 1, -3, 8}, {4, 4, 1, -7, -4, 6, 1, 2, -2, -2, -4, 1, 8, 3, -5, 1, -6, -1, -1, -5, -3, -5, -3, -2, -1, 6, -3, 7, -3, -3, -1, 5, 3, 3, -5, -1, 3, 0, 0, 2, -5, 9, -6, -4, -2, -2, -5, 2, -3, 3, 0, 1, 3, -2, -7, -1, -6, 0, 8, 3, -1, 3, 2, -3, -5, -1, -1, 4, -1, -9, -6, 5, 1, 7}, {1, -8, -5, 4, -5, 1, 0, 10, -2, 1, 2, -2, -5, 2, -6, -2, 5, -5, 7, -1, -2, 2, -2, -2, -3, -7, -1, -7, -2, 4, 2, 1, -1, -2, -2, 1, -1, -3, 5, -2, 3, -5, 5, -2, -3, 4, -2, -4, 0, 0, 3, 1, 1, 0, 5, 4, 3, -2, 3, -4, 4, -3, -2, 4, -2, -7, 3, 6, 4, -2, 4, -3, 3, 1}, {6, -9, 3, 3, 3, 1, 1, 2, -4, 1, 5, -7, 4, -2, 3, -2, -4, 4, 4, -2, -4, -5, -1, 0, 0, -3, -4, 1, -2, 0, 5, -2, -3, 7, -5, -2, 1, 8, 0, 1, -2, 0, -3, -1, -2, 2, -9, 2, 0, 3, 3, -5, -5, -1, 3, -3, 3, -9, 5, 2, 4, -7, 8, -5, -2, 0, 1, 4, 2, -2, 2, -2, 0, -1}, {-2, 2, -1, 2, 6, 3, -4, -3, -6, 6, 2, 2, 7, -3, 3, -5, -3, -2, -1, -1, 5, -2, 5, 3, -4, 1, -1, -1, -5, 11, 7, 1, 4, -1, -2, -9, -7, -4, -6, 3, 5, -1, 0, 3, 6, 3, -3, 0, 2, 2, 6, -6, 6, -1, 0, 1, 8, -5, 2, 7, 1, 1, 3, -1, -9, -1, 6, 4, 4, 4, 7, 2, 5, 3}, {-4, 2, 0, 6, -3, 1, 5, 3, 0, 3, -5, -2, 4, 0, -2, -3, -6, -3, -6, -1, 6, -4, -2, 6, 3, -2, 5, 1, -3, 0, 2, -3, -8, -3, -1, 1, -3, -11, 0, -1, 3, -5, 4, 2, -6, -3, 3, 0, 1, 0, -4, 5, -2, 5, 5, -5, 5, -2, -10, 1, -3, 3, 0, -6, 6, 0, 0, 2, 0, 6, -2, -1, 6, 4}, {4, 2, -2, -4, 5, -3, -1, -3, 2, 4, -9, -5, -6, -5, 4, 1, 0, -2, 1, 1, 3, -2, -9, -4, 0, 3, -1, -7, -9, -4, -6, 0, 1, -1, 4, 8, 0, -4, 3, 3, -7, 1, -7, 2, 0, -6, 4, -2, -3, 2, 1, -2, -1, -2, 5, 5, 4, 0, -2, 5, -2, 7, -2, 3, 4, -2, -11, 2, 0, 3, 1, 7, -2, 4}, {0, 0, 1, 8, 5, -3, 2, 2, 0, -1, 1, 7, -1, 2, 2, -1, -5, 0, 6, 0, -1, 1, -4, 6, 0, 1, 3, 3, 5, 1, 1, -4, -1, -8, -2, -4, 3, 4, -1, -1, 1, -7, 1, 1, -5, -2, -4, 0, 5, 7, 5, 2, -3, -4, 7, -1, 4, -1, -1, 2, 1, 2, -5, 1, 3, 5, 2, -7, -7, -1, 2, 3, -2, 3}, {0, -1, 12, -3, 5, -1, 2, -5, 1, -4, 2, 9, 3, -6, 2, 1, -3, -2, 2, 0, -3, -8, -3, -6, 5, 5, 4, 4, 6, 6, -2, -8, 2, -1, -4, -3, 0, 2, 0, 4, 0, 6, -8, 0, 7, -2, -1, -4, -4, 3, -8, 2, 1, 2, 8, 3, 2, -5, 3, 3, 4, -7, -4, -2, 1, 8, 1, 1, 3, -4, 3, -6, -1, 1}, {-5, 5, -4, -2, -8, -5, 1, 2, -1, 1, 3, -2, -8, 6, 1, 2, -4, 4, 0, -3, 1, 2, 0, 5, 4, -1, 0, -9, -3, 4, -4, -3, 4, -5, 0, 4, 4, -3, 3, -5, 7, -4, -3, 3, 5, 1, 8, 1, 0, 1, 3, 0, 7, -1, 1, 4, 9, 0, 4, -1, 7, 4, 0, -3, 8, -2, 1, -3, 1, 1, 3, 1, -1, -4}, {2, 5, -4, 5, -1, 2, -5, 4, -6, -1, -8, -1, -4, -7, 3, -3, 7, -5, 2, 4, 6, -6, 2, 3, 0, 0, 7, -2, -4, 3, -1, 2, 1, -1, -5, -6, -2, -4, 3, -2, -9, -3, -2, -3, -2, 5, -1, 1, 5, 2, -3, 7, -3, -2, 0, 3, 4, 0, -2, 6, -4, -4, 0, 6, -6, 0, 4, 6, -5, 1, 2, 2, 1, 0}, {4, 2, -1, 2, -1, 5, -10, -2, -6, -1, -1, 3, -8, 3, -1, -5, 3, -6, 1, 4, -1, 0, -5, 7, -1, 4, 4, 2, -2, -2, 0, -6, -1, -4, -5, -5, 6, 1, 4, -4, -6, 0, 3, -6, 5, -2, -7, 4, 3, 1, -1, 4, -3, 7, -1, -3, -3, 0, -1, 1, 1, -4, -2, -2, -1, 3, -1, -2, -5, -6, -3, -2, -2, 3}, {7, 5, -1, 0, -1, 2, -2, -3, 2, -3, 3, -4, 0, 6, -7, 9, 0, 1, -3, -4, -3, -2, 4, -7, 0, -6, -3, 1, -3, -6, 3, 8, -4, 4, -10, 6, -3, 5, -4, 1, 2, 7, -5, 5, -2, 0, 3, 0, -5, -3, 0, -1, 2, 0, 0, -6, 1, -1, 2, 0, -2, 6, 0, -4, -1, -5, -3, -2, 6, 4, 3, 6, 8, 4}, {1, 3, -8, 6, -3, -1, 0, -1, -1, 0, 3, -3, 4, 10, 3, -3, 1, 0, -6, 5, -3, 0, 7, -4, 5, 3, -7, 4, 2, -2, 4, -1, -3, 7, -7, -10, -7, 1, 2, 3, 0, 5, 2, 3, -6, 0, -2, 1, -6, 0, -2, -2, -4, -2, -8, 1, 3, 0, 10, 8, 0, 2, 0, 5, -2, 4, 8, 2, 1, 9, -6, 6, -8, -1}};

void random_base(uint8_t *r);
void printf_base(uint8_t key[], uint8_t r, char *c);
void base_to_key(int8_t base[], uint8_t key[]);
void abs_key(uint8_t key[], uint8_t abs[]);

// (each entry corresponds to the number of degree-(l_i) to be required in the action: this the one given in Onuki et al. work)
static int8_t B[] =	{
 2,  3, 3, 3, 3, 3,  3,  3, 
 3,  3, 3, 4, 4, 4,  4,  4, 
 4,  4, 4, 4, 4, 4,  4,  4, 
 4,  4, 5, 5, 5, 5,  5,  5, 
 5,  6, 6, 6, 6, 6,  7,  7, 
 7,  7, 7, 7, 7, 7,  7,  7, 
 7,  7, 8, 9, 9, 9, 10, 10,
10, 10, 9, 8, 8, 8,  7,  7,
 7,  7, 7, 6, 5, 1,  2,  2,
 2,  2
};

// (NUMBER_OF_BATCHES) different subsets (i.e., batches)
static uint8_t BATCH_0[] = { 0, 3, 6, 9, 12, 15, 18, 21, 24, 27, 30, 33, 36, 39, 42, 45, 48, 51, 54, 57, 60, 63, 66, 69, 72 };
static uint8_t BATCH_1[] = { 1, 4, 7, 10, 13, 16, 19, 22, 25, 28, 31, 34, 37, 40, 43, 46, 49, 52, 55, 58, 61, 64, 67, 70, 73 };
static uint8_t BATCH_2[] = { 2, 5, 8, 11, 14, 17, 20, 23, 26, 29, 32, 35, 38, 41, 44, 47, 50, 53, 56, 59, 62, 65, 68, 71 };

static uint8_t SIZE_OF_EACH_BATCH[NUMBER_OF_BATCHES] = {25, 25, 24};
static uint8_t *BATCHES[NUMBER_OF_BATCHES] = { BATCH_0, BATCH_1, BATCH_2 };

static uint8_t LAST_ISOGENY[NUMBER_OF_BATCHES] = { 72, 73, 71 };
static uint16_t NUMBER_OF_ISOGENIES = 404;

// The complement of each batch
static uint8_t SIZE_OF_EACH_COMPLEMENT_BATCH[NUMBER_OF_BATCHES] = { 49, 49, 50 };
static uint8_t COMPLEMENT_OF_EACH_BATCH[NUMBER_OF_BATCHES][N] = {
{  1,  2,  
   4,  5, 
   7,  8,
  10, 11,
  13, 14,
  16, 17,
  19, 20,
  22, 23,
  25, 26,
  28, 29,
  31, 32,
  34, 35,
  37, 38,
  40, 41,
  43, 44,
  46, 47,
  49, 50,
  52, 53,
  55, 56,
  58, 59,
  61, 62,
  64, 65,
  67, 68,
  70, 71,
  73, 
  74, 74,
  74, 74,
  74, 74,
  74, 74,
  74, 74,
  74, 74,
  74, 74,
  74, 74,
  74, 74,
  74, 74,
  74, 74,
  74, 74,
  74  
},
{  0,  2,
   3,  5,
   6,  8,
   9, 11,
  12, 14,
  15, 17,
  18, 20,
  21, 23,
  24, 26,
  27, 29,
  30, 32,
  33, 35,
  36, 38,
  39, 41,
  42, 44,
  45, 47,
  48, 50,
  51, 53,
  54, 56,
  57, 59,
  60, 62,
  63, 65,
  66, 68,
  69, 71,
  72, 
  74, 74,
  74, 74,
  74, 74,
  74, 74,
  74, 74,
  74, 74,
  74, 74,
  74, 74,
  74, 74,
  74, 74,
  74, 74,
  74, 74,
  74
},
{  0,  1,
   3,  4,
   6,  7,
   9, 10,
  12, 13,
  15, 16,
  18, 19,
  21, 22,
  24, 25,
  27, 28,
  30, 31,
  33, 34,
  36, 37,
  39, 40,
  42, 43,
  45, 46,
  48, 49,
  51, 52,
  54, 55,
  57, 58,
  60, 61,
  63, 64,
  66, 67,
  69, 70,
  72, 73, 
  74, 74,
  74, 74,
  74, 74,
  74, 74,
  74, 74,
  74, 74,
  74, 74,
  74, 74,
  74, 74,
  74, 74,
  74, 74,
  74, 74
}
};
#endif