/**
 * \file Needleman-Wunsch-recmemo.h
 * \brief recursive implementation with memoization of Needleman-Wunsch global alignment algorithm that computes the distance between two genetic sequences 
 * \version 0.1
 * \date 03/10/2022 
 * \author Jean-Louis Roch (Ensimag, Grenoble-INP - University Grenoble-Alpes) jean-louis.roch@grenoble-inp.fr
 */

#include <stdlib.h> /* for size_t */

/*
 * Costs for operations on canonical bases
 * Three  operations: insertion and sustitution of one base by an another 
 * Note= substitution of an unknown base N by another one (known or unknown) as the same cost than substitution between 2 different known bases
 */
/** \def SUBSTITUTION_COST
 *  \brief Cost of substitution of one canonical base by another
 */
#define SUBSTITUTION_COST 1

/** \def SUBSTITUTION_UNKNOWN_COST
 *  \brief Cost of substitution of an unknown base (N) by another one (canonical or unknown)
 */
#define SUBSTITUTION_UNKNOWN_COST 1 /* Cost for sustitition of an Unknown bas N by another on -known or unkown- */

/** \def INSERTION_COST
 *  \brief Cost of insertion of a canonical base 
 */
#define INSERTION_COST 2

/********************************************************************************
 * Recursive implementation of NeedlemanWunsch with memoization
 */
/**
 * \fn long EditDistance_NW_Rec(char* A, size_t lengthA, char* B, size_t lengthB);
 * \brief computes the edit distance between A[0 .. lengthA-1] and B[0 .. lengthB-1]
 * \param A  : array of char represneting a genetic sequence A 
 * \param lengthA :  number of elements in A 
 * \param B  : array of char represneting a genetic sequence B
 * \param lengthB :  number of elements in B 
 * \return :  edit distance between A and B }
 *
 * editDistance_RecMemo is a memoized recursive immplementatioin of Needleman-Wunsch algorithm.
 * It allocates the data structure for memoization table and calls the internal recursive function _editDistance_memo
 * that fills in the memoization table.
 * 
 * If lengthA < lengthB, the sequences A and B are swapped.
 *
 */
long EditDistance_NW_Rec(char *A, size_t lengthA, char *B, size_t lengthB);

/********************************************************************************
 * Iterative implementation of NeedlemanWunsch with memoization
 */
/**
 * \fn long EditDistance_NW_Rec(char* A, size_t lengthA, char* B, size_t lengthB);
 * \brief computes the edit distance between A[0 .. lengthA-1] and B[0 .. lengthB-1]
 * \param A  : array of char represneting a genetic sequence A 
 * \param lengthA :  number of elements in A 
 * \param B  : array of char represneting a genetic sequence B
 * \param lengthB :  number of elements in B 
 * \return :  edit distance between A and B }
 *
 * editDistance_iteratif : Using an array of (N+1) elements when traversed M times (with M the length of the largest sequence 
 * and N that of the smallest). After each traversal, Bellman's equation is applied with a storage of the prev_value"which is
 * the old value of the neighbor to the left of our iterator before the last traversal.
 * Our array is renewed after each scan, until the end and the value of tab[N] is thus the value sought.
 */
long EditDistance_NW_iteratif(char *A, size_t lengthA, char *B, size_t lengthB);

/**
 * \fn long EditDistance_cache_aware(char* A, size_t lengthA, char* B, size_t lengthB, int Z);
 * \brief computes the edit distance between A[0 .. lengthA-1] and B[0 .. lengthB-1]
 * \param A  : array of char represneting a genetic sequence A 
 * \param lengthA :  number of elements in A 
 * \param B  : array of char represneting a genetic sequence B
 * \param lengthB :  number of elements in B 
 * \param Z : cache size
 * \return :  edit distance between A and B }
 *
 * editDistance_cache_aware :the same principle as the iterative version, except that here we take the size
 * of the cache Z into consideration. The table is divides into mini tables of nbr_case elements.
 * Hence the usefulness of the col array of M+1 elements which makes it possible to store the ancient values needed for
 * Bellman's esquation.
 */
long EditDistance_NW_cache_aware(char *A, size_t lengthA, char *B, size_t lengthB, int Z);

/**
 * \fn long EditDistance_cache_oblivious(char* A, size_t lengthA, char* B, size_t lengthB, int seuil);
 * \brief computes the edit distance between A[0 .. lengthA-1] and B[0 .. lengthB-1]
 * \param A  : array of char represneting a genetic sequence A 
 * \param lengthA :  number of elements in A 
 * \param B  : array of char represneting a genetic sequence B
 * \param lengthB :  number of elements in B 
 * \param seuil : a threshold for the calculation
 * \return :  edit distance between A and B }
 *
 * editDistance_cache_oblivious :From the iterative version, we take a threshold for the number of boxes
 * that we can calculate, we do a recursion: that is to say that we divide our N in two, and then we call
 * the function for half of the sequence and then the second half. and still using our col array of M+1
 * elements which allows us to store the values. this recursion is done using the fonction cache_oblivious_helper
 */
long EditDistance_NW_cache_oblivious(char *A, size_t lengthA, char *B, size_t lengthB, int seuil);

/**
 * \fn void cache_oblivious_helper(char *X, size_t M, char *Y, size_t N, long *col, int seuil, long debut_seq, long fin_seq);
 * \brief helps in findint the edit distance between A[0 .. lengthA-1] and B[0 .. lengthB-1]
 * \param X  : array of char represneting a genetic sequence X 
 * \param M :  number of elements in X 
 * \param Y  : array of char represneting a genetic sequence Y
 * \param N :  number of elements in Y
 * \param seuil : a threshold for the calculation
 * }
 *
 * cache_oblivious_helper : it helps in the recursion needed for the cache oblivious version
 */
static void cache_oblivious_helper(char *X, size_t M, char *Y, size_t N, long *col, int seuil, long debut_seq, long fin_seq);
