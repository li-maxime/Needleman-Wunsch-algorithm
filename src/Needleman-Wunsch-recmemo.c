/**
 * \file Needleman-Wunsch-recmemo.c
 * \brief recursive implementation with memoization of Needleman-Wunsch global alignment algorithm that computes the distance between two genetic sequences 
 * \version 0.1
 * \date 03/10/2022 
 * \author Jean-Louis Roch (Ensimag, Grenoble-INP - University Grenoble-Alpes) jean-louis.roch@grenoble-inp.fr
 *
 * Documentation: see Needleman-Wunsch-recmemo.h
 * Costs of basic base opertaions (SUBSTITUTION_COST, SUBSTITUTION_UNKNOWN_COST, INSERTION_COST) are
 * defined in Needleman-Wunsch-recmemo.h
 */

#include "Needleman-Wunsch-recmemo.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h> /* for strchr */
// #include <ctype.h> /* for toupper */

#include "characters_to_base.h" /* mapping from char to base */

/*****************************************************************************/

/* Context of the memoization : passed to all recursive calls */
/** \def NOT_YET_COMPUTED
 * \brief default value for memoization of minimal distance (defined as an impossible value for a distance, -1).
 */
#define NOT_YET_COMPUTED -1L

/** \struct NW_MemoContext
 * \brief data for memoization of recursive Needleman-Wunsch algorithm 
*/
struct NW_MemoContext
{
	char *X;	 /*!< the longest genetic sequences */
	char *Y;	 /*!< the shortest genetic sequences */
	size_t M;	 /*!< length of X */
	size_t N;	 /*!< length of Y,  N <= M */
	long **memo; /*!< memoization table to store memo[0..M][0..N] (including stopping conditions phi(M,j) and phi(i,N) */
};

/*
 *  static long EditDistance_NW_RecMemo(struct NW_MemoContext *c, size_t i, size_t j) 
 * \brief  EditDistance_NW_RecMemo :  Private (static)  recursive function with memoization \
 * direct implementation of Needleman-Wursch extended to manage FASTA sequences (cf TP description)
 * \param c : data passed for recursive calls that includes the memoization array 
 * \param i : starting position of the left sequence :  c->X[ i .. c->M ] 
 * \param j : starting position of the right sequence :  c->Y[ j .. c->N ] 
 */
static long EditDistance_NW_RecMemo(struct NW_MemoContext *c, size_t i, size_t j)
/* compute and returns phi(i,j) using data in c -allocated and initialized by EditDistance_NW_Rec */
{
	if (c->memo[i][j] == NOT_YET_COMPUTED)
	{
		long res;
		char Xi = c->X[i];
		char Yj = c->Y[j];
		if (i == c->M) /* Reach end of X */
		{
			if (j == c->N)
				res = 0; /* Reach end of Y too */
			else
				res = (isBase(Yj) ? INSERTION_COST : 0) + EditDistance_NW_RecMemo(c, i, j + 1);
		}
		else if (j == c->N) /* Reach end of Y but not end of X */
		{
			res = (isBase(Xi) ? INSERTION_COST : 0) + EditDistance_NW_RecMemo(c, i + 1, j);
		}
		else if (!isBase(Xi)) /* skip ccharacter in Xi that is not a base */
		{
			ManageBaseError(Xi);
			res = EditDistance_NW_RecMemo(c, i + 1, j);
		}
		else if (!isBase(Yj)) /* skip ccharacter in Yj that is not a base */
		{
			ManageBaseError(Yj);
			res = EditDistance_NW_RecMemo(c, i, j + 1);
		}
		else
		{			   /* Note that stopping conditions (i==M) and (j==N) are already stored in c->memo (cf EditDistance_NW_Rec) */
			long min = /* initialization  with cas 1*/
				(isUnknownBase(Xi) ? SUBSTITUTION_UNKNOWN_COST
								   : (isSameBase(Xi, Yj) ? 0 : SUBSTITUTION_COST)) +
				EditDistance_NW_RecMemo(c, i + 1, j + 1);
			{
				long cas2 = INSERTION_COST + EditDistance_NW_RecMemo(c, i + 1, j);
				if (cas2 < min)
					min = cas2;
			}
			{
				long cas3 = INSERTION_COST + EditDistance_NW_RecMemo(c, i, j + 1);
				if (cas3 < min)
					min = cas3;
			}
			res = min;
		}
		c->memo[i][j] = res;
	}
	return c->memo[i][j];
}

/* EditDistance_NW_Rec :  is the main function to call, cf .h for specification 
 * It allocates and initailizes data (NW_MemoContext) for memoization and call the 
 * recursivefunction EditDistance_NW_RecMemo 
 * See .h file for documentation
 */
long EditDistance_NW_Rec(char *A, size_t lengthA, char *B, size_t lengthB)
{
	_init_base_match();
	struct NW_MemoContext ctx;
	if (lengthA >= lengthB) /* X is the longest sequence, Y the shortest */
	{
		ctx.X = A;
		ctx.M = lengthA;
		ctx.Y = B;
		ctx.N = lengthB;
	}
	else
	{
		ctx.X = B;
		ctx.M = lengthB;
		ctx.Y = A;
		ctx.N = lengthA;
	}
	size_t M = ctx.M;
	size_t N = ctx.N;
	{ /* Allocation and initialization of ctx.memo to NOT_YET_COMPUTED*/
		/* Note: memo is of size (N+1)*(M+1) but is stored as (M+1) distinct arrays each with (N+1) continuous elements 
       * It would have been possible to allocate only one big array memezone of (M+1)*(N+1) elements 
       * and then memo as an array of (M+1) pointers, the memo[i] being the address of memzone[i*(N+1)].
       */
		ctx.memo = (long **)malloc((M + 1) * sizeof(long *));
		if (ctx.memo == NULL)
		{
			perror("EditDistance_NW_Rec: malloc of ctx_memo");
			exit(EXIT_FAILURE);
		}
		for (int i = 0; i <= M; ++i)
		{
			ctx.memo[i] = (long *)malloc((N + 1) * sizeof(long));
			if (ctx.memo[i] == NULL)
			{
				perror("EditDistance_NW_Rec: malloc of ctx_memo[i]");
				exit(EXIT_FAILURE);
			}
			for (int j = 0; j <= N; ++j)
				ctx.memo[i][j] = NOT_YET_COMPUTED;
		}
	}

	/* Compute phi(0,0) = ctx.memo[0][0] by calling the recursive function EditDistance_NW_RecMemo */
	long res = EditDistance_NW_RecMemo(&ctx, 0, 0);

	{ /* Deallocation of ctx.memo */
		for (int i = 0; i <= M; ++i)
			free(ctx.memo[i]);
		free(ctx.memo);
	}
	return res;
}

/* EditDistance_NW_iteratif : la version itérative de l'algorithme.
 * See .h file for documentation
 */
long EditDistance_NW_iteratif(char *A, size_t lengthA, char *B, size_t lengthB)
{
	_init_base_match();
	struct NW_MemoContext ctx;
	if (lengthA >= lengthB) /* X is the longest sequence, Y the shortest */
	{
		ctx.X = A;
		ctx.M = lengthA;
		ctx.Y = B;
		ctx.N = lengthB;
	}
	else
	{
		ctx.X = B;
		ctx.M = lengthB;
		ctx.Y = A;
		ctx.N = lengthA;
	}
	size_t M = ctx.M;
	size_t N = ctx.N;
	long tab[N + 1];
	long min;
	int delta;
	int prev_value;
	tab[0] = 0;
	//on initialise notre tableau
	for (int j = 1; j < N + 1; j++)
	{
		tab[j] = tab[j - 1] + (isBase(ctx.Y[N - j]) * 2);
	}
	//pour chaque parcours, on remet a jour le tableau
	for (int i = 1; i < M + 1; i++)
	{
		//prev_value nous permet de stocker l'ancienne valeur du voisin gauche de notre element
		prev_value = tab[0];
		tab[0] = tab[0] + 2 * isBase(ctx.X[M - i]);
		for (int j = 1; j < N + 1; j++)
		{
			if (!isBase(ctx.Y[N - j]))
			{
				prev_value = tab[j];
				tab[j] = tab[j - 1];
			}
			else if (!isBase(ctx.X[M - i]))
			{
				prev_value = tab[j];
			}
			else
			{
				min = (tab[j] < tab[j - 1]) ? tab[j] : tab[j - 1];
				min += INSERTION_COST;
				// delta est égale au cout du mismatch si c'est le cas
				delta = (int)!isSameBase(ctx.X[M - i], ctx.Y[N - j]);
				delta += prev_value;
				prev_value = tab[j];
				//on retrouve le min pour la nouvelle valeur
				tab[j] = (min < delta) ? min : delta;
			}
		}
	}
	return tab[N];
}


/* EditDistance_NW_cache_aware : la version cache aware de l'algorithme.
 * See .h file for documentation
 */
long EditDistance_NW_cache_aware(char *A, size_t lengthA, char *B, size_t lengthB, int Z)
{
	_init_base_match();
	struct NW_MemoContext ctx;
	if (lengthA >= lengthB) /* X is the longest sequence, Y the shortest */
	{
		ctx.X = A;
		ctx.M = lengthA;
		ctx.Y = B;
		ctx.N = lengthB;
	}
	else
	{
		ctx.X = B;
		ctx.M = lengthB;
		ctx.Y = A;
		ctx.N = lengthA;
	}
	size_t M = ctx.M;
	int nb_case = (int)Z / (5 * sizeof(long));
	long N = ctx.N;
	long tab[nb_case + 1];
	long col[M + 1];
	long min;
	int delta;
	long int bordure;
	int prev_value;
	tab[0] = 0;
	col[0] = 0;
	// on initialise notre col qui permet de stocker l'ancienne valeur du dernier element calculé pour chaque parcours
	for (int i = 1; i < M + 1; i++)
	{
		col[i] = col[i - 1] + (isBase(ctx.X[M - i]) * 2);
	}
	//le calcul se fait par nb_cases elements au fur et à mesure jusqu'à atteindre N
	while (N > 0)
	{
		bordure = (N < nb_case) ? N : nb_case;
		tab[0] = col[0];
		// on initialise le tableau 
		for (int j = 1; j < bordure + 1; j++)
		{
			tab[j] = tab[j - 1] + (isBase(ctx.Y[N - j]) * 2);
		}
		// mise à jour de la valeur de col[0] 
		col[0] = tab[bordure];

		//pour chaque parcours, on remet a jour le tableau
		for (int i = 1; i < M + 1; i++)
		{
			prev_value = tab[0];
			tab[0] = col[i];
			for (int j = 1; j < bordure + 1; j++)
			{
				if (!isBase(ctx.Y[N - j]))
				{
					prev_value = tab[j];
					tab[j] = tab[j - 1];
				}
				else if (!isBase(ctx.X[M - i]))
				{
					prev_value = tab[j];
				}
				else
				{
					min = (tab[j] < tab[j - 1]) ? tab[j] : tab[j - 1];
					min += INSERTION_COST;
					delta = (int)!isSameBase(ctx.X[M - i], ctx.Y[N - j]);
					delta += prev_value;
					prev_value = tab[j];
					tab[j] = (min < delta) ? min : delta;
				}
			}
			//mise à jour de la valeur de col[i] pour les prochains calculs 
			col[i] = tab[bordure];
		}
		N = N - nb_case;
	}
	//col[M] represente la dernière valeur calculé à la fin de la sequence Y 
	return col[M];
}

/* EditDistance_NW_cache_oblivious : la version cache oblivious de l'algorithme.
 * See .h file for documentation
 */
long EditDistance_NW_cache_oblivious(char *A, size_t lengthA, char *B, size_t lengthB, int seuil)
{
	_init_base_match();
	struct NW_MemoContext ctx;
	if (lengthA >= lengthB) // X is the longest sequence, Y the shortest
	{
		ctx.X = A;
		ctx.M = lengthA;
		ctx.Y = B;
		ctx.N = lengthB;
	}
	else
	{
		ctx.X = B;
		ctx.M = lengthB;
		ctx.Y = A;
		ctx.N = lengthA;
	}
	size_t M = ctx.M;
	long N = ctx.N;
	long col[M + 1];
	col[0] = 0;
	//on initialise le tableau colonne;
	for (int i = 1; i < M + 1; i++)
	{
		col[i] = col[i - 1] + (isBase(ctx.X[M - i]) * 2);
	}
	//on appelle la fonction qui fera les calculs en prenant en considération le seuil
	cache_oblivious_helper(ctx.X, ctx.M, ctx.Y, ctx.N, col, seuil, 0, N);
	return col[M];
}

/* EditDistance_NW_oblivious_helper : une fonction utilisé pour la version cache oblivious.
 * See .h file for documentation
 */
void cache_oblivious_helper(char *X, size_t M, char *Y, size_t N, long *col, int seuil, long debut_seq, long fin_seq)
{
	//la taille du sous tableau 
	long taille = fin_seq - debut_seq;
	if (taille > seuil)
	{
		/* notre taille supérieure au seuil 
		 * on divise notre taille en deux et on appelle par recursivité pour les deux sous tableaux
		 */
		long milieu = taille / 2;
		milieu += debut_seq;
		cache_oblivious_helper(X, M, Y, N, col, seuil, debut_seq, milieu);
		cache_oblivious_helper(X, M, Y, N, col, seuil, milieu, fin_seq);
	}
	/* notre taille est inférieure au seuil: 
	* la suite revient aux étapes de la version itérative pour le sous tableau
	*/
	else
	{
		long tab[taille + 1];
		long prev_value;
		tab[0] = col[0];
		long min;
		long delta;
		// on initialise notre sous-tableau
		for (int j = 1; j < taille + 1; j++)
		{
			tab[j] = tab[j - 1] + (isBase(Y[N - j - debut_seq]) * 2);
		}
		col[0] = tab[fin_seq];
		for (int i = 1; i < M + 1; i++)
		{
			prev_value = tab[0];
			tab[0] = col[i];
			for (int j = 1; j < taille + 1; j++)
			{
				if (!isBase(Y[N - j - debut_seq]))
				{
					prev_value = tab[j];
					tab[j] = tab[j - 1];
				}
				else if (!isBase(X[M - i]))
				{
					prev_value = tab[j];
				}
				else
				{
					min = (tab[j] < tab[j - 1]) ? tab[j] : tab[j - 1];
					min += INSERTION_COST;
					delta = (int)!isSameBase(X[M - i], Y[N - j - debut_seq]);
					delta += prev_value;
					prev_value = tab[j];
					tab[j] = (min < delta) ? min : delta;
				}
			}
			// mise à jour de col[i]
			col[i] = tab[taille];
		}
	}
}
