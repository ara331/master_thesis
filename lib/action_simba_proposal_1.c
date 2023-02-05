#include "edwards_curve.h"

void random_key(uint8_t key[])
{
	uint8_t i, tmp;
	int8_t exp, sgn;
	for(i = 0; i < N; i++)
	{
		// exp is randomly selected from |[ 0, 2B ]|
		randombytes(&tmp, 1);
		while ( issmaller((int32_t)B[i] << 1, (int32_t)tmp) == -1 )	// constant-time comparison
			randombytes(&tmp, 1);

		exp = (int8_t)tmp;
		// Mapping integers from |[ 0, 2B |] into |[ -B, B ]|
		exp = exp - (int8_t)B[i];
		sgn = exp >> 7;	// sign of exp

		// Next, to write  key[i] = e || ((1 + sgn)/2)
		cmov(&exp, -exp, sgn == -1);
		key[i] = (exp << 1) ^ (1 & (1 + sgn));
	};
};

void printf_key(uint8_t key[], char *c)
{
	int i;
	printf("%s := ", c);
	printf("{\t  %3d", (int)( (2*(key[0] & 0x1) - 1) * (key[0] >> 1) ));

	for(i = 1; i < N; i++)
	{
		printf(", %3d", (int)( (2*(key[i] & 0x1) - 1) * (key[i] >> 1) ) );
		if( (i % 18) == 17 )
			printf("\n\t\t");
	};

	printf("};\n");

};

/* ----------------------------------------------------------------------------------------------- *
   action_evaluation()
   inputs: a the secret key, the Edwards curve constants A[0]:=a, and A[1]:=(a - d);
   output: the isogenous Edwards curve constants C[0]:=a' and C[1]:=(a' - d') determined by the action 
           evaluated at the secret key and public curve A
   
    NOTE: As far as we've understood how simba works; this next code implements simba approach. The 
          action computed by the next code uses two torsion points T_{+} and T_{-}, which their affine
          y-coordinates (in the isomorphic Montgomery curve) belongs to Fp and Fp2\Fp, respectively.
 * ----------------------------------------------------------------------------------------------- */
// void action_evaluation(proj C, const uint8_t key[], const proj A)
void action_evaluation(proj C1, proj C2, const uint8_t key[], const proj A1, const proj A2)
{
	// --------------------------------------------------------------------------------------------------------
	// SIMBA parameters
	// Batches
	uint8_t batches[NUMBER_OF_BATCHES][SIZE_OF_EACH_BATCH[0]];
	uint8_t size_of_each_batch[NUMBER_OF_BATCHES];

	for(uint8_t i = 0; i < NUMBER_OF_BATCHES; i++)
		memcpy(batches[i], BATCHES[i], sizeof(uint8_t) * SIZE_OF_EACH_BATCH[i]);

	memcpy(size_of_each_batch, SIZE_OF_EACH_BATCH, sizeof(uint8_t) * NUMBER_OF_BATCHES);
	// Complement of each batch
	uint8_t complement_of_each_batch[NUMBER_OF_BATCHES][N];
	uint8_t size_of_each_complement_batch[NUMBER_OF_BATCHES];

	memcpy(complement_of_each_batch, COMPLEMENT_OF_EACH_BATCH, sizeof(uint8_t) * NUMBER_OF_BATCHES * N);
	memcpy(size_of_each_complement_batch, SIZE_OF_EACH_COMPLEMENT_BATCH, sizeof(uint8_t) * NUMBER_OF_BATCHES);
	// --------------------------------------------------------------------------------------------------------

	// --------------------------------------------------------------------------------------------------------
	// Copy of public and private data (the private key is modified each iteration)
	uint8_t tmp_e[N];
	memcpy(tmp_e, key, sizeof(uint8_t) * N);	// exponents

	proj current_A[2], current_T[4];
	point_copy(current_A[0], A1);			// initial Edwards curve constants a and (a -d)
    point_copy(current_A[1], A2);
	// --------------------------------------------------------------------------------------------------------

	// --------------------------------------------------------------------------------------------------------
	// Vairables required for running SIMBA
	int8_t ec = 0;
	uint16_t count = 0;
	proj G[4], K[(LARGE_L >> 1) + 1];		// Current kernel
	uint8_t finished[N];				// flag that determines if the maximum number of isogeny constructions has been reached
	memset(finished, 0, sizeof(uint8_t) * N);

	int8_t counter[N];				// This variable determines how many isogeny construcctions has been perfomed
	memset(counter, 0, sizeof(int8_t) * N);
	memcpy(counter, B, sizeof(int8_t) * N);		// At the beginning, we must perfomed b_i isogeny constructions for each l_i
	uint64_t isog_counter = 0;			// Total number of isogeny construction perfomed

	uint8_t last_isogeny[NUMBER_OF_BATCHES];
	//index for skipping point evaluations (the last one of each batch)
	memcpy(last_isogeny, LAST_ISOGENY, sizeof(uint8_t) * NUMBER_OF_BATCHES);
	uint32_t bc;
	// --------------------------------------------------------------------------------------------------------

	// --------------------------------------------------------------------------------------------------------
	// Main loop
	uint8_t m = 0, i, j;
	uint64_t number_of_batches = NUMBER_OF_BATCHES;

	while (isog_counter < NUMBER_OF_ISOGENIES)
	{
		m = (m + 1) % number_of_batches;
		
		if(count == MY*number_of_batches) {  	//merge the batches after my rounds
			m = 0;
			size_of_each_complement_batch[m] = 0;
			size_of_each_batch[m] = 0;
			number_of_batches = 1;

			for(i = 0; i < N; i++) {
				if( counter[i] == 0 )
				{
					// l_i reached
					complement_of_each_batch[m][size_of_each_complement_batch[m]] = i;
					size_of_each_complement_batch[m] += 1;
				}
				else
				{
					last_isogeny[0] = i;
					// l_i not reached
					batches[m][size_of_each_batch[m]] = i;
					size_of_each_batch[m] += 1;
				};
			}
		}

		// ＋とーでそれぞれ点を入手
		elligator(current_T[1], current_T[0], current_A[0]);
		//追加
		elligator(current_T[3], current_T[2], current_A[1]);

		// 4つの点を4倍する．
		// T_{-}
		yDBL(current_T[0], current_T[0], current_A[0]); // mult. by [2]
		yDBL(current_T[0], current_T[0], current_A[0]); // mult. by [2]
		// T_{+}
		yDBL(current_T[1], current_T[1], current_A[0]); // mult. by [2]
		yDBL(current_T[1], current_T[1], current_A[0]); // mult. by [2]
		
		//追加
		// T_{-}
		yDBL(current_T[2], current_T[2], current_A[1]); // mult. by [2]
		yDBL(current_T[2], current_T[2], current_A[1]); // mult. by [2]
		// T_{+}
		yDBL(current_T[3], current_T[3], current_A[1]); // mult. by [2]
		yDBL(current_T[3], current_T[3], current_A[1]); // mult. by [2]

		// Now, it is required to multiply by the complement of the batch
		for(i = 0; i < size_of_each_complement_batch[m]; i++)
		{
			yMUL(current_T[0], current_T[0], current_A[0], complement_of_each_batch[m][i]);	// Corresponding with T_{-}
			yMUL(current_T[1], current_T[1], current_A[0], complement_of_each_batch[m][i]);	// Corresponding with T_{+}

			//追加
			yMUL(current_T[2], current_T[2], current_A[1], complement_of_each_batch[m][i]);	// Corresponding with T_{-}
			yMUL(current_T[3], current_T[3], current_A[1], complement_of_each_batch[m][i]);	// Corresponding with T_{+}
		};
		for(i = 0; i < size_of_each_batch[m]; i++)
		{
			if( finished[batches[m][i]] == 1 )
			{ 
				//depends only on randomness
				continue;
			}
			else
			{
				ec = lookup(batches[m][i], tmp_e);	// To get current e_i in constant-time
				//先にbcがいる．(本当か？そもそも，bcの定義？)
				bc = isequal(ec >> 1, 0) & 1;		// Bit that determines if a dummy operation will be perfomed

				//ここから試し
				fp_cswap(current_T[0][0], current_T[1][0], (ec & 1));		// constant-time swap: T_{+} or T_{-}, that is the question.
				fp_cswap(current_T[0][1], current_T[1][1], (ec & 1));
				fp_cswap(current_T[2][0], current_T[3][0], (ec & 1));		// constant-time swap: T_{+} or T_{-}, that is the question.
				fp_cswap(current_T[2][1], current_T[3][1], (ec & 1));		// constant-time swap: T_{+} or T_{-}, that is the question.
				
				//追加　Gを入れ替えて今から計算しないといけないので，Aも入れ替えないといけない．
				fp_cswap(current_A[0][0], current_A[1][0], bc);	// constant-time swap: dummy or not dummy, that is the question.
				fp_cswap(current_A[0][1], current_A[1][1], bc);	// constant-time swap: dummy or not dummy, that is the question.
				
				//追加 bcによって曲線を入れ替えてしまったので，これが必要か．　
				fp_cswap(current_T[0][0], current_T[2][0], bc);		// constant-time swap: dummy or not dummy, that is the question.
				fp_cswap(current_T[0][1], current_T[2][1], bc);		// constant-time swap: dummy or not dummy, that is the question.
				fp_cswap(current_T[1][0], current_T[3][0], bc);		// constant-time swap: dummy or not dummy, that is the question.
				fp_cswap(current_T[1][1], current_T[3][1], bc);		// constant-time swap: dummy or not dummy, that is the question.

				// ここまで試し
				// Now, a degree-(l_{batches[m][i]}) will be constructed. Let l = l_{batches[m][i]}.
				point_copy(G[0], current_T[0]);	// order-l point determined by T_{-}
                    
				for (j = (i + 1); j < size_of_each_batch[m]; j++)
				{
					if( finished[batches[m][j]] == 0 )
					{
						//depends only on randomness
						yMUL(G[0], G[0], current_A[0], batches[m][j]);	// Corresponding with T_{-}
					};
				};

                //使わない曲線上の点をl倍しなければ，次のループでうまく使えないのでは？
				//Aを入れ替えたので，Tも入れ替える必要がある．
				yMUL(current_T[2], current_T[2], current_A[1], batches[m][i]);	// [l]T[1]
				yMUL(current_T[3], current_T[3], current_A[1], batches[m][i]);	// [l]T[1]

				if (isinfinity(G[0]) != 1)	// Depending on randomness
				{

					if ( isequal(batches[m][i], last_isogeny[m]) == 0)	// constant-time ask: just for avoiding the last isogeny evaluation
					{
                        yMUL(current_T[1], current_T[1], current_A[0], batches[m][i]);	// [l]T[1]
                        yISOG(K, current_A[0], G[0], current_A[0], batches[m][i]);

						yEVAL(current_T[0], current_T[0], K, batches[m][i]);	// evaluation of T[0]
						yEVAL(current_T[1], current_T[1], K, batches[m][i]);	// evaluation of T[0]

						

					}else{
                        yISOG(K, current_A[0], G[0], current_A[0], batches[m][i]);
                    };

					tmp_e[batches[m][i]] = (((ec >> 1) - (bc ^ 1)) << 1) ^ (ec & 0x1);
					counter[batches[m][i]] -= 1;
					isog_counter += 1;
				}
				else
				{
					// We must perform two scalar multiplications by l.                    
					yMUL(current_T[1], current_T[1], current_A[0], batches[m][i]);
				};

				fp_cswap(current_T[0][0], current_T[2][0], bc);		// constant-time swap: dummy or not dummy, that is the question.
				fp_cswap(current_T[0][1], current_T[2][1], bc);		// constant-time swap: dummy or not dummy, that is the question.
				fp_cswap(current_T[1][0], current_T[3][0], bc);		// constant-time swap: dummy or not dummy, that is the question.
				fp_cswap(current_T[1][1], current_T[3][1], bc);		// constant-time swap: dummy or not dummy, that is the question.

				fp_cswap(current_T[0][0], current_T[1][0], (ec & 1));		// constant-time swap: dummy or not dummy, that is the question.
				fp_cswap(current_T[0][1], current_T[1][1], (ec & 1));		// constant-time swap: dummy or not dummy, that is the question.
				//追加
				fp_cswap(current_T[2][0], current_T[3][0], (ec & 1));		// constant-time swap: dummy or not dummy, that is the question.
				fp_cswap(current_T[2][1], current_T[3][1], (ec & 1));		// constant-time swap: dummy or not dummy, that is the question.

				//追加
				fp_cswap(current_A[0][0], current_A[1][0], bc);	// constant-time swap: dummy or not dummy, that is the question.
				fp_cswap(current_A[0][1], current_A[1][1], bc);	// constant-time swap: dummy or not dummy, that is the question.
				
				if( counter[batches[m][i]] == 0 )
				{	
					//depends only on randomness

					finished[batches[m][i]] = 1;
					complement_of_each_batch[m][size_of_each_complement_batch[m]] = batches[m][i];
					size_of_each_complement_batch[m] += 1;
				};
			};
		};
		count += 1;
	};

	
	// --------------------------------------------------------------------------------------------------------	
	point_copy(C1, current_A[0]);
    point_copy(C2, current_A[1]);
};
