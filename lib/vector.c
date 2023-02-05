#include "edwards_curve_2.h"

#include "gmp.h"
#include "classgroup.h"

#define CCCDRS
//#define MCR

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

#ifdef MCR

/* ----------------------------------------------------------------------------------------------- *
   action_evaluation()
   inputs: a the secret key, the Edwards curve constants A[0]:=a, and A[1]:=(a - d);
   output: the isogenous Edwards curve constants C[0]:=a' and C[1]:=(a' - d') determined by the action 
           evaluated at the secret key and public curve A
   
    NOTE: As far as we've understood how simba works; this next code implements simba approach. The 
          action computed by the next code uses two torsion points T_{+} and T_{-}, which their affine
          y-coordinates (in the isomorphic Montgomery curve) belongs to Fp and Fp2\Fp, respectively.
 * ----------------------------------------------------------------------------------------------- */
void action_evaluation(int8_t mid_r[], proj C, uint16_t start, const uint8_t key[], const proj A)
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
	int8_t exp;

	proj current_A[2], current_T[4];
	point_copy(current_A[0], A);			// initial Edwards curve constants a and (a -d)
	// --------------------------------------------------------------------------------------------------------

	// --------------------------------------------------------------------------------------------------------
	// Vairables required for running SIMBA
	int8_t ec = 0, mask;
	uint16_t count = 0;
	proj G[4], K[(LARGE_L >> 1) + 1], Z;		// Current kernel
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
	uint32_t si;

	//追加．途中で止めてみる．その時のノルムを知りたい．
	
	proj mid_A;
	int8_t mid_counter[N];
	memset(mid_counter, 0, sizeof(int8_t) * N);

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

		// Before constructing isogenies, we must to search for suitable points
		elligator(current_T[1], current_T[0], current_A[0]);

		// Next, it is required to multiply the point by 4 and each l_i that doesn't belong to the current batch
		// T_{-}
		yDBL(current_T[0], current_T[0], current_A[0]); // mult. by [2]
		yDBL(current_T[0], current_T[0], current_A[0]); // mult. by [2]
		// T_{+}
		yDBL(current_T[1], current_T[1], current_A[0]); // mult. by [2]
		yDBL(current_T[1], current_T[1], current_A[0]); // mult. by [2]
		// Now, it is required to multiply by the complement of the batch
		for(i = 0; i < size_of_each_complement_batch[m]; i++)
		{
			yMUL(current_T[0], current_T[0], current_A[0], complement_of_each_batch[m][i]);	// Corresponding with T_{-}
			yMUL(current_T[1], current_T[1], current_A[0], complement_of_each_batch[m][i]);	// Corresponding with T_{+}
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
				// Now, a degree-(l_{batches[m][i]}) will be constructed. Let l = l_{batches[m][i]}.
				point_copy(G[0], current_T[0]);	// order-l point determined by T_{-}
				point_copy(G[1], current_T[1]);	// order-l point determined by T_{+}
				point_copy(G[2], current_T[0]);	// T_{-}
				point_copy(G[3], current_T[1]);	// T_{+}
                
				ec = lookup(batches[m][i], tmp_e);	// To get current e_i in constant-time
				//printf("%d,%d\n",ec,batches[m][i]);
				fp_cswap(G[0][0], G[1][0], (ec & 1));		// constant-time swap: T_{+} or T_{-}, that is the question.
				fp_cswap(G[0][1], G[1][1], (ec & 1));		// constant-time swap: T_{+} or T_{-}, that is the question.
				fp_cswap(G[2][0], G[3][0], (ec & 1));		// constant-time swap: T_{+} or T_{-}, that is the question.
				fp_cswap(G[2][1], G[3][1], (ec & 1));		// constant-time swap: T_{+} or T_{-}, that is the question.

				fp_cswap(current_T[0][0], current_T[1][0], (ec & 1));		// constant-time swap: T_{+} or T_{-}, that is the question.
				fp_cswap(current_T[0][1], current_T[1][1], (ec & 1));		// constant-time swap: T_{+} or T_{-}, that is the question.
                    
				for (j = (i + 1); j < size_of_each_batch[m]; j++)
				{
					if( finished[batches[m][j]] == 0 )
					{
						//depends only on randomness
						yMUL(G[0], G[0], current_A[0], batches[m][j]);	// Corresponding with T_{-}
					};
				};

				if (isinfinity(G[0]) != 1)	// Depending on randomness
				{
					bc = isequal(ec >> 1, 0) & 1;		// Bit that determines if a dummy operation will be perfomed
					
					//fp_cswap(G[0][0], G[2][0], bc);		// constant-time swap: dummy or not dummy, that is the question.
					//fp_cswap(G[0][1], G[2][1], bc);		// constant-time swap: dummy or not dummy, that is the question.


					yISOG(K, current_A[1], G[0], current_A[0], batches[m][i]);

					if ( isequal(batches[m][i], last_isogeny[m]) == 0)	// constant-time ask: just for avoiding the last isogeny evaluation
					{
						mask = isequal(L[batches[m][i]], 3);		// Just for catching the case l = 3. This ask is done in constant-time
						si = (L[batches[m][i]] >> 1);			// (l - 1) / 2

						yMUL(current_T[1], current_T[1], current_A[0], batches[m][i]);	// [l]T[1]

						yEVAL(current_T[2], current_T[0], K, batches[m][i]);	// evaluation of T[0]
						yEVAL(current_T[3], current_T[1], K, batches[m][i]);	// evaluation of T[0]

						yADD(Z, K[(si + mask) - 1], G[0], K[(si + mask) - 2]);		// [(l + 1)/2]G[0]
						fp_cswap(Z[0], K[si][0], mask ^ 1);				// constant-time swap: catching degree-3 isogeny case
						fp_cswap(Z[1], K[si][1], mask ^ 1);				// constant-time swap: catching degree-3 isogeny case
						yADD(current_T[0], K[si], K[si - 1], G[0]);			// [l]T[0]

						// fp_cswap(current_T[0][0], current_T[2][0], bc ^ 1);		// constant-time swap: dummy or not dummy, that is the question.
						// fp_cswap(current_T[0][1], current_T[2][1], bc ^ 1);		// constant-time swap: dummy or not dummy, that is the question.
						// fp_cswap(current_T[1][0], current_T[3][0], bc ^ 1);		// constant-time swap: dummy or not dummy, that is the question.
						// fp_cswap(current_T[1][1], current_T[3][1], bc ^ 1);		// constant-time swap: dummy or not dummy, that is the question.
						fp_cswap(current_T[0][0], current_T[2][0], 1);
						fp_cswap(current_T[0][1], current_T[2][1], 1);
						fp_cswap(current_T[1][0], current_T[3][0], 1);
						fp_cswap(current_T[1][1], current_T[3][1], 1);

					};

					//fp_cswap(current_A[0][0], current_A[1][0], bc ^ 1);	// constant-time swap: dummy or not dummy, that is the question.
					//fp_cswap(current_A[0][1], current_A[1][1], bc ^ 1);	// constant-time swap: dummy or not dummy, that is the question.

					fp_cswap(current_A[0][0], current_A[1][0], 1);
					fp_cswap(current_A[0][1], current_A[1][1], 1);

					tmp_e[batches[m][i]] = (((ec >> 1) - (bc ^ 1)) << 1) ^ (ec & 0x1);
					counter[batches[m][i]] -= 1;
					isog_counter += 1;
				}
				else
				{
					// We must perform two scalar multiplications by l.                    
					yMUL(current_T[1], current_T[1], current_A[0], batches[m][i]);
				};

				fp_cswap(current_T[0][0], current_T[1][0], (ec & 1));		// constant-time swap: dummy or not dummy, that is the question.
				fp_cswap(current_T[0][1], current_T[1][1], (ec & 1));		// constant-time swap: dummy or not dummy, that is the question.
				
				if( counter[batches[m][i]] == 0 )
				{	
					//depends only on randomness

					finished[batches[m][i]] = 1;
					complement_of_each_batch[m][size_of_each_complement_batch[m]] = batches[m][i];
					size_of_each_complement_batch[m] += 1;
				};

				//追加．条件に達したときに中間の結果をもらう．
				if (isog_counter==start){
					point_copy(mid_A,current_A[0]);
					memcpy(mid_counter, counter, sizeof(int8_t) * N);
					int8_t k,sgn,genzai;
					// printf("中間のcounter=[");
					// for (int ii=0;ii<=73;ii++){
					// 	printf("%d,",mid_counter[ii]);
					// 	//printf("%d,",sgn);
					// }
					// printf("]\n");
					//printf("補正のcounter=[");
					for (int ii=0;ii<=73;ii++){
						k = (int)( (2*(key[ii] & 0x1) - 1) * (key[ii] >> 1) );
						//秘密鍵が負の時はsgn=-1, 正の時はsgn=0
						sgn = (int8_t)k >> 7;
						genzai = (2*sgn+1) * (B[ii]-mid_counter[ii]);
						//いろいろ条件分けしたけれど，これでいいらしい．
						exp=k-genzai;
						sgn = (int8_t)exp >> 7;	// sign of exp
						cmov(&exp,-exp,sgn==-1);
						//printf("%d,",sgn);
						mid_r[ii] = (exp << 1) ^ (1 & (1 + sgn));
					}
					return;
					//printf("]\n");
                    // fp ss;
                    // fp_inv(mid_A[1]);
                    // fp_mul(ss, mid_A[0], mid_A[1]);
                    // fp_print(ss, NUMBER_OF_WORDS, 0, "曲線");
					//fp_print(mid_A[0], NUMBER_OF_WORDS, 0, "E_alice_a 中間");
	        		//fp_print(mid_A[1], NUMBER_OF_WORDS, 0, "E_alice_ad 中間");
				}
			};
		};
		count += 1;
	};

	
	// --------------------------------------------------------------------------------------------------------	
	point_copy(C, current_A[0]);
};
#endif

#ifdef CCCDRS
void action_evaluation(int8_t mid_r[], proj C, uint16_t start, const uint8_t key[], const proj A)
{
	// --------------------------------------------------------------------------------------------------------
    //step0：変数の宣言など

	// SIMBA parameters（一応simbaは使えるようにしている．ただし最適は知るところではない．）
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
    //秘密鍵のコピー
	uint8_t tmp_e[N];
	memcpy(tmp_e, key, sizeof(uint8_t) * N);	// exponents

	proj current_A[2], current_T[4];
    //初期曲線のコピー
	point_copy(current_A[0], A);			// initial Edwards curve constants a and (a -d)
	// --------------------------------------------------------------------------------------------------------

	// --------------------------------------------------------------------------------------------------------
	// Vairables required for running SIMBA
	int8_t ec = 0;
	uint16_t count = 0;
	proj G[4], K[(LARGE_L >> 1) + 1], Z;		// Current kernel

    //終わったかどうかを検知する．（今回は，二つのベクトルを使うので，どうしたものか，みたいなところはある．今までは違う政府のものが現れたらどうしてた？）
	uint8_t finished[N];				// flag that determines if the maximum number of isogeny constructions has been reached
	memset(finished, 0, sizeof(uint8_t) * N);

	int8_t counter[N];				// This variable determines how many isogeny construcctions has been perfomed
    //途中の進捗状況を表す変数．
	memset(counter, 0, sizeof(int8_t) * N);
	memcpy(counter, B, sizeof(int8_t) * N);		// At the beginning, we must perfomed b_i isogeny constructions for each l_i
	uint64_t isog_counter = 0;			// Total number of isogeny construction perfomed

    //追加．途中で止めてみる．その時のノルムを知りたい．
	proj mid_A;
	int8_t mid_counter[N];
	memset(mid_counter, 0, sizeof(int8_t) * N);


	uint8_t last_isogeny[NUMBER_OF_BATCHES];
	//index for skipping point evaluations (the last one of each batch)
	memcpy(last_isogeny, LAST_ISOGENY, sizeof(uint8_t) * NUMBER_OF_BATCHES);
	uint32_t bc;

	// --------------------------------------------------------------------------------------------------------
    //step 1：babaiを動かして，新しいベクトルを得る．

    //初期化
    init_classgroup();

    mpz_t c;
    mpz_init_set_str(c,class_number,10);
    mpz_t sum;
    mpz_init(sum);

    int8_t k1,k2,genzai,sgn1, sgn2;
					
    for (int ii=0; ii<=73; ii++){
        mpz_t d,e,mul;

        
        mpz_init(mul);
        mpz_init_set_str(d,dlogs[ii],10);
        k1 = (int)( (2*(tmp_e[ii] & 0x1) - 1) * (tmp_e[ii] >> 1) );
        mpz_init_set_si(e,k1);
        mpz_mul(mul,d,e);
        mpz_add(sum,sum,mul);
        
        mpz_mod(sum,sum,c);
        
    }
    mpz_neg(sum,sum);
    
    mpz_mod(sum,sum,c);

    int8_t vec[N];

    //イデアルをノルムの小さなベクトルに変換．
    mod_cn_2_vec(sum,vec);

    //---------------------------------------------------------------------------------------------------------
    //step2：新しいベクトルを作る．

	int8_t exp1,exp2;

	
    for (int ii=0; ii<=73; ii++){

        k1 = (int)( (2*(tmp_e[ii] & 0x1) - 1) * (tmp_e[ii] >> 1) );
        exp1=k1+vec[ii];
		sgn1 = exp1 >> 7;	// sign of exp
		// Next, to write  key[i] = e || ((1 + sgn)/2)
		cmov(&exp1, -exp1, sgn1 == -1);		
		tmp_e[ii] = (exp1 << 1) ^ (1 & (1 + sgn1));
        
    }

	// --------------------------------------------------------------------------------------------------------
	//step3：前半のisogeny計算
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
		

		// Before constructing isogenies, we must to search for suitable points
		elligator(current_T[1], current_T[0], current_A[0]);

		// Next, it is required to multiply the point by 4 and each l_i that doesn't belong to the current batch
		// T_{-}
		yDBL(current_T[0], current_T[0], current_A[0]); // mult. by [2]
		yDBL(current_T[0], current_T[0], current_A[0]); // mult. by [2]
		// T_{+}
		yDBL(current_T[1], current_T[1], current_A[0]); // mult. by [2]
		yDBL(current_T[1], current_T[1], current_A[0]); // mult. by [2]
		// Now, it is required to multiply by the complement of the batch
		for(i = 0; i < size_of_each_complement_batch[m]; i++)
		{
			yMUL(current_T[0], current_T[0], current_A[0], complement_of_each_batch[m][i]);	// Corresponding with T_{-}
			yMUL(current_T[1], current_T[1], current_A[0], complement_of_each_batch[m][i]);	// Corresponding with T_{+}
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

				fp_cswap(current_T[0][0], current_T[1][0], (ec & 1));		// constant-time swap: T_{+} or T_{-}, that is the question.
				fp_cswap(current_T[0][1], current_T[1][1], (ec & 1));		// constant-time swap: T_{+} or T_{-}, that is the question.
				
				point_copy(G[0], current_T[0]);	// order-l point determined by T_{-}
				
				for (j = (i + 1); j < size_of_each_batch[m]; j++)
				{
					if( finished[batches[m][j]] == 0 )
					{
						//depends only on randomness
						yMUL(G[0], G[0], current_A[0], batches[m][j]);	// Corresponding with T_{-}
					};
				};

				//if ((isinfinity(G[0]) != 1) && (isinfinity(G[1]) != 1))	// Depending on randomness
				if ((isinfinity(G[0]) != 1))	// Depending on randomness
				{
					bc = isequal(ec >> 1, 0) & 1;		// Bit that determines if a dummy operation will be perfomed

					//yISOG(K, current_A[0], G[0], current_A[0], batches[m][i]);

					//dummy-freeのアルゴリズム
					if ( isequal(batches[m][i], last_isogeny[m]) == 0)	// constant-time ask: just for avoiding the last isogeny evaluation
					{
						//どんな時でも使う．
						//yEVAL(current_T[0], current_T[0], K, batches[m][i]);	// evaluation of T[0]
						//yEVAL(current_T[1], current_T[1], K, batches[m][i]);	// evaluation of T[1]
						yMUL(current_T[0], current_T[0], current_A[0], batches[m][i]);
						yMUL(current_T[1], current_T[1], current_A[0], batches[m][i]);	// [l]T[1]
					};

					tmp_e[batches[m][i]] = ((((ec >> 1) - (bc ^ 1)) ^ bc) << 1) ^ ((ec & 0x1) ^ bc);
					counter[batches[m][i]] -= 1;
					isog_counter += 1;
					
					if (isog_counter==start){
						point_copy(mid_A,current_A[0]);
						
						
						for (int ii=0;ii<=73;ii++){
							k1=(int)( (2*(tmp_e[ii] & 0x1) - 1) * (tmp_e[ii] >> 1) );
							exp1=k1-vec[ii];
							sgn1 = (int8_t)exp1 >> 7;	// sign of exp
							cmov(&exp1,-exp1,sgn1==-1);

							// // Next, to write  key[i] = e || ((1 + sgn)/2)
							// cmov(&exp, -exp, sgn == -1);
							mid_r[ii] = (exp1 << 1) ^ (1 & (1 + sgn1));
							mid_counter[ii]=B[ii]-counter[ii];
						}

					}
				}
				else
				{
					// We must perform two scalar multiplications by l.                    
					yMUL(current_T[1], current_T[1], current_A[0], batches[m][i]);
				};

				fp_cswap(current_T[0][0], current_T[1][0], (ec & 1));		// constant-time swap: dummy or not dummy, that is the question.
				fp_cswap(current_T[0][1], current_T[1][1], (ec & 1));		// constant-time swap: dummy or not dummy, that is the question.
				
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

	
	int8_t k3;

	point_copy(current_A[1],mid_A);

	int8_t ec1 = 0, ec2=0,bc1=0,bc2=0;

	int8_t B1[N];
	memset(B1, 0, sizeof(int8_t) * N);
	//memset(mid_r, 0, sizeof(int8_t) * N);

	uint64_t number_of_isogenies_1=0;
	uint64_t number_of_isogenies_2=0;

	
	for (int ii=0; ii<=73; ii++){
		k1 = (int)( (2*(tmp_e[ii] & 0x1) - 1) * (tmp_e[ii] >> 1) );
		sgn1 = k1 >> 7;	// sign of exp
		k2 = (int)( (2*(mid_r[ii] & 0x1) - 1) * (mid_r[ii] >> 1) );
		sgn2 = k2 >> 7;	// sign of exp
		// Next, to write  key[i] = e || ((1 + sgn)/2)
		cmov(&k1, -k1, sgn1 == -1);
		cmov(&k2, -k2, sgn2 == -1);
		k3=B[ii]-k1-k2-mid_counter[ii];
		sgn1=k3 >> 7;

		number_of_isogenies_1=number_of_isogenies_1+k1;
		number_of_isogenies_2=number_of_isogenies_2+k2+(sgn1+1)*(k3 + (k3%2));
		B1[ii]=B1[ii]+k1+k2+(sgn1+1)*(k3 + (k3%2));

	}
	for (int ii=0; ii<=73; ii++){
		mid_r[ii]=B1[ii];

	}
	// // --------------------------------------------------------------------------------------------------------
    // //step0：変数の宣言など


	// // SIMBA parameters
	// // Batches
	// uint8_t batches[NUMBER_OF_BATCHES][SIZE_OF_EACH_BATCH[0]];
	// uint8_t size_of_each_batch[NUMBER_OF_BATCHES];

	// for(uint8_t i = 0; i < NUMBER_OF_BATCHES; i++)
	// 	memcpy(batches[i], BATCHES[i], sizeof(uint8_t) * SIZE_OF_EACH_BATCH[i]);

	// memcpy(size_of_each_batch, SIZE_OF_EACH_BATCH, sizeof(uint8_t) * NUMBER_OF_BATCHES);
	// // Complement of each batch
	// uint8_t complement_of_each_batch[NUMBER_OF_BATCHES][N];
	// uint8_t size_of_each_complement_batch[NUMBER_OF_BATCHES];

	// memcpy(complement_of_each_batch, COMPLEMENT_OF_EACH_BATCH, sizeof(uint8_t) * NUMBER_OF_BATCHES * N);
	// memcpy(size_of_each_complement_batch, SIZE_OF_EACH_COMPLEMENT_BATCH, sizeof(uint8_t) * NUMBER_OF_BATCHES);
	// // --------------------------------------------------------------------------------------------------------

	// // --------------------------------------------------------------------------------------------------------
	// // Copy of public and private data (the private key is modified each iteration)
    // //秘密鍵のコピー
	// uint8_t tmp_e[N];
	// memcpy(tmp_e, key, sizeof(uint8_t) * N);	// exponents

	// proj current_A[2], current_T[4];
    // //初期曲線のコピー
	// point_copy(current_A[0], A);			// initial Edwards curve constants a and (a -d)
	// // --------------------------------------------------------------------------------------------------------

	// // --------------------------------------------------------------------------------------------------------
	// // Vairables required for running SIMBA
	// int8_t ec = 0, mask;
	// uint16_t count = 0;
	// proj G[4], K[(LARGE_L >> 1) + 1], Z;		// Current kernel

    // //終わったかどうかを検知する．（今回は，二つのベクトルを使うので，どうしたものか，みたいなところはある．今までは違う政府のものが現れたらどうしてた？）
	// uint8_t finished[N];				// flag that determines if the maximum number of isogeny constructions has been reached
	// memset(finished, 0, sizeof(uint8_t) * N);

	// int8_t counter[N];				// This variable determines how many isogeny construcctions has been perfomed
    // //途中の進捗状況を表す変数．
	// memset(counter, 0, sizeof(int8_t) * N);
	// //memcpy(counter, B, sizeof(int8_t) * N);		// At the beginning, we must perfomed b_i isogeny constructions for each l_i
	// uint64_t isog_counter = 0;			// Total number of isogeny construction perfomed

    // //追加．途中で止めてみる．その時のノルムを知りたい．
	// proj mid_A;
	// int8_t mid_counter[N];
	// memset(mid_counter, 0, sizeof(int8_t) * N);

	// uint8_t last_isogeny[NUMBER_OF_BATCHES];
	// //index for skipping point evaluations (the last one of each batch)
	// memcpy(last_isogeny, LAST_ISOGENY, sizeof(uint8_t) * NUMBER_OF_BATCHES);
	// uint32_t bc;

	// // --------------------------------------------------------------------------------------------------------
    // //step 1：babaiを動かして，新しいベクトルを得る．

    // //よくわからん初期化
    // init_classgroup();

    // mpz_t c;
    // mpz_init_set_str(c,class_number,10);
    // mpz_t sum;
    // mpz_init(sum);

    // int8_t k,k1,k2,genzai,sgn1, sgn2;
					
    // for (int ii=0; ii<=73; ii++){
    //     mpz_t d,e,mul;

        
    //     mpz_init(mul);
    //     mpz_init_set_str(d,dlogs[ii],10);
    //     k = (int)( (2*(tmp_e[ii] & 0x1) - 1) * (tmp_e[ii] >> 1) );
    //     mpz_init_set_si(e,k);
    //     mpz_mul(mul,d,e);
    //     mpz_add(sum,sum,mul);
        
    //     mpz_mod(sum,sum,c);
        
    // }
    // //sumは，∑秘密*dlogで，g^sumのactionを表している．
    // //sum≢0 mod hならば，秘密はrelation latticeのベクトル．
    // //つまり，e+r=0を作りたいのなら，sum+dlog74*x=0なので，x=-sum*dlog74^-1
    // mpz_neg(sum,sum);
    // //gmp_printf("%Zd\n",sum);
    // mpz_t dlog73,x;
    // mpz_init(x);
    // mpz_init_set_str(dlog73,dlogs[72],10);
    // mpz_invert(x,dlog73,c);
    // mpz_mul(x,x,sum);
    // mpz_mul(x,dlog73,x);
    // mpz_mod(x,x,c);

    // //gmp_printf("%Zd\n",x);

    // int8_t vec[N];

    // //babaiを動かす．
    // mod_cn_2_vec(x,vec);

    // // printf("[");
    // // for(int ii=0;ii<=73;ii++){
    // //     printf("%d,",vec[ii]);
    // // }
    // // printf("]\n");
    // //---------------------------------------------------------------------------------------------------------
    // //step2：新しいベクトルを作る．
	// //いったん二つのベクトルを足し合わせる．
	// //5に併せずに，340回で切る．
	// //bc=1ならば，ベクトルに±1を付け足す．それで，奇数の時は増えて，偶数の時は減って同じになる．
	// //とりあえずdummyfreeのアルゴリズムを理解する．
	// int8_t exp1,exp2;
	// int8_t B2[N];
	// memset(B2, 0, sizeof(int8_t) * N);
	
    // for (int ii=0; ii<=73; ii++){

    //     k = (int)( (2*(tmp_e[ii] & 0x1) - 1) * (tmp_e[ii] >> 1) );
    //     exp1=k+vec[ii];
	// 	//exp = (int8_t)tmp;
	// 	// // Mapping integers from |[ 0, 2B |] into |[ -B, B ]|
	// 	// exp = exp - (int8_t)B[i];
	// 	sgn1 = exp1 >> 7;	// sign of exp
	// 	// Next, to write  key[i] = e || ((1 + sgn)/2)
	// 	cmov(&exp1, -exp1, sgn1 == -1);		
	// 	tmp_e[ii] = (exp1 << 1) ^ (1 & (1 + sgn1));
        
    // }
	// //int8_t counter1[N];
	// //memset(counter1, 0, sizeof(int8_t) * N);		// At the beginning, we must perfomed b_i isogeny constructions for each l_i

	// // printf("秘密鍵+α=[");
    // // for(int ii=0;ii<=73;ii++){
	// // 	k = (int)( (2*(tmp_e[ii] & 0x1) - 1) * (tmp_e[ii] >> 1) );
    // //     printf("%d,",k);
    // // }
    // // printf("]\n");
	// uint8_t tmp_z[N];
	// memcpy(tmp_z, tmp_e, sizeof(uint8_t) * N);

	



	// // --------------------------------------------------------------------------------------------------------
	// //step3：前半のisogeny計算
	// uint8_t m = 0, i, j;
	// uint64_t number_of_batches = NUMBER_OF_BATCHES;
	// uint32_t si;


	// while (isog_counter < NUMBER_OF_ISOGENIES)
	// {
	// 	m = (m + 1) % number_of_batches;
		
	// 	if(count == MY*number_of_batches) {  	//merge the batches after my rounds
	// 		m = 0;
	// 		size_of_each_complement_batch[m] = 0;
	// 		size_of_each_batch[m] = 0;
	// 		number_of_batches = 1;

	// 		for(i = 0; i < N; i++) {
	// 			if( counter[i] == B[i] )
	// 			{
	// 				// l_i reached
	// 				complement_of_each_batch[m][size_of_each_complement_batch[m]] = i;
	// 				size_of_each_complement_batch[m] += 1;
	// 			}
	// 			else
	// 			{
	// 				last_isogeny[0] = i;
	// 				// l_i not reached
	// 				batches[m][size_of_each_batch[m]] = i;
	// 				size_of_each_batch[m] += 1;
	// 			};
	// 		}
	// 	}
		

	// 	// Before constructing isogenies, we must to search for suitable points
	// 	elligator(current_T[1], current_T[0], current_A[0]);

	// 	// Next, it is required to multiply the point by 4 and each l_i that doesn't belong to the current batch
	// 	// T_{-}
	// 	yDBL(current_T[0], current_T[0], current_A[0]); // mult. by [2]
	// 	yDBL(current_T[0], current_T[0], current_A[0]); // mult. by [2]
	// 	// T_{+}
	// 	yDBL(current_T[1], current_T[1], current_A[0]); // mult. by [2]
	// 	yDBL(current_T[1], current_T[1], current_A[0]); // mult. by [2]
	// 	// Now, it is required to multiply by the complement of the batch
	// 	for(i = 0; i < size_of_each_complement_batch[m]; i++)
	// 	{
	// 		yMUL(current_T[0], current_T[0], current_A[0], complement_of_each_batch[m][i]);	// Corresponding with T_{-}
	// 		yMUL(current_T[1], current_T[1], current_A[0], complement_of_each_batch[m][i]);	// Corresponding with T_{+}
	// 	};

	// 	for(i = 0; i < size_of_each_batch[m]; i++)
	// 	{
	// 		if( finished[batches[m][i]] == 1 )
	// 		{ 
	// 			//depends only on randomness
	// 			continue;
	// 		}
	// 		else
	// 		{
	// 			// Now, a degree-(l_{batches[m][i]}) will be constructed. Let l = l_{batches[m][i]}.
	// 			point_copy(G[0], current_T[0]);	// order-l point determined by T_{-}
	// 			point_copy(G[1], current_T[1]);	// order-l point determined by T_{+}
                
	// 			ec = lookup(batches[m][i], tmp_e);	// To get current e_i in constant-time
	// 			fp_cswap(G[0][0], G[1][0], (ec & 1));		// constant-time swap: T_{+} or T_{-}, that is the question.
	// 			fp_cswap(G[0][1], G[1][1], (ec & 1));		// constant-time swap: T_{+} or T_{-}, that is the question.
	// 			//fp_cswap(G[2][0], G[3][0], (ec & 1));		// constant-time swap: T_{+} or T_{-}, that is the question.
	// 			//fp_cswap(G[2][1], G[3][1], (ec & 1));		// constant-time swap: T_{+} or T_{-}, that is the question.

	// 			fp_cswap(current_T[0][0], current_T[1][0], (ec & 1));		// constant-time swap: T_{+} or T_{-}, that is the question.
	// 			fp_cswap(current_T[0][1], current_T[1][1], (ec & 1));		// constant-time swap: T_{+} or T_{-}, that is the question.
                    
	// 			for (j = (i + 1); j < size_of_each_batch[m]; j++)
	// 			{
	// 				if( finished[batches[m][j]] == 0 )
	// 				{
	// 					//depends only on randomness
	// 					yMUL(G[0], G[0], current_A[0], batches[m][j]);	// Corresponding with T_{-}
	// 				};
	// 			};

	// 			//本当に両方必要か？諸説あり．しかし，そういうパターンもそうない気もする．
	// 			//if ((isinfinity(G[0]) != 1) && (isinfinity(G[1]) != 1))	// Depending on randomness
	// 			if ((isinfinity(G[0]) != 1))	// Depending on randomness
	// 			{
	// 				bc = isequal(ec >> 1, 0) & 1;		// Bit that determines if a dummy operation will be perfomed
					
	// 				//fp_cswap(G[0][0], G[2][0], bc);		// constant-time swap: dummy or not dummy, that is the question.
	// 				//fp_cswap(G[0][1], G[2][1], bc);		// constant-time swap: dummy or not dummy, that is the question.

					
	// 				//tmp_e[batches[m][i]] = (((ec >> 1) - (bc ^ 1)) << 1) ^ (ec & 0x1);
	// 				tmp_e[batches[m][i]] = ((((ec >> 1) - (bc ^ 1)) ^ bc) << 1) ^ ((ec & 0x1) ^ bc);
	// 				//counter[batches[m][i]] -= 1;
	// 				counter[batches[m][i]] += 1;
	// 				isog_counter += 1;
					
	// 				// printf("[");
	// 				// for(int ii=0;ii<=73;ii++){
	// 				// 	k = (int)( (2*(tmp_e[ii] & 0x1) - 1) * (tmp_e[ii] >> 1) );
	// 				// 	printf("%d,",k);
	// 				// }
	// 				// printf("]\n");
	// 				//追加．条件に達したときに中間の結果をもらう．
	// 				if (isog_counter==start){
	// 					point_copy(mid_A,current_A[0]);

	// 					for (int ii=0;ii<=73;ii++){
	// 						k1=(int)( (2*(tmp_e[ii] & 0x1) - 1) * (tmp_e[ii] >> 1) );
	// 						exp1=k1-vec[ii];
	// 						sgn1 = (int8_t)exp1 >> 7;	// sign of exp
	// 						cmov(&exp1,-exp1,sgn1==-1);

	// 						exp2=B[ii]-counter[ii];
	// 						sgn1 = (int8_t)exp2 >> 7;
	// 						cmov(&exp2,-exp2,sgn1==-1);

	// 						// // Next, to write  key[i] = e || ((1 + sgn)/2)
	// 						// cmov(&exp, -exp, sgn == -1);
	// 						mid_r[ii] = exp1+exp2;
							
	// 					}

	// 				}
						
					
	// 			}
	// 			else
	// 			{
	// 				// We must perform two scalar multiplications by l.                    
	// 				yMUL(current_T[1], current_T[1], current_A[0], batches[m][i]);
	// 			};

	// 			fp_cswap(current_T[0][0], current_T[1][0], (ec & 1));		// constant-time swap: dummy or not dummy, that is the question.
	// 			fp_cswap(current_T[0][1], current_T[1][1], (ec & 1));		// constant-time swap: dummy or not dummy, that is the question.
				
	// 			if( counter[batches[m][i]] == B[batches[m][i]] )
	// 			{	
	// 				//depends only on randomness

	// 				finished[batches[m][i]] = 1;
	// 				complement_of_each_batch[m][size_of_each_complement_batch[m]] = batches[m][i];
	// 				size_of_each_complement_batch[m] += 1;
	// 			};
                
	// 		};
	// 	};
	// 	count += 1;
	// };
	// point_copy(C, current_A[0]);
};
#endif