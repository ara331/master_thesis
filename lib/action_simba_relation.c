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
        //結局expは正の値になる。
		cmov(&exp, -exp, sgn == -1);
        //負なら後半は0、正なら後半は1
        //2*expの値は偶数。正なら2*epx+1,負なら2*expの値を取る。
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

void random_base(uint8_t *r)
{
    uint8_t i, tmp;
    // Random number from |[0, B]|
    randombytes(&tmp, 1);
    while ( issmaller((int32_t)N, (int32_t)tmp) == -1 )	// constant-time comparison
        randombytes(&tmp, 1);
    *r = tmp;
    //printf("%2d\n",*r);
}

void base_to_key(int8_t base[],uint8_t key[])
{
    uint8_t i;
	int8_t exp, sgn;
	for(i = 0; i < N; i++)
	{
        exp = base[i];
		sgn = exp >> 7;	// sign of exp

		// Next, to write  key[i] = e || ((1 + sgn)/2)
        //結局値は正にする。
		cmov(&exp, -exp, sgn == -1);
        //負なら後半は0、正なら後半は1
        //2*expの値は偶数。正なら2*epx+1,負なら2*expの値を取る。
		key[i] = (exp << 1) ^ (1 & (1 + sgn));
	};
}

void abs_key(uint8_t key[], uint8_t abs[])
{
    
    int8_t k,sgn;
    uint8_t i;
    for (i=0; i<N; i++)
    {
        //keyを整数値に変換
        k = (int)( (2*(key[i] & 0x1) - 1) * (key[i] >> 1) );
        sgn = k >> 7;	// sign of exp

		// Next, to write  key[i] = e || ((1 + sgn)/2)
		cmov(&k, -k, sgn == -1);
		abs[i] = k;
    }
}

void printf_base(uint8_t key[], uint8_t r , char *c)
{
    int i;
	printf("%s (%2d th):= ", c, r);
	// printf("{\t  %2d", base[0]);
    // //printf("{\t  %3d", (int)( (2*(base[0] & 0x1) - 1) * (base[0] >> 1) ));

	// for(i = 1; i < 74; i++)
	// {
	// 	printf(", %2d", base[i]);
    //     //printf(", %3d", (int)( (2*(base[i] & 0x1) - 1) * (base[i] >> 1) ) );
	// 	if( (i % 18) == 17 )
	// 		printf("\n\t\t");
	// };

	// printf("};\n"); 

    // int i;
	// printf("%s := ", c);
	printf("{\t  %3d", (int)( (2*(key[0] & 0x1) - 1) * (key[0] >> 1) ));

	for(i = 1; i < N; i++)
	{
		printf(", %3d", (int)( (2*(key[i] & 0x1) - 1) * (key[i] >> 1) ) );
		if( (i % 18) == 17 )
			printf("\n\t\t");
	};

	printf("};\n");
}


/* ----------------------------------------------------------------------------------------------- *
   action_evaluation()
   inputs: a the secret key, the Edwards curve constants A[0]:=a, and A[1]:=(a - d);
   output: the isogenous Edwards curve constants C[0]:=a' and C[1]:=(a' - d') determined by the action 
           evaluated at the secret key and public curve A
   
    NOTE: As far as we've understood how simba works; this next code implements simba approach. The 
          action computed by the next code uses two torsion points T_{+} and T_{-}, which their affine
          y-coordinates (in the isomorphic Montgomery curve) belongs to Fp and Fp2\Fp, respectively.
    ここがgroup actionの中身ととらえる。
    relation latticeはヘッダーファイルにおいてここで引用する。
 * ----------------------------------------------------------------------------------------------- */
void action_evaluation(proj C, const uint8_t key[], const proj A)
{
    //relation latticeを決定する。
    uint8_t r;
    uint64_t number_of_isogenies_2;
    random_base(&r);
    int8_t Base[N];
    uint8_t key2[N];
    memcpy(Base,Relation_lattice[r], sizeof(int8_t)*N);
    base_to_key(Base,key2);
    number_of_isogenies_2 = Number_of_isogenies_2[r];
    // printf_base(key2,r,"random_base");
    // printf_key(key2,"random_base");

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
	point_copy(current_A[0], A);			// initial Edwards curve constants a and (a -d)
    //printf("前：%ld \n",current_A[0][0][0]);
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
    //ここは簡単なやり方なら、eの絶対値を入れて、簡潔なコードにするなら|e|+|f|を入れる。
	//memcpy(counter, B, sizeof(int8_t) * N);		// At the beginning, we must perfomed b_i isogeny constructions for each l_i
    abs_key(tmp_e,counter);
    

	uint64_t isog_counter = 0;			// Total number of isogeny construction perfomed

	uint8_t last_isogeny[NUMBER_OF_BATCHES];
	//index for skipping point evaluations (the last one of each batch)
	memcpy(last_isogeny, LAST_ISOGENY, sizeof(uint8_t) * NUMBER_OF_BATCHES);
	uint32_t bc;

    // for (int iii=0; iii<N; iii++){
    //     printf("%2d,",counter[iii]);
    // }
    // printf("\n");

	// --------------------------------------------------------------------------------------------------------

	// --------------------------------------------------------------------------------------------------------
	// Main loop
	uint8_t m = 0, i, j;
	uint64_t number_of_batches = NUMBER_OF_BATCHES;
	uint32_t si;

    uint64_t number_of_isogenies_1=0;
    for (int ii=0; ii<N; ii++)
    {
        number_of_isogenies_1=number_of_isogenies_1+(uint64_t)(counter[ii]);
        if (counter[ii]==0){
            finished[ii]=1;
        }
    }
    // for (int ii=0; ii<N; ii++){
    //     printf("%2d,",finished[ii]);
    // }
    // printf("\n\n");

    // printf("number_of_isogenies(secret):%2ld \n",number_of_isogenies_1);

	//while (isog_counter < NUMBER_OF_ISOGENIES)
    while (isog_counter < number_of_isogenies_1)
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
				//point_copy(G[2], current_T[0]);	// T_{-}
				//point_copy(G[3], current_T[1]);	// T_{+}
                
                //ecが正か負か（偶数か奇数か）
				ec = lookup(batches[m][i], tmp_e);	// To get current e_i in constant-time
				fp_cswap(G[0][0], G[1][0], (ec & 1));		// constant-time swap: T_{+} or T_{-}, that is the question.
				fp_cswap(G[0][1], G[1][1], (ec & 1));		// constant-time swap: T_{+} or T_{-}, that is the question.
				//fp_cswap(G[2][0], G[3][0], (ec & 1));		// constant-time swap: T_{+} or T_{-}, that is the question.
				//fp_cswap(G[2][1], G[3][1], (ec & 1));		// constant-time swap: T_{+} or T_{-}, that is the question.

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
                    //0にいっちするなら1しないなら0
				    //bc = isequal(ec >> 1, 0) & 1;		// Bit that determines if a dummy operation will be perfomed
                    
					//既に終わっているなら、入れ替える。
					// fp_cswap(G[0][0], G[2][0], bc);		// constant-time swap: dummy or not dummy, that is the question.
					// fp_cswap(G[0][1], G[2][1], bc);		// constant-time swap: dummy or not dummy, that is the question.

					//yISOG(K, current_A[1], G[0], current_A[0], batches[m][i]);
                    //printf("前：%ld \n",current_A[0][0][0]);
                    yISOG(K, current_A[1], G[0], current_A[0], batches[m][i]);
                    //printf("後：%ld \n",current_A[0][0][0]);

					if ( isequal(batches[m][i], last_isogeny[m]) == 0)	// constant-time ask: just for avoiding the last isogeny evaluation
					{
						// mask = isequal(L[batches[m][i]], 3);		// Just for catching the case l = 3. This ask is done in constant-time
						// si = (L[batches[m][i]] >> 1);			// (l - 1) / 2

                        //こいつはisogenyをとおしてもl倍されないので、さきにl倍する。
                        yMUL(current_T[1], current_T[1], current_A[0], batches[m][i]);	// [l]T[1]

						//yEVAL(current_T[2], current_T[0], K, batches[m][i]);	// evaluation of T[0]
                        yEVAL(current_T[0], current_T[0], K, batches[m][i]);	// evaluation of T[0]
						//yEVAL(current_T[3], current_T[1], K, batches[m][i]);	// evaluation of T[0]
                        yEVAL(current_T[1], current_T[1], K, batches[m][i]);	// evaluation of T[0]

						// yADD(Z, K[(si + mask) - 1], G[0], K[(si + mask) - 2]);		// [(l + 1)/2]G[0]
						// fp_cswap(Z[0], K[si][0], mask ^ 1);				// constant-time swap: catching degree-3 isogeny case
						// fp_cswap(Z[1], K[si][1], mask ^ 1);				// constant-time swap: catching degree-3 isogeny case
						// yADD(current_T[0], K[si], K[si - 1], G[0]);			// [l]T[0]

                        //まだ終わってないなら入れ替える
						// fp_cswap(current_T[0][0], current_T[2][0], bc ^ 1);		// constant-time swap: dummy or not dummy, that is the question.
						// fp_cswap(current_T[0][1], current_T[2][1], bc ^ 1);		// constant-time swap: dummy or not dummy, that is the question.
						// fp_cswap(current_T[1][0], current_T[3][0], bc ^ 1);		// constant-time swap: dummy or not dummy, that is the question.
						// fp_cswap(current_T[1][1], current_T[3][1], bc ^ 1);		// constant-time swap: dummy or not dummy, that is the question.

					};

                    //まだ終わってないなら入れ替える
					fp_cswap(current_A[0][0], current_A[1][0], 1);	// constant-time swap: dummy or not dummy, that is the question.
					fp_cswap(current_A[0][1], current_A[1][1], 1);	// constant-time swap: dummy or not dummy, that is the question.

					tmp_e[batches[m][i]] = (((ec >> 1) - (1)) << 1) ^ (ec & 0x1);
					counter[batches[m][i]] -= 1;
					isog_counter += 1;
				}
				else
				{
					// We must perform two scalar multiplications by l.                    
					//yMUL(current_T[1], current_T[1], current_A[0], batches[m][i]);
                    yMUL(current_T[1], current_T[1], current_A[0], batches[m][i]);
				};

                //正か負でずれてたやつを元に戻す。
				fp_cswap(current_T[0][0], current_T[1][0], (ec & 1));		// constant-time swap: dummy or not dummy, that is the question.
				fp_cswap(current_T[0][1], current_T[1][1], (ec & 1));		// constant-time swap: dummy or not dummy, that is the question.
				
				if( counter[batches[m][i]] == 0 )
				{	
					//depends only on randomness

					finished[batches[m][i]] = 1;
					complement_of_each_batch[m][size_of_each_complement_batch[m]] = batches[m][i];
					size_of_each_complement_batch[m] += 1;
                    // for (int ii=0; ii<N; ii++){
                    //     printf("%2d,",finished[ii]);
                    // }
                    // printf("\n\n");
                    // printf_key(tmp_e,"tmp_e");
                    // printf("\n %2ld \n",isog_counter);
				};
			};
		};
		count += 1;
	};

	
	// --------------------------------------------------------------------------------------------------------	
	// proj C_;
    // point_copy(C_, current_A[0]);

    // for (int ii=0; ii<N; ii++){
    //     printf("%2d,",finished[ii]);
    // }
    // printf("\n");
    
    isog_counter = 0;
    count = 0;
    m = 0;
    ec=0;
    uint8_t tmp_z[N];
	memcpy(tmp_z, key2, sizeof(uint8_t) * N);	// exponents
    memset(finished, 0, sizeof(uint8_t) * N);
    memset(counter, 0, sizeof(int8_t) * N);
	abs_key(tmp_z,counter);
    //printf("finished secret key\n");
    // for (int ii=0; ii<N; ii++){
    //     printf("%2d,",counter[ii]);
    // }

    for (int ii=0; ii<N; ii++)
    {
        if (counter[ii]==0){
            finished[ii]=1;
        }
    }

    // for (int ii=0; ii<N; ii++){
    //     printf("%2d,",finished[ii]);
    // }

    // printf("%2ld \n",number_of_isogenies_2);

    // --------------------------------------------------------------------------------------------------------
	// SIMBA parameters
	// Batches
    number_of_batches = NUMBER_OF_BATCHES;

	for(uint8_t i = 0; i < NUMBER_OF_BATCHES; i++)
		memcpy(batches[i], BATCHES[i], sizeof(uint8_t) * SIZE_OF_EACH_BATCH[i]);

	memcpy(size_of_each_batch, SIZE_OF_EACH_BATCH, sizeof(uint8_t) * NUMBER_OF_BATCHES);
	// Complement of each batch

	memcpy(complement_of_each_batch, COMPLEMENT_OF_EACH_BATCH, sizeof(uint8_t) * NUMBER_OF_BATCHES * N);
	memcpy(size_of_each_complement_batch, SIZE_OF_EACH_COMPLEMENT_BATCH, sizeof(uint8_t) * NUMBER_OF_BATCHES);
	// --------------------------------------------------------------------------------------------------------


	//index for skipping point evaluations (the last one of each batch)
	memcpy(last_isogeny, LAST_ISOGENY, sizeof(uint8_t) * NUMBER_OF_BATCHES);
    

    while (isog_counter < number_of_isogenies_2)
	{
		m = (m + 1) % number_of_batches;
		
        //既にバッチは消え去っていると考えられる。
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
				// point_copy(G[2], current_T[0]);	// T_{-}
				// point_copy(G[3], current_T[1]);	// T_{+}
                
				ec = lookup(batches[m][i], tmp_z);	// To get current e_i in constant-time
				fp_cswap(G[0][0], G[1][0], (ec & 1));		// constant-time swap: T_{+} or T_{-}, that is the question.
				fp_cswap(G[0][1], G[1][1], (ec & 1));		// constant-time swap: T_{+} or T_{-}, that is the question.
				// fp_cswap(G[2][0], G[3][0], (ec & 1));		// constant-time swap: T_{+} or T_{-}, that is the question.
				// fp_cswap(G[2][1], G[3][1], (ec & 1));		// constant-time swap: T_{+} or T_{-}, that is the question.

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
					// bc = isequal(ec >> 1, 0) & 1;		// Bit that determines if a dummy operation will be perfomed
					
					// fp_cswap(G[0][0], G[2][0], bc);		// constant-time swap: dummy or not dummy, that is the question.
				    // fp_cswap(G[0][1], G[2][1], bc);		// constant-time swap: dummy or not dummy, that is the question.

					yISOG(K, current_A[1], G[0], current_A[0], batches[m][i]);

					if ( isequal(batches[m][i], last_isogeny[m]) == 0)	// constant-time ask: just for avoiding the last isogeny evaluation
					{
						// mask = isequal(L[batches[m][i]], 3);		// Just for catching the case l = 3. This ask is done in constant-time
						// si = (L[batches[m][i]] >> 1);			// (l - 1) / 2

						yMUL(current_T[1], current_T[1], current_A[0], batches[m][i]);	// [l]T[1]

						yEVAL(current_T[0], current_T[0], K, batches[m][i]);	// evaluation of T[0]
                        //yEVAL(current_T[0], current_T[0], K, batches[m][i]);	// evaluation of T[0]
						yEVAL(current_T[1], current_T[1], K, batches[m][i]);	// evaluation of T[0]
                        //yEVAL(current_T[1], current_T[1], K, batches[m][i]);	// evaluation of T[0]

					    // yADD(Z, K[(si + mask) - 1], G[0], K[(si + mask) - 2]);		// [(l + 1)/2]G[0]
						// fp_cswap(Z[0], K[si][0], mask ^ 1);				// constant-time swap: catching degree-3 isogeny case
						// fp_cswap(Z[1], K[si][1], mask ^ 1);				// constant-time swap: catching degree-3 isogeny case
						// yADD(current_T[0], K[si], K[si - 1], G[0]);			// [l]T[0]

						// fp_cswap(current_T[0][0], current_T[2][0], bc ^ 1);		// constant-time swap: dummy or not dummy, that is the question.
						// fp_cswap(current_T[0][1], current_T[2][1], bc ^ 1);		// constant-time swap: dummy or not dummy, that is the question.
						// fp_cswap(current_T[1][0], current_T[3][0], bc ^ 1);		// constant-time swap: dummy or not dummy, that is the question.
						// fp_cswap(current_T[1][1], current_T[3][1], bc ^ 1);		// constant-time swap: dummy or not dummy, that is the question.

					};

					//fp_cswap(current_A[0][0], current_A[1][0], bc ^ 1);	// constant-time swap: dummy or not dummy, that is the question.
					//fp_cswap(current_A[0][1], current_A[1][1], bc ^ 1);	// constant-time swap: dummy or not dummy, that is the question.
					fp_cswap(current_A[0][0], current_A[1][0], 1);
					fp_cswap(current_A[0][1], current_A[1][1], 1);

					tmp_z[batches[m][i]] = (((ec >> 1) - (bc ^ 1)) << 1) ^ (ec & 0x1);
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
                    // for (int ii=0; ii<N; ii++){
                    //     printf("%2d,",finished[ii]);
                    // }
                    // printf("\n\n");
                    // printf_key(tmp_e,"tmp_e");
                    // printf("\n %2ld \n",isog_counter);
				};
			};
		};
		count += 1;
	};

    // for (int ii=0; ii<N; ii++){
    //     printf("%2d,",finished[ii]);
    // }
    // printf("\n");
    point_copy(C, current_A[0]);
    // printf("%ld \n",C[0][0]);
    // printf("%ld \n",C_[0][0]);
};