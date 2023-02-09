#include "edwards_curve_2.h"

#include "gmp.h"
#include "classgroup.h"

// #define VECTOR1
#define VECTOR2

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

#ifdef VECTOR2

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
	int8_t ec = 0, ec1=0, ec2=0;
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
	uint32_t bc1,bc2;

	

	// --------------------------------------------------------------------------------------------------------
    //step 1：babaiを動かして，新しいベクトルを得る．

    //初期化
    init_classgroup();

    mpz_t c;
    mpz_init_set_str(c,class_number,10);
    mpz_t sum;
    mpz_init(sum);

    int8_t k1,k2,k3,genzai,sgn1, sgn2;
					
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
    //mpz_neg(sum,sum);
    
    mpz_mod(sum,sum,c);

    int8_t vec[N];

    //イデアルをノルムの小さなベクトルに変換．
    mod_cn_2_vec(sum,vec);

    //---------------------------------------------------------------------------------------------------------
    //step2：新しいベクトルを作る．

	int8_t exp1;
    for (int ii=0; ii<=73; ii++){
		k1 = (int)( (2*(tmp_e[ii] & 0x1) - 1) * (tmp_e[ii] >> 1) );
		exp1=(int) -vec[ii];
		sgn1 = exp1 >> 7;	// sign of exp
		
		// Next, to write  key[i] = e || ((1 + sgn)/2)
		if (exp1!=0){
			cmov(&exp1, -exp1, sgn1 == -1);
		}
				
		vec[ii] = (exp1 << 1) ^ (1 & (1 + sgn1));
    }
	

	

	// --------------------------------------------------------------------------------------------------------
	//step3：前半のisogeny計算
	uint8_t m = 0, i, j;
	uint64_t number_of_batches = NUMBER_OF_BATCHES;
	
	uint8_t tmp_z[N];
	memcpy(tmp_z, vec, sizeof(uint8_t) * N);	// exponents

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
				ec1 = lookup(batches[m][i], tmp_e);	// To get current e_i in constant-time	
				ec2 = lookup(batches[m][i], tmp_z);

				bc1 = isequal(ec1 >> 1, 0) & 1;		// Bit that determines if a dummy operation will be perfomed	
				bc2	= isequal(ec2 >> 1, 0) & 1;

				ec = bc1 * ec2 + (bc1 ^ 1) * ec1;

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
					yISOG(K, current_A[0], G[0], current_A[0], batches[m][i]);

					//dummy-freeのアルゴリズム
					if ( isequal(batches[m][i], last_isogeny[m]) == 0)	// constant-time ask: just for avoiding the last isogeny evaluation
					{
						//どんな時でも使う．
						yEVAL(current_T[0], current_T[0], K, batches[m][i]);	// evaluation of T[0]
						yEVAL(current_T[1], current_T[1], K, batches[m][i]);	// evaluation of T[1]

						yMUL(current_T[1], current_T[1], current_A[0], batches[m][i]);	// [l]T[1]
					};

					// tmp_e[batches[m][i]] = ((((ec >> 1) - (bc ^ 1)) ^ bc) << 1) ^ ((ec & 0x1) ^ bc);
					//tmp_e[batches[m][i]] = ((((ec1 >> 1) - (bc ^ 1)) ^bc)<< 1) ^ (ec1 & 0x1);
					//tmp_z[batches[m][i]] = (((ec2 >> 1) - bc) << 1) ^ ((ec2 & 0x1) ^ bc);

					tmp_e[batches[m][i]] = (((ec1 >> 1) - (bc1 ^ 1)) << 1) ^ (ec1 & 0x1);
					tmp_z[batches[m][i]] = ((((ec2 >> 1) - (bc1 & (bc2 ^ 1))) ^ (bc1 & bc2)) << 1) ^ ((ec2 & 0x1) ^ (bc2 & bc1));
					counter[batches[m][i]] -= 1;
					isog_counter += 1;
					
					//追加．条件に達したときに中間の結果をもらう．
					if (isog_counter==start){
						point_copy(mid_A,current_A[0]);
						
						
						for (int ii=0;ii<=73;ii++){
							k1=(int)( (2*(tmp_e[ii] & 0x1) - 1) * (tmp_e[ii] >> 1) );
							k2=(int)( (2*(tmp_z[ii] & 0x1) - 1) * (tmp_z[ii] >> 1) );
							k3=(int)( (2*(vec[ii] & 0x1) - 1) * (vec[ii] >> 1) );
							exp1=k1-k3+k2;
							// exp1=k1-vec[ii];
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
	int8_t B1[N];
	memset(B1, 0, sizeof(int8_t) * N);



	
	for (int ii=0; ii<=73; ii++){
		k1 = (int)( (2*(tmp_e[ii] & 0x1) - 1) * (tmp_e[ii] >> 1) );
		
		k3 = (int)( (2*(tmp_z[ii] & 0x1) - 1) * (tmp_z[ii] >> 1) );
		exp1=k1+k3;
		// exp1=k1-vec[ii];
		sgn1 = (int8_t)exp1 >> 7;	// sign of exp
		cmov(&exp1,-exp1,sgn1==-1);
		tmp_e[ii] = (exp1 << 1) ^ (1 & (1 + sgn1));

		k2 = (int)( (2*(mid_r[ii] & 0x1) - 1) * (mid_r[ii] >> 1) );
		sgn2 = k2 >> 7;	// sign of exp
		cmov(&k2, -k2, sgn2 == -1);
		k3=B[ii]-exp1-k2-mid_counter[ii];
		sgn1=k3 >> 7;

		B1[ii]=B1[ii]+exp1+k2+(sgn1+1)*(k3 + (k3%2));

	}
	for (int ii=0; ii<=73; ii++){
		mid_r[ii]=B1[ii];

	}
};
#endif

#ifdef VECTOR1
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
};
#endif