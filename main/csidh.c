#include<time.h>

#include "fp.h"
#include "edwards_curve.h"

// Measuring the perfomance
static uint64_t get_cycles()
{
   uint32_t lo, hi;
   asm volatile("rdtsc":"=a"(lo),"=d"(hi));
   return ((uint64_t)hi<<32) | lo;
};

#ifdef PROPOSAL_1
static uint8_t csidh(proj out1, proj out2, const uint8_t sk[], const proj in1, const proj in2, const uint8_t shared_flag)
{
	FP_ADD_COMPUTED = 0;
	FP_SQR_COMPUTED = 0;
	FP_MUL_COMPUTED = 0;

	if (!validate(in1)) {
		printf("supersingularではない\n");
		return 0;
	};
    if (!validate(in2)) {
        printf("supersingularではない\n");
		return 0;
	};

	action_evaluation(out1, out2, sk, in1, in2);
    if (shared_flag==1){
        //二つの出力を先に逆元計算で統一したものとし，それを和にしたものをout1
        fp_inv(out1[1]);
        fp_mul(out1[0],out1[0],out1[1]);
        fp_inv(out2[1]);
        fp_mul(out2[0],out2[0],out2[1]);
        fp_add(out1[0],out1[0],out2[0]);

    }
	return 1;
};
int main()
{
	
	uint64_t c0, c1;

	fp u;
	uint8_t sk_alice[N],	// secret key of Alice
	         sk_bob[N];	// secret key of Bob

	
	// ---
	proj random_E1, random_E2, tmp_E;
	uint8_t key[N];
	random_key(key);
	
	csidh(random_E1, random_E2, key, E, E, 0);

	// ---
	printf("NOTE: all the arithmetic is in Montgomery domain. ");
	printf("In addition, the ordering of the prime factors l_i's is:\n");
	printf("L := {\t  %3d", L[0]);

	for(int i = 1; i < N; i++)
	{
		printf(", %3d", L[i]);
		if( (i % 18) == 17 )
			printf("\n\t");
	};
	printf("};\n");

	printf("Moreover, the secret exponents belong to { -B_i, ..., 0, ..., B_i } where\n");
	printf("B := {\t  %3d", B[0]);
	for(int i = 1; i < N; i++)
	{
		printf(", %3d", B[i]);
		if( (i % 18) == 17 )
			printf("\n\t");
	};
	printf("};\n");
	printf("The dummy-free case implies that the secret exponents satisfy (e_i mod 2) = (b_i mod 2).\n");

	printf("-----------------------------------------------------------------------------------------------------\n");
	printf("First step: Alice and Bob random generate a secret key, and both compute and send their public curves\n\n");
	// Alice: random key generation
	random_key(sk_alice);
	printf_key(sk_alice, "sk_alice");

	proj E_alice1, E_alice2;

	// fp ss_E;
	// fp_inv(E[1]);
	// fp_mul(ss_E, E[0], E[1]);
	// fp_print(ss_E, NUMBER_OF_WORDS, 0, "ss_E");
	// fp_inv(E[1]);

	c0 = get_cycles();
	assert(csidh(E_alice1, E_alice2, sk_alice, E, E, 0));
	c1 = get_cycles();
	fp_print(E_alice1[0], NUMBER_OF_WORDS, 0, "E_alice1_a ");
	fp_print(E_alice1[1], NUMBER_OF_WORDS, 0, "E_alice1_ad");

    fp_print(E_alice2[0], NUMBER_OF_WORDS, 0, "E_alice2_a ");
	fp_print(E_alice2[1], NUMBER_OF_WORDS, 0, "E_alice2_ad");

	// fp_inv(E_alice[1]);
	// fp_mul(ss_E, E_alice[0], E_alice[1]);
	// fp_print(ss_E, NUMBER_OF_WORDS, 0, "ss_Alice");
	// fp_inv(E_alice[1]);

	printf("clock cycles: %3.03lf\n", ( 1.0 * (c1 - c0)) / (1000000.0));
	printf("Number of field operations computed: (%lu)M + (%lu)S + (%lu)a\n\n", FP_MUL_COMPUTED, FP_SQR_COMPUTED, FP_ADD_COMPUTED);

	// Bob: random key generation
	random_key(sk_bob);
	printf_key(sk_bob, "sk_bob");

	proj E_bob1, E_bob2;
	c0 = get_cycles();
	assert(csidh(E_bob1, E_bob2, sk_bob, E, E, 0));
	c1 = get_cycles();
	fp_print(E_bob1[0], NUMBER_OF_WORDS, 0, "E_bob1_a ");
	fp_print(E_bob1[1], NUMBER_OF_WORDS, 0, "E_bob1_ad");

    fp_print(E_bob2[0], NUMBER_OF_WORDS, 0, "E_bob2_a ");
	fp_print(E_bob2[1], NUMBER_OF_WORDS, 0, "E_bob2_ad");

	// fp_inv(E_bob[1]);
	// fp_mul(ss_E, E_bob[0], E_bob[1]);
	// fp_print(ss_E, NUMBER_OF_WORDS, 0, "ss_Bob");
	// fp_inv(E_bob[1]);

	printf("clock cycles: %3.03lf\n", ( 1.0 * (c1 - c0)) / (1000000.0));
	printf("Number of field operations computed: (%lu)M + (%lu)S + (%lu)a\n", FP_MUL_COMPUTED, FP_SQR_COMPUTED, FP_ADD_COMPUTED);


	printf("\n");	
	printf("------------------------------------------------------------------------------------------------------------\n");
	printf("Second step: Alice a Bob compute the shared secret by using the public curves of Bob and Alice, respectively\n");
	// Alice: shared secret
	proj ss_alice1, ss_alice2;
	c0 = get_cycles();
	assert(csidh(ss_alice1, ss_alice2, sk_alice, E_bob1, E_bob2, 1));
	c1 = get_cycles();
	fp_print(ss_alice1[0], NUMBER_OF_WORDS, 0, "ss_alice ");
	// fp_print(ss_alice1[1], NUMBER_OF_WORDS, 0, "ss_alice_ad");
	printf("clock cycles: %3.03lf\n", ( 1.0 * (c1 - c0)) / (1000000.0));
	printf("Number of field operations computed: (%lu)M + (%lu)S + (%lu)a\n", FP_MUL_COMPUTED, FP_SQR_COMPUTED, FP_ADD_COMPUTED);

	printf("\n");
	// Bob: shared secret
	proj ss_bob1, ss_bob2;
	c0 = get_cycles();
	assert(csidh(ss_bob1, ss_bob2, sk_bob, E_alice1, E_alice2, 1));
	c1 = get_cycles();
	fp_print(ss_bob1[0], NUMBER_OF_WORDS, 0, "ss_bob ");
	// fp_print(ss_bob1[1], NUMBER_OF_WORDS, 0, "ss_bob_ad");
	printf("clock cycles: %3.03lf\n", ( 1.0 * (c1 - c0)) / (1000000.0));
	printf("Number of field operations computed: (%lu)M + (%lu)S + (%lu)a\n", FP_MUL_COMPUTED, FP_SQR_COMPUTED, FP_ADD_COMPUTED);
	
	printf("\n");
	printf("------------------------------------------------------------------------------------------------------------\n");
	printf("At the end of the protocol, Alice and Bob have different but isomorphic Edwards curves. In other words, the\n");
	printf("Montgomery curve isomorphic to each one is the same. Thus, (ss_alice_a / ss_alice_ad) = (ss_bob_a / ss_bob_ad).\n");

	// fp ss_a, ss_b;
	// fp_inv(ss_alice[1]);
	// fp_mul(ss_a, ss_alice[0], ss_alice[1]);

	// fp_inv(ss_bob[1]);
	// fp_mul(ss_b, ss_bob[0], ss_bob[1]);
	// fp_print(ss_a, NUMBER_OF_WORDS, 0, "ss_a");
	// fp_print(ss_b, NUMBER_OF_WORDS, 0, "ss_b");

	if( compare(ss_alice1[0], ss_bob1[0], NUMBER_OF_WORDS) != 0 )
	{
		printf("\x1b[31m    _ ___    __ _     _        __ __ __ _  __ _____    __    _  _  __ _  \x1b[0m\n");
		printf("\x1b[31m|\\|/ \\ |    |_ / \\| ||_||     (_ |_ /  |_)|_ /   |    (_ |_||_||_)|_ | \\ \x1b[0m\n");
		printf("\x1b[31m| |\\_/ |    |__\\_X|_|| ||__   __)|__\\__| \\|__\\__ |    __)| || || \\|__|_/ \x1b[0m\n");
	}
	else
	{
		printf("\x1b[32m __ _     _        __ __ __ _  __ _____    __    _  _  __ _ \x1b[0m\n");
		printf("\x1b[32m|_ / \\| ||_||     (_ |_ /  |_)|_ /   |    (_ |_||_||_)|_ | \\ \x1b[0m\n");
		printf("\x1b[32m|__\\_X|_|| ||__   __)|__\\__| \\|__\\__ |    __)| || || \\|__|_/ \x1b[0m\n");
	}
	printf("\n");

	return 0;
};
#else


static uint8_t csidh(proj out, const uint8_t sk[], const proj in)
{
	FP_ADD_COMPUTED = 0;
	FP_SQR_COMPUTED = 0;
	FP_MUL_COMPUTED = 0;

	if (!validate(in)) {
		printf("supersingularではない\n");
		return 0;
	};

	action_evaluation(out, sk, in);
	return 1;
};

int main()
{
	
	uint64_t c0, c1;

	fp u;
	uint8_t sk_alice[N],	// secret key of Alice
	         sk_bob[N];	// secret key of Bob

	
	// ---
	proj random_E, tmp_E;
	uint8_t key[N];
	random_key(key);
	
	csidh(random_E, key, E);

	// ---
	printf("NOTE: all the arithmetic is in Montgomery domain. ");
	printf("In addition, the ordering of the prime factors l_i's is:\n");
	printf("L := {\t  %3d", L[0]);

	for(int i = 1; i < N; i++)
	{
		printf(", %3d", L[i]);
		if( (i % 18) == 17 )
			printf("\n\t");
	};
	printf("};\n");

	printf("Moreover, the secret exponents belong to { -B_i, ..., 0, ..., B_i } where\n");
	printf("B := {\t  %3d", B[0]);
	for(int i = 1; i < N; i++)
	{
		printf(", %3d", B[i]);
		if( (i % 18) == 17 )
			printf("\n\t");
	};
	printf("};\n");
	printf("The dummy-free case implies that the secret exponents satisfy (e_i mod 2) = (b_i mod 2).\n");

	printf("-----------------------------------------------------------------------------------------------------\n");
	printf("First step: Alice and Bob random generate a secret key, and both compute and send their public curves\n\n");
	// Alice: random key generation
	random_key(sk_alice);
	printf_key(sk_alice, "sk_alice");

	proj E_alice;

	// fp ss_E;
	// fp_inv(E[1]);
	// fp_mul(ss_E, E[0], E[1]);
	// fp_print(ss_E, NUMBER_OF_WORDS, 0, "ss_E");
	// fp_inv(E[1]);

	c0 = get_cycles();
	assert(csidh(E_alice, sk_alice, E));
	c1 = get_cycles();
	fp_print(E_alice[0], NUMBER_OF_WORDS, 0, "E_alice_a ");
	fp_print(E_alice[1], NUMBER_OF_WORDS, 0, "E_alice_ad");

	// fp_inv(E_alice[1]);
	// fp_mul(ss_E, E_alice[0], E_alice[1]);
	// fp_print(ss_E, NUMBER_OF_WORDS, 0, "ss_Alice");
	// fp_inv(E_alice[1]);

	printf("clock cycles: %3.03lf\n", ( 1.0 * (c1 - c0)) / (1000000.0));
	printf("Number of field operations computed: (%lu)M + (%lu)S + (%lu)a\n\n", FP_MUL_COMPUTED, FP_SQR_COMPUTED, FP_ADD_COMPUTED);

	// Bob: random key generation
	random_key(sk_bob);
	printf_key(sk_bob, "sk_bob");

	proj E_bob;
	c0 = get_cycles();
	assert(csidh(E_bob, sk_bob, E));
	c1 = get_cycles();
	fp_print(E_bob[0], NUMBER_OF_WORDS, 0, "E_bob_a ");
	fp_print(E_bob[1], NUMBER_OF_WORDS, 0, "E_bob_ad");

	// fp_inv(E_bob[1]);
	// fp_mul(ss_E, E_bob[0], E_bob[1]);
	// fp_print(ss_E, NUMBER_OF_WORDS, 0, "ss_Bob");
	// fp_inv(E_bob[1]);

	printf("clock cycles: %3.03lf\n", ( 1.0 * (c1 - c0)) / (1000000.0));
	printf("Number of field operations computed: (%lu)M + (%lu)S + (%lu)a\n", FP_MUL_COMPUTED, FP_SQR_COMPUTED, FP_ADD_COMPUTED);


	printf("\n");	
	printf("------------------------------------------------------------------------------------------------------------\n");
	printf("Second step: Alice a Bob compute the shared secret by using the public curves of Bob and Alice, respectively\n");
	// Alice: shared secret
	proj ss_alice;
	c0 = get_cycles();
	assert(csidh(ss_alice, sk_alice, E_bob));
	c1 = get_cycles();
	fp_print(ss_alice[0], NUMBER_OF_WORDS, 0, "ss_alice_a ");
	fp_print(ss_alice[1], NUMBER_OF_WORDS, 0, "ss_alice_ad");
	printf("clock cycles: %3.03lf\n", ( 1.0 * (c1 - c0)) / (1000000.0));
	printf("Number of field operations computed: (%lu)M + (%lu)S + (%lu)a\n", FP_MUL_COMPUTED, FP_SQR_COMPUTED, FP_ADD_COMPUTED);

	printf("\n");
	// Bob: shared secret
	proj ss_bob;
	c0 = get_cycles();
	assert(csidh(ss_bob, sk_bob, E_alice));
	c1 = get_cycles();
	fp_print(ss_bob[0], NUMBER_OF_WORDS, 0, "ss_bob_a ");
	fp_print(ss_bob[1], NUMBER_OF_WORDS, 0, "ss_bob_ad");
	printf("clock cycles: %3.03lf\n", ( 1.0 * (c1 - c0)) / (1000000.0));
	printf("Number of field operations computed: (%lu)M + (%lu)S + (%lu)a\n", FP_MUL_COMPUTED, FP_SQR_COMPUTED, FP_ADD_COMPUTED);
	
	printf("\n");
	printf("------------------------------------------------------------------------------------------------------------\n");
	printf("At the end of the protocol, Alice and Bob have different but isomorphic Edwards curves. In other words, the\n");
	printf("Montgomery curve isomorphic to each one is the same. Thus, (ss_alice_a / ss_alice_ad) = (ss_bob_a / ss_bob_ad).\n");

	fp ss_a, ss_b;
	fp_inv(ss_alice[1]);
	fp_mul(ss_a, ss_alice[0], ss_alice[1]);

	fp_inv(ss_bob[1]);
	fp_mul(ss_b, ss_bob[0], ss_bob[1]);
	fp_print(ss_a, NUMBER_OF_WORDS, 0, "ss_a");
	fp_print(ss_b, NUMBER_OF_WORDS, 0, "ss_b");

	if( compare(ss_a, ss_b, NUMBER_OF_WORDS) != 0 )
	{
		printf("\x1b[31m    _ ___    __ _     _        __ __ __ _  __ _____    __    _  _  __ _  \x1b[0m\n");
		printf("\x1b[31m|\\|/ \\ |    |_ / \\| ||_||     (_ |_ /  |_)|_ /   |    (_ |_||_||_)|_ | \\ \x1b[0m\n");
		printf("\x1b[31m| |\\_/ |    |__\\_X|_|| ||__   __)|__\\__| \\|__\\__ |    __)| || || \\|__|_/ \x1b[0m\n");
	}
	else
	{
		printf("\x1b[32m __ _     _        __ __ __ _  __ _____    __    _  _  __ _ \x1b[0m\n");
		printf("\x1b[32m|_ / \\| ||_||     (_ |_ /  |_)|_ /   |    (_ |_||_||_)|_ | \\ \x1b[0m\n");
		printf("\x1b[32m|__\\_X|_|| ||__   __)|__\\__| \\|__\\__ |    __)| || || \\|__|_/ \x1b[0m\n");
	}
	printf("\n");

	return 0;
};
#endif

