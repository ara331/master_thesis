#include <math.h>
#include<time.h>

#include "fp.h"
#include "edwards_curve.h"

// const size_t stacksz = 0x8000;  /* 32k */

// Measuring the perfomance
static uint64_t get_cycles()
{
   uint32_t lo, hi;
   asm volatile("rdtsc":"=a"(lo),"=d"(hi));
   return ((uint64_t)hi<<32) | lo;
};

#ifdef PROPOSAL_1
static uint8_t csidh(proj out1, proj out2, const uint8_t sk[], const proj in1, const proj in2)
{
	FP_ADD_COMPUTED = 0;
	FP_SQR_COMPUTED = 0;
	FP_MUL_COMPUTED = 0;

	if (!validate(in1)) {
		return 0;
	};
    if (!validate(in2)) {
		return 0;
	};
	action_evaluation(out1, out2, sk, in1, in2);
	return 1;
};

unsigned long its = 5000;

int main()
{
	unsigned int i;

	double cc_min = 0xFFFFFFFFFFFFFFFF, cc_max = 0;
	double cc_mean = 0, cc_variance = 0;

	double cc_sample[its], cc_0, cc_1;

	// ---
	uint8_t key[N];
	proj random_E1, random_E2;

	point_copy(random_E1, E);
    point_copy(random_E2, E);

	for(i = 0; i < its; ++i)
	{

		if (its >= 100 && i % (its / 100) == 0) {
			printf("Doing %lu iterations of action with validation key:\t", its);
            printf("%2lu%%", 100 * i / its);
            fflush(stdout);
            printf("\r\x1b[K");
        }

		//printf("%d\n",i);

		random_key(key);
		cc_0 = get_cycles();
		assert(csidh(random_E1, random_E2, key, random_E1, random_E2));
		cc_1 = get_cycles();

		//fp_print(random_E[0],8,0,"a");
		//fp_print(random_E[1],8,0,"ad");

		
		cc_sample[i] = (double)(cc_1 - cc_0) / (1000000.0) ;
		// ---


		/**************************************/
		if(cc_min > cc_sample[i])
			cc_min = cc_sample[i];

		/**************************************/
		if(cc_max < cc_sample[i])
			cc_max = cc_sample[i];
		
		/**************************************/
		cc_mean += (double)cc_sample[i];
	};


	cc_mean = cc_mean / ((double)its * 1.0);

	for (i = 0; i < its; ++i)
	{
		cc_variance += (cc_sample[i] - cc_mean)*(cc_sample[i] - cc_mean);
	};

	cc_variance = cc_variance / ((double)its - 1.0);
		
	printf("\x1b[01;33mIterations: %lu\x1b[0m\n\n", its);

	printf("\x1b[33mAverage number of clock cycles: \x1b[32m %f \x1b[0m\n", cc_mean);
	printf("\x1b[33mStandar deviation of the number of clock cycles: \x1b[32m %f \x1b[0m\n", sqrt(cc_variance));
	printf("\x1b[33mMinimum of the number clock cycles: \x1b[32m %f \x1b[0m\n", cc_min);
	printf("\x1b[33mMaximum of the number of clock cycles: \x1b[32m %f \x1b[0m\n", cc_max);
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
		return 0;
	};
	action_evaluation(out, sk, in);
	return 1;
};

unsigned long its = 5000;

int main()
{
	unsigned int i;

	double cc_min = 0xFFFFFFFFFFFFFFFF, cc_max = 0;
	double cc_mean = 0, cc_variance = 0;

	double cc_sample[its], cc_0, cc_1;

	// ---
	uint8_t key[N];
	proj random_E;

	point_copy(random_E, E);

	// // メモリに関して？？
	// size_t bytes = 0;
    // unsigned char *stack;

	// __asm__ __volatile__ ("mov %%rsp, %0" : "=m"(stack));
    // stack -= stacksz;

	for(i = 0; i < its; ++i)
	{

		if (its >= 100 && i % (its / 100) == 0) {
			printf("Doing %lu iterations of action with validation key:\t", its);
            printf("%2lu%%", 100 * i / its);
            fflush(stdout);
            printf("\r\x1b[K");
        }

		//printf("%d\n",i);

		// メモリに関して？？
		/* spray stack */
        // unsigned char canary;
        // randombytes(&canary, 1);
        // for (size_t j = 0; j < stacksz; ++j)
        //     stack[j] = canary;

		random_key(key);
		cc_0 = get_cycles();
		assert(csidh(random_E, key, random_E));
		cc_1 = get_cycles();

		// //fp_print(random_E[0],8,0,"a");
		// //fp_print(random_E[1],8,0,"ad");

		// assert(*stack == canary);   /* make sure we sprayed enough */
        // for (size_t j = 0; j < stacksz - bytes; ++j)
        //     if (stack[j] != canary)
        //         bytes = stacksz - j;

		
		cc_sample[i] = (double)(cc_1 - cc_0) / (1000000.0) ;
		// ---


		/**************************************/
		if(cc_min > cc_sample[i])
			cc_min = cc_sample[i];

		/**************************************/
		if(cc_max < cc_sample[i])
			cc_max = cc_sample[i];
		
		/**************************************/
		cc_mean += (double)cc_sample[i];
	};


	cc_mean = cc_mean / ((double)its * 1.0);

	for (i = 0; i < its; ++i)
	{
		cc_variance += (cc_sample[i] - cc_mean)*(cc_sample[i] - cc_mean);
	};

	cc_variance = cc_variance / ((double)its - 1.0);
		
	printf("\x1b[01;33mIterations: %lu\x1b[0m\n\n", its);

	printf("\x1b[33mAverage number of clock cycles: \x1b[32m %f \x1b[0m\n", cc_mean);
	printf("\x1b[33mStandar deviation of the number of clock cycles: \x1b[32m %f \x1b[0m\n", sqrt(cc_variance));
	printf("\x1b[33mMinimum of the number clock cycles: \x1b[32m %f \x1b[0m\n", cc_min);
	printf("\x1b[33mMaximum of the number of clock cycles: \x1b[32m %f \x1b[0m\n", cc_max);
	//printf("stack memory usage: %lu b\n", bytes);
	printf("\n");

	return 0;
};

#endif