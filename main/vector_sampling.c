#include<time.h>

#include "fp.h"
#include "edwards_curve_2.h"


/* defaults */

#ifndef BENCH_ITS
    #define BENCH_ITS 5000
#endif



const unsigned long its = BENCH_ITS;


// Measuring the perfomance
static uint64_t get_cycles()
{
   uint32_t lo, hi;
   asm volatile("rdtsc":"=a"(lo),"=d"(hi));
   return ((uint64_t)hi<<32) | lo;
};

static uint8_t csidh(uint8_t mid_r[], proj out, uint16_t start, const uint8_t sk[], const proj in)
{
	FP_ADD_COMPUTED = 0;
	FP_SQR_COMPUTED = 0;
	FP_MUL_COMPUTED = 0;

	if (!validate(in)) {
		return 0;
	};
	

	action_evaluation(mid_r, out, start, sk, in);
	return 1;
};

static void print_montgomery(proj E,char *m)
{
    //E[0]がaでE[1]がa-d
    fp d;
    proj E0;
    point_copy(E0,E);
    //d=a+(d-a)
    //(d=d-a)
    fp_sub(d,p,E0[1]);
    //d=d+a
    fp_add(d,d,E0[0]);
    fp A;
    //A/C=4A/4Cを求める
    //A=(a+d)/2
    //A=(a+d)
    fp_add(A,E0[0],d);
    //A=2(a+d)
    fp_add(A,A,A);
    //4C=a-d
    //(4C)^-1=(a-d)^-1
    fp_inv(E0[1]);
    fp_mul(A,A,E0[1]);
    fp_print(A, NUMBER_OF_WORDS, 0, m);
}


int main()
{

	uint64_t c0, c1;

	fp u;
	uint8_t priv[N];	// secret key of Bob
    int8_t mid_r[N];
    proj random_E;

    uint16_t mid_number_of_isogenies=305;

    int norm0=0;
    int norm1=0;

    for (uint16_t start=mid_number_of_isogenies;start<=mid_number_of_isogenies+50;start++){
        norm0=0;
        norm1=0;
        printf("%d個で止める",start);
        for (unsigned long i = 0; i < its; ++i) {

            if (its >= 100 && i % (its / 100) == 0) {
                printf("%2lu%%", 100 * i / its);
                fflush(stdout);
                printf("\r\x1b[K");
            }

            /**************************************/ 
            //csidhの秘密鍵をランダムに出力            
            random_key(priv);

            csidh(mid_r, random_E, start, priv, E);
            
            int n0=0;
            int n1=0;
            int8_t k1,k2,sgn;
            for (int ii=0;ii<=73;ii++){
                k1 = (int)( (2*(priv[ii] & 0x1) - 1) * (priv[ii] >> 1) );
                //k2 = (int)( (2*(mid_r[ii] & 0x1) - 1) * (mid_r[ii] >> 1) );
                k2=mid_r[ii];
                n0=n0+abs(k1);
                n1=n1+k2;
            }
            //秘密鍵の平均値
            norm0=norm0+n0;
            //途中で止まった場所と，本当の曲線の距離
            norm1=norm1+n1;            

            /**************************************/

        }
        //printf("iterations: %lu\n", its);
        printf("止めた場所： %d\n", start);
        //printf("秘密鍵のnormの平均: %lf⇒%lf\n", (double) norm0/ (double)its, ((double) norm0/ (double)its)/74);
        printf("途中で止まった場所からのベクトルのnormの平均: %lf⇒%lf\n", (double) norm1/ (double)its, ((double) norm1/ (double)its)/74);

    }



    

    return 0;
}

	
	