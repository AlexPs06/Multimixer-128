#include <iostream>
#include <stdint.h>
#include <wmmintrin.h>
#include <immintrin.h>
#include <emmintrin.h>
#include <pmmintrin.h>
#include <cstring>
#define ALIGN(n) __attribute__ ((aligned(n)))
#define pipeline 1
#define size_message 16777216   // Tamaño del mensaje a procesar


#define EXPAND_ASSIST(v1,v2,v3,v4,shuff_const,aes_const)                    \
    v2 = _mm_aeskeygenassist_si128(v4,aes_const);                           \
    v3 = _mm_castps_si128(_mm_shuffle_ps(_mm_castsi128_ps(v3),              \
                                         _mm_castsi128_ps(v1), 16));        \
    v1 = _mm_xor_si128(v1,v3);                                              \
    v3 = _mm_castps_si128(_mm_shuffle_ps(_mm_castsi128_ps(v3),              \
                                         _mm_castsi128_ps(v1), 140));       \
    v1 = _mm_xor_si128(v1,v3);                                              \
    v2 = _mm_shuffle_epi32(v2,shuff_const);                                 \
    v1 = _mm_xor_si128(v1,v2)

using namespace std;


void Multimixer(const uint8_t* input,uint8_t* tag,__m128i * keys, const uint64_t lenght);
void AES_128_Key_Expansion(const unsigned char *userkey, void *key);
static inline void AES_encrypt(__m128i tmp, __m128i *out,__m128i *key, int rounds);
static inline __m128i gf_reduce_128(__m128i x, __m128i y);
static inline __m128i AES_Encrypt_rounds(__m128i tmp, __m128i *key, int rounds);
static inline __m128i AES_Encrypt_rounds_static_keys(__m128i tmp, __m128i key, int rounds);


char infoString[]= "Para-Hash AVX128 i7-11700";  /* Each AE implementation must have a global one */

#define Nr 10   // Número de rondas para AES-128



static const __m128i RK_ZERO = _mm_setzero_si128();
static const __m128i RK_ONES = _mm_set1_epi8(0x01);


void imprimiArreglo(int tam, unsigned char *in )
{

    for (int i = 0; i<tam; i++){
        printf("%02x", in[i] );
    }
    printf("\n" );

}

void imprimiArreglo_32bits(int tam, unsigned int *in )
{

    for (int i = 0; i<tam; i++){
        printf("%08x ", in[i] );
    }
    printf("\n" );

}

void imprimiArreglo_64bits(int tam, unsigned long int *in )
{

    for (int i = 0; i<tam; i++){
        printf("%lx ", in[i] );
    }
    printf("\n" );

}


static inline void update_function(__m128i X,
                    __m128i Y,
                    __m128i * output)
{

        __m128i U,V,R[4];
        __m128i checksum_x = _mm_setzero_si128();
        __m128i checksum_y = _mm_setzero_si128();


	for (size_t i = 0; i < 4; i++)
	{
    	checksum_x = _mm_add_epi32(checksum_x, _mm_shuffle_epi32(X, _MM_SHUFFLE((0+i)%4, (1+i)%4, (2+i)%4, (3+i)%4))); 
    	checksum_y = _mm_add_epi32(checksum_y, _mm_shuffle_epi32(Y, _MM_SHUFFLE((0+i)%4, (1+i)%4, (2+i)%4, (3+i)%4))); 
	}
	

    U = _mm_sub_epi32(checksum_x,_mm_shuffle_epi32(X,_MM_SHUFFLE(2, 1, 0, 3))); 
    V = _mm_sub_epi32(checksum_y,_mm_shuffle_epi32(Y,_MM_SHUFFLE(3, 2, 1, 0)));  


    R[0] = _mm_mul_epu32(X,Y);
    R[1] = _mm_mul_epu32(_mm_srli_epi64(X,32),_mm_srli_epi64(Y,32));
    R[2] = _mm_mul_epu32(U,V);
	R[3] = _mm_mul_epu32(_mm_srli_epi64(U,32),_mm_srli_epi64(V,32));

    

	output[0] = _mm_add_epi64(R[0],output[0]);
	output[1] = _mm_add_epi64(R[1],output[1]);
	output[2] = _mm_add_epi64(R[2],output[2]);
	output[3] = _mm_add_epi64(R[3],output[3]);

    
}



static inline __m128i gf_reduce_128(__m128i x, __m128i y)
{

    // Obtained the low and high parts
    __m128i lo = _mm_unpacklo_epi64(x, y);
    __m128i hi = _mm_unpackhi_epi64(x, y);

    /* First fold: reduce x^64 terms */
    lo = _mm_xor_si128(lo,hi);
    lo = _mm_xor_si128(lo,_mm_slli_epi64(hi,1));
    lo = _mm_xor_si128(lo,_mm_slli_epi64(hi,3));
    lo = _mm_xor_si128(lo,_mm_slli_epi64(hi,4));

    /* Bits that overflowed beyond bit 63 */
   __m128i carry;
    carry  = _mm_srli_epi64(hi, 63);
    carry ^= _mm_srli_epi64(hi, 62);
    carry ^= _mm_srli_epi64(hi, 60);
    carry ^= _mm_srli_epi64(hi, 59);

    /* First fold: reduce x^64 terms */
    lo = _mm_xor_si128(lo,carry);
    lo = _mm_xor_si128(lo,_mm_slli_epi64(carry,1));
    lo = _mm_xor_si128(lo,_mm_slli_epi64(carry,3));
    lo = _mm_xor_si128(lo,_mm_slli_epi64(carry,4));

    return lo;
}


void Multimixer(const uint8_t* input,
                 uint8_t* tag,
                 __m128i * keys,
                 const uint64_t lenght)
{
 
    uint64_t i = 0;

    /*
     * Accumulators for the hash computation.
     * Each entry stores a 128-bit value split into two 64-bit lanes.
     */
    __m128i output[4];

    __m128i final_output[4];
    
    /*
     * Initialize accumulators to zero.
     */
    for (i = 0; i < 4; i++) {
        output[i] = _mm_setzero_si128();
    }

    /*
     * Compute the number of 128-bit input blocks.
     * Only full blocks are processed.
     */
    uint64_t size = 0;
    if (lenght % 16 == 0)
        size = lenght / 16; 

    /*
    * Load input blocks into NEON registers.
    */
    const __m128i* ptr = (const __m128i*)input;

    /*
     * Main processing loop:
     * Blocks are processed in pairs (X, Y).
     */
    
    for (i = 0; i < size - 1; i = i + 2) {

        /*
         * Load keys blocks into NEON registers.
         */
        __m128i generate_key_x = keys[i];
        __m128i generate_key_y = keys[i+1];
        
        // __m128i temp_1 = _mm_set_epi32(0xb7875539 ,0x79a122bd ,0x7af81a11, 0xa28b6d0);
        // __m128i temp_2 = _mm_set_epi32(0x70bd189c ,0x3a23a5e6 ,0x1d63e341, 0x33480528);
        // generate_key_x = _mm_set_epi32(0xcd1d8297 ,0x8a7729e5 ,0xbc06411c, 0xcbcb5d32);
        // generate_key_y = _mm_set_epi32(0x33b692cf ,0x73482e9d ,0x95a937fb, 0x2319ac9b);
        
        // 0xb7875539 ,0x79a122bd ,0x7af81a11, 0xa28b6d0 
        // 0x70bd189c ,0x3a23a5e6 ,0x1d63e341, 0x33480528
        // 0xcd1d8297 ,0x8a7729e5 ,0xbc06411c, 0xcbcb5d32
        // 0x33b692cf ,0x73482e9d ,0x95a937fb, 0x2319ac9b

        /*
         * XOR input blocks with the generated masks and
         * reinterpret them as 32-bit word vectors.
         */
        __m128i X = _mm_add_epi32(ptr[i], generate_key_x); 
        __m128i Y = _mm_add_epi32(ptr[i+1], generate_key_y);
        // __m128i X = _mm_add_epi32(temp_1, generate_key_x); 
        // __m128i Y = _mm_add_epi32(temp_2, generate_key_y);

        /*
         * Update the internal hash state using polynomial
         * multiplications and integer accumulation.
         */
        update_function(X, Y, output);

        
        

    }

    
    // printf("-------------------------------------------------\n");
    final_output[0]= _mm_castpd_si128(_mm_shuffle_pd(_mm_castsi128_pd(output[0]), _mm_castsi128_pd(output[1]), 0x00));
    final_output[1]= _mm_castpd_si128(_mm_shuffle_pd(_mm_castsi128_pd(output[0]), _mm_castsi128_pd(output[1]), 0x03));
    final_output[2]= _mm_castpd_si128(_mm_shuffle_pd(_mm_castsi128_pd(output[2]), _mm_castsi128_pd(output[3]), 0x00));
    final_output[3]= _mm_castpd_si128(_mm_shuffle_pd(_mm_castsi128_pd(output[2]), _mm_castsi128_pd(output[3]), 0x03));
    // imprimiArreglo_64bits(2, (unsigned long int*)&final_output[0]);
    // imprimiArreglo_64bits(2, (unsigned long int*)&final_output[1]);
    // imprimiArreglo_64bits(2, (unsigned long int*)&final_output[2]);
    // imprimiArreglo_64bits(2, (unsigned long int*)&final_output[3]);
    // exit(1);
    /*
     * Final tag generation.
     * Currently, the tag is obtained by reinterpreting the
     * reduced output values as a byte array.
     */
    memcpy(tag, final_output, 64);  

}

void generate_keys(__m128i * roundKeys, uint64_t length, __m128i * obtained_keys){
    
    int64_t i = 0;
     /*
     * Define a constant counter increment (used as domain separator /
     * block index for key generation).
     */
    uint32_t constant = 1;

    /*
     * Vectorized version of the constant and the running index.
     */
    __m128i const_vec = _mm_set1_epi32(constant);
    __m128i index     = _mm_set1_epi32(constant);

    /*
     * Compute the number of 128-bit input blocks.
     * Only full blocks are processed.
     */
    int64_t size = length / 16;  // truncation intentional

    
    /*
     * Main processing loop:
     * Blocks are processed in pairs (X, Y).
     */
    for (i = 0; i < size - 1; i = i + 2) {

  
        /*
        * Extra inside index for a better pipeline for the processor.
        */
        __m128i idx0 = index;
        index = _mm_add_epi32(index, const_vec);
        __m128i idx1 = index;
        index = _mm_add_epi32(index, const_vec);


        
        /*
         * Generate a pseudo-random mask for block X
         * using AES with roundKeys_1 and the current index.
         */
        __m128i generate_key_x =
            AES_Encrypt_rounds(idx0,
                               roundKeys, 8);

       
        
        /*
         * Generate a pseudo-random mask for block Y.
         */
        __m128i generate_key_y =
            AES_Encrypt_rounds(idx1,
                               roundKeys, 8);
        
        
        /*
         * Save the keys in memory.
         */
        obtained_keys[i] = generate_key_x;
        obtained_keys[i+1] = generate_key_y;
    }
}


static inline __m128i AES_Encrypt_rounds(__m128i tmp, __m128i *key, int rounds){
	int j;


	tmp = _mm_xor_si128 (tmp,key[0]);
	for (j=1; j<rounds; j++){
        tmp = _mm_aesenc_si128 (tmp,key[j]);
    }  
	tmp = _mm_aesenc_si128 (tmp,key[j]);
   
    return tmp;
}

static inline __m128i AES_Encrypt_rounds_static_keys(__m128i tmp, __m128i key, int rounds){
	int j;
	tmp = _mm_xor_si128 (tmp,key);
	for (j=1; j<rounds; j++)  tmp = _mm_aesenc_si128 (tmp,key);
	tmp = _mm_aesenc_si128 (tmp,key);
    return tmp;
}

static inline void AES_encrypt(__m128i tmp, __m128i *out,__m128i *key, int rounds){
	int j;
	tmp = _mm_xor_si128 (tmp,key[0]);
	for (j=1; j<rounds; j++)  tmp = _mm_aesenc_si128 (tmp,key[j]);
	tmp = _mm_aesenclast_si128 (tmp,key[j]);
	_mm_store_si128 ((__m128i*)out,tmp);
}


void AES_128_Key_Expansion(const unsigned char *userkey, void *key)
{
    __m128i x0,x1,x2;
    __m128i *kp = (__m128i *)key;
    kp[0] = x0 = _mm_loadu_si128((__m128i*)userkey);
    x2 = _mm_setzero_si128();
    EXPAND_ASSIST(x0,x1,x2,x0,255,1);   kp[1]  = x0;
    EXPAND_ASSIST(x0,x1,x2,x0,255,2);   kp[2]  = x0;
    EXPAND_ASSIST(x0,x1,x2,x0,255,4);   kp[3]  = x0;
    EXPAND_ASSIST(x0,x1,x2,x0,255,8);   kp[4]  = x0;
    EXPAND_ASSIST(x0,x1,x2,x0,255,16);  kp[5]  = x0;
    EXPAND_ASSIST(x0,x1,x2,x0,255,32);  kp[6]  = x0;
    EXPAND_ASSIST(x0,x1,x2,x0,255,64);  kp[7]  = x0;
    EXPAND_ASSIST(x0,x1,x2,x0,255,128); kp[8]  = x0;
    EXPAND_ASSIST(x0,x1,x2,x0,255,27);  kp[9]  = x0;
    EXPAND_ASSIST(x0,x1,x2,x0,255,54);  kp[10] = x0;
}