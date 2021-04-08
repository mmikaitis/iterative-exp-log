#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <stdbool.h>
#include <stdlib.h>

// Number of bits in the fractional part of the internal representation (<= 60).
#define FRACT_BITS 40

#define CALC_LOG 0
#define CALC_EXP 1

// Carry-save number representation.
typedef struct {
  int64_t s;            // Intermediate sum.
  int64_t c;            // Intermediate carries.
} cs_t;

double exp_log_iterative(uint64_t x, unsigned int n, bool exp_not_log,
                         bool print);

// Functions that simulate various hardware arithmetic units.
cs_t CSA_64(uint64_t x, uint64_t y, uint64_t z);
cs_t CSA_64_4to2(uint64_t x, uint64_t y, uint64_t z, uint64_t o,
                 int cin0, int cin1);
cs_t right_shifter(cs_t x, int shift_by);
uint64_t cs_to_binary(cs_t x);

// Conversion between binary64 and fixed-point number formats.
double fixed_to_double(uint64_t x);
uint64_t double_to_fixed(double x);

// Tables of constants used in the iterative algorithm.
static int64_t log_table [64];
static int64_t log_table_neg [64];


int main() {
  
  // Calculate exponential of an arbitrary value.
  double x = 0.5;
  exp_log_iterative(double_to_fixed(x), 16, CALC_EXP, 1);

  // Calculate logarithm of an arbitrary value.
  x = 1.5;
  exp_log_iterative(double_to_fixed(x), 16, CALC_LOG, 1);

  return (0);
}


/*
  Iterative exp(x)/ln(x) algorithm in carry-save representation.

  Algorithm from p. 139, Chapter 8 of Elementary Functions: Algorithms
  and Implementation (3rd edition) by J.-M. Muller.

  Input x is the number to be exponentiated (in the interval [-1.2, ~0.86) and
  in a chosen format); n is the number of iterations.

  For logarithm, x must be in the interval [~0.4, ~3.4].

*/
double exp_log_iterative(uint64_t x, unsigned int n, bool exp_not_log,
                         bool print) {

  cs_t E = {.s = 0, .c = 0};
  cs_t L = {.s = 0, .c = 0};

  // Initialize E_1 and L_1.
  if (exp_not_log) {
    E.s = 1152921504606846976 >> (60 - FRACT_BITS); // 1.0
    L.s = x;
  }
  else {
    E.s = x;
    L.s = 0;
  }

  if (print) {
    if (exp_not_log)
      printf("Exponential of %f", fixed_to_double(x));
    else
      printf("Logarithm of %f", fixed_to_double(x));
    printf("\n");
    printf("========================================================\n");

    printf(" i            E_n                    L_n            d \n");
    printf(" 0 %23.20f %23.20f \n", fixed_to_double(cs_to_binary(E)),
           fixed_to_double(cs_to_binary(L)));
  }

  // Start from iteration 2 as ln(1 - 2^0) is not defined.
  for (int i = 2; i <= n; i++) {
    // Calculate L* or lambda*.
    cs_t tmp = {.s = 0, .c = 0};
    if (exp_not_log) {
      tmp.s = L.s << (i-1);
      tmp.c = L.c << (i-1);
    } else {
      // Get En - 1.
      cs_t tmp0 = CSA_64 (E.s, E.c, 0xffffffffffffffff << FRACT_BITS);
      tmp.s = tmp0.s << (i-1);
      tmp.c = tmp0.c << (i-1);
    } 

    // Leave 3 integer and 1 fractional bits.
    cs_t EL_star_cs = right_shifter(tmp, FRACT_BITS-1);
    EL_star_cs.s = EL_star_cs.s & 0xF;
    EL_star_cs.c = EL_star_cs.c & 0xF;

    // Convert from carry-save to non-redundant
    // representation
    uint64_t EL_star_b = cs_to_binary(EL_star_cs);
    EL_star_b &= 0xF;

    // Choose d_n.
    int d = 0;
    if (exp_not_log)
      if (EL_star_b == 0x0 || EL_star_b == 0x1 || EL_star_b == 0x2
          || EL_star_b == 0x3)
        d = 1;
      else if (EL_star_b == 0xA || EL_star_b == 0xB || EL_star_b == 0xC
               || EL_star_b == 0xD)
        d = -1;
      else if (EL_star_b == 0xE || EL_star_b == 0xF)
        d = 0;
      else
        printf("WARNING: Impossible case at iteration %d. \n", i);
    else
      if (EL_star_b == 0x0 || EL_star_b == 0xF || (EL_star_b==0x1 && i==1))
        d = 0;
      else if (EL_star_b == 0xA || EL_star_b == 0xB || EL_star_b == 0xC
               || EL_star_b == 0xD || EL_star_b == 0xE)
        d = 1;
      else if (EL_star_b == 0x1 || EL_star_b == 0x2 || EL_star_b == 0x3
               || EL_star_b == 0x4 || EL_star_b == 0x5 || EL_star_b == 0x6
               || EL_star_b == 0x7 || EL_star_b == 0x8 || EL_star_b == 0x9)
        d = -1;
      else
        printf("WARNING: Impossible case at iteration %d. \n", i);

    // Calculate E_n+1 and L_n+1.
    if (d == -1) {
      L  = CSA_64(L.s, L.c, log_table_neg[i-1] >> (60-FRACT_BITS));
      cs_t temp0 = right_shifter(E, i-1);
      E = CSA_64_4to2(E.s, E.c, ~temp0.s, ~temp0.c, 1, 1);
    }
    else if (d == 1) {
      L = CSA_64(L.s, L.c, (log_table[i-1]) >> (60-FRACT_BITS));
      cs_t temp0 = right_shifter(E, i-1);
      cs_t temp1 = CSA_64 (E.s, E.c, temp0.s);
      E = CSA_64(temp1.s, temp1.c, temp0.c);
    }
        
    if (print)
      printf("%2d %23.20f %23.20f %2d \n", i-1,
             fixed_to_double(cs_to_binary(E)),
             fixed_to_double(cs_to_binary(L)),
             d);
  }

  // Print out some information about the calculation.
  if (print) {
    double approx, ref;
    if (exp_not_log) {
      approx = fixed_to_double(cs_to_binary(E));
      ref = exp(fixed_to_double(x));
    }
    else {
      approx = fixed_to_double(cs_to_binary(L));
      ref = log(fixed_to_double(x));
    }
 
    printf("\n");
    printf("Approximation:        %33.30f \n", approx);
    printf("Double precision ref: %33.30f \n", ref);
    printf("Abs. error:           %33.30f \n", ref-approx);
    printf("Iterations performed: %3d \n", n);
    printf("Bits in the fraction: %3d \n", FRACT_BITS);
    printf("Machine epsilon:      %33.30f \n", pow(2, -FRACT_BITS));
    printf("========================================================\n");
    printf("\n");
  }

  // Use a ripple carry 64-bit adder to convert to non-redundant repr.
  if (exp_not_log)
    return fixed_to_double(cs_to_binary(E));
  else
    return fixed_to_double(cs_to_binary(L));
}

// Carry-save adder made out of 3:2 compressors (full-adders).
// Adds three 64-bit numbers and produces two 64-bit numbers - intermediate
// sum and carry.
cs_t CSA_64(uint64_t x, uint64_t y, uint64_t z) {
  cs_t result = {0, 0};
  result.s = x ^ y ^ z;
  result.c = ((x & y) | (x & z) | (z & y)) << 1;

  return result;
}

// 4:2 carry-save adder. Adds two carry-save numbers.
cs_t CSA_64_4to2(uint64_t x, uint64_t y, uint64_t z, uint64_t o,
                 int cin0, int cin1) {
  cs_t result = {0, 0};
  uint64_t majority = (((x & y) | (y & z) | (x & z)) << 1) + cin0;
  uint64_t odd_parity = ((x ^ y) ^ (z ^ o));

  result.s = odd_parity ^ majority;
  result.c = (((odd_parity & majority) + (~odd_parity & o)) << 1) + cin1;

  return result;
}

// Binary adder.
uint64_t cs_to_binary(cs_t x) {
  return x.s + x.c;
}

// Right shifter for carry-save numbers.
// Note the issue shown in Tenca et al 2006 which would require a more complex
// logic in the shifter. https://doi.org/10.1109/TC.2006.70.
// As far as we are aware, the issue does not appear in this exp/log algorithm.
cs_t right_shifter(cs_t x, int shift_by) {
  x.s = x.s >> shift_by;
  x.c = x.c >> shift_by;

  return x;
}

// Convert a fixed-point number to the binary64 representation.
double fixed_to_double(uint64_t x) {
  int negative = (x & ((uint64_t)0x1 << (FRACT_BITS+4))) != 0;
  if (negative)
    x = ~x + 1;
  uint64_t integer = x >> FRACT_BITS;
  double exponent =
    (double)(x - (integer << FRACT_BITS)) / ((uint64_t)1 << FRACT_BITS);

  double result = integer + exponent;
  if (negative)
    result = -result;
    
  return result;
}

// Convert a binary64 number to a fixed-point representation.
uint64_t double_to_fixed (double x) {
  bool negative = x < 0;
  if (negative) x = -1 * x;
  uint64_t temp = floor(x);
  uint64_t result = (uint64_t)((temp << FRACT_BITS)
                               + ((uint64_t)1 << FRACT_BITS) * (x - temp));
  if (negative)
    result = -result;
  
  return result;
}

// Two's complement log_table [n] = log (1 + 2^(-n)) in s3.60 format.
static int64_t log_table [64] = {
  17647599783384385637u,
  17979274631203908867u,
  18189477074785057738u,
  18310949479023432097u,
  18376848643508726477u,
  18411266766700229631u,
  18428868963640801336u,
  18437771876642035831u,
  18442249247335610912u,
  18444494470059998150u,
  18445618723200870878u,
  18446181261150360911u,
  18446462633086987946u,
  18446603344810431893u,
  18446673707112770223u,
  18446708889874322774u,
  18446726481657723563u,
  18446735277650083669u,
  18446739675671429099u,
  18446741874688393213u,
  18446742974198448128u,
  18446743523953868800u,
  18446743798831677440u,
  18446743936270606336u,
  18446744004990076928u,
  18446744039349813760u,
  18446744056529682560u,
  18446744065119617056u,
  18446744069414584328u,
  18446744071562067970u,
  18446744072635809792u,
  18446744073172680704u,
  18446744073441116160u,
  18446744073575333888u,
  18446744073642442752u,
  18446744073675997184u,
  18446744073692774400u,
  18446744073701163008u,
  18446744073705357312u,
  18446744073707454464u,
  18446744073708503040u,
  18446744073709027328u,
  18446744073709289472u,
  18446744073709420544u,
  18446744073709486080u,
  18446744073709518848u,
  18446744073709535232u,
  18446744073709543424u,
  18446744073709547520u,
  18446744073709549568u,
  18446744073709550592u,
  18446744073709551104u,
  18446744073709551360u,
  18446744073709551488u,
  18446744073709551552u,
  18446744073709551584u,
  18446744073709551600u,
  18446744073709551608u,
  18446744073709551612u,
  18446744073709551614u,
  18446744073709551615u,
  0,
  0,
  0
};

// log_table [n] = log (1 - 2^(-n)) * (-1) in s3.60 format.
static int64_t log_table_neg [64] = {
  0,
  799144290325165979,
  331674847819523230,
  153951214096912252,
  74407848895029353,
  36603757030154788,
  18156619410792733,
  9042567959264482,
  4512418694204213,
  2254001704453199,
  1126450020832802,
  563087437130417,
  281509342042454,
  140746078989035,
  70370891748697,
  35184908970667,
  17592320263509,
  8796126576811,
  4398054899733,
  2199025352707,
  1099512152064,
  549755944960,
  274877939712,
  137438961664,
  68719478784,
  34359738880,
  17179869312,
  8589934624,
  4294967304,
  2147483650,
  1073741825,
  536870912,
  268435456,
  134217728,
  67108864,
  33554432,
  16777216,
  8388608,
  4194304,
  2097152,
  1048576,
  524288,
  262144,
  131072,
  65536,
  32768,
  16384,
  8192,
  4096,
  2048,
  1024,
  512,
  256,
  128,
  64,
  32,
  16,
  8,
  4,
  2,
  1,
  1,
  0,
  0
};
