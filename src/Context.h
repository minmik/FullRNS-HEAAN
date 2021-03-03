/*
* Copyright (c) by CryptoLab inc.
* This program is licensed under a
* Creative Commons Attribution-NonCommercial 3.0 Unported License.
* You should have received a copy of the license along with this
* work.  If not, see <http://creativecommons.org/licenses/by-nc/3.0/>.
*/

#ifndef HEAANNTT_CONTEXT_H_
#define HEAANNTT_CONTEXT_H_

#include <complex>
#include <chrono>
#include <map>

#include "Common.h"
#include "Numb.h"


#define Q0_BIT_SIZE 61 ///< /// Q0 = q0 = 2^61 (fixed)
#define M_PI       3.14159265358979323846

using namespace std;

static string LOGARITHM = "Logarithm"; ///< log(x)
static string EXPONENT  = "Exponent"; ///< exp(x)
static string SIGMOID   = "Sigmoid"; ///< sigmoid(x) = exp(x) / (1 + exp(x))

class Context {
public:

	// Encryption parameters
	long logN; ///< Logarithm of Ring Dimension ex) 15
	long logNh; ///< Logarithm of Ring Dimension - 1 ex) 15 - 1 = 14
	long L; ///< Maximum Level that we want to support ex) 10
	long K; ///< The number of special modulus (usually L + 1) ex) 10 + 1 = 11
	/// q_0 q_1 q_2 ... q_(L-1)
	/// p_0 p_1 p_2 ... p_(K-1)

	long N; ///< 2^logN ex) 2^15
	long M; ///< 2 * N ex) 2^16
	long Nh; ///< N / 2 ex) 2^14

	long logp; ///< scaling factor's log ex) 55
	long p; ///< scaling factor ex) 2^55

	long h; ///< Hamming distance ex) 64 (default)
	double sigma; ///< error sampling standard deviation ex) 3.2(default)

	uint64_t* qVec; ///< vector of q
	uint64_t* pVec; ///< vector of p

	uint64_t* qrVec; // Barrett reduction
	uint64_t* prVec; // Barrett recution

	long* qTwok; // Barrett reduction
	long* pTwok; // Barrett reduction

	uint64_t* qkVec; // Montgomery reduction
	uint64_t* pkVec; // Montgomery reduction

	uint64_t* qdVec; // q * 2
	uint64_t* pdVec; // p * 2

	uint64_t* qInvVec; // q^-1 mod 2^64
	uint64_t* pInvVec; // p^-1 mod 2^64

    // for r such that
    // r^((q - 1) / prime) = 1 mod q
    // qRoot is such that
    // qRoot = r^((q - 1) / M) mod q
	uint64_t* qRoots; // M-th root of unity of q
	uint64_t* pRoots; // M-th root of unity of p

    // qRoots^-1 mod q
	uint64_t* qRootsInv; // qRoots^-1 mod q
	uint64_t* pRootsInv; // pRoots^-1 mod p

    // bitreversed order
	// 1 qRoot qRoots^2 ... qRoots^(N-1)
	// all mod q
	uint64_t** qRootPows; // 1, qRoot, qRoots^2, ..., qRoots^(N-1) mod q
	uint64_t** pRootPows; // 1, pRoot, pRoots^2, ..., pRoots^(N-1) mod p

    // bitreversed order
	// 1 qRootsIns qRootsInv^2 ... qRootsInv^(N-1)
	// all mod q
	uint64_t** qRootPowsInv; // 1, qRootsInv, qRootsInv^2, ..., qRootsInv^(N-1) mod q
	uint64_t** pRootPowsInv; // 1, pRootsIns, pRootsInv^2, ..., pRootsInv^(N-1) mod p

    // N^-1 mod q
	uint64_t* NInvModq; // N^-1 mod q
	uint64_t* NInvModp; // N^-1 mod p

    // bitreversed order
	// (qRoots^{1, 2, ..., N - 1} * 2^32 mod q) * 2^32 mod q
	uint64_t** qRootScalePows; // (qRoots^{1, 2, ..., N - 1} * 2^32 mod q) * 2^32 mod q
	uint64_t** pRootScalePows; // (pRoots^{1, 2, ..., N - 1} * 2^32 mod p) * 2^32 mod p

    // bitreversed order
	// floor(qRoots^{1, 2, ..., N - 1} * 2^64 / q)
	uint64_t** qRootScalePowsOverq; // floor(qRoots^{1, 2, ..., N - 1} * 2^64 / q)
	uint64_t** pRootScalePowsOverp; // floor(pRoots^{1, 2, ..., N - 1} * 2^64 / p)

    // bitreversed order
	// (qRootsInv^{1, 2, ..., N - 1} * 2^32 mod q) * 2^32 mod q
	uint64_t** qRootScalePowsInv; // (qRootsInv^{1, 2, ..., N - 1} * 2^32 mod q) * 2^32 mod q
	uint64_t** pRootScalePowsInv; // (pRootsInv^{1, 2, ..., N - 1} * 2^32 mod p) * 2^32 mod p

    // 2^32 * N^-1 mod q
	uint64_t* NScaleInvModq; // 2^32 * N^-1 mod q
	uint64_t* NScaleInvModp; // 2^32 * N^-1 mod p

    // for basis {q_0, q_1, q_2, ... q_l} (l <= L)
	// qHatModq[l][i] = (mult of all q_0 through q_l except q_i) mod q_i
	uint64_t** qHatModq; // qHatModq[l][i] = (mult of all q_0 through q_l except q_i) mod q_i
	
	// for fixed basis {p_0, p_1, ..., p_{K-1}}
	// pHatModp[k] = (mult of all p_0 through p_{K-1} except p_k) mod p_k
	uint64_t* pHatModp; // pHatModp[k] = (mult of all p_0 through p_{K-1} except p_k) mod p_k

    // multicative inverse of qHatModq in mod q_i
	uint64_t** qHatInvModq; // [l][i] multicative inverse of qHatModq in mod q_i
	uint64_t* pHatInvModp; // [k] multicative inverse of pHatModp in mod p_k

    // for basis {q_0, q_1, q_2, ... q_l} (l <= L)
	// qHatModq[l][i][k] = (mult of all q_0 through q_l except q_i) mod p_k
	uint64_t*** qHatModp; // qHatModq[l][i][k] = (mult of all q_0 through q_l except q_i) mod p_k

    // for fixed basis {p_0, p_1, ..., p_{K-1}}
	// pHatModq[k][i] = (mult of all p_0 through p_{K-1} except p_k) mod q_i
	uint64_t** pHatModq; // pHatModq[k][i] = (mult of all p_0 through p_{K-1} except p_k) mod q_i

    // large P = mult of all p_k's
	// P mod q_i and their inverse
	uint64_t* PModq; // [i] large P mod q_i
	uint64_t* PInvModq; // [i] large P^-1 mod q_i

    // for basis {q_0, q_1, q_2, ... q_l} (l <= L)
    // large Q = mult of all q_i's = q_0 * q_1 * ... * q_l
	// QModp[l][k] = Q mod p_k and their inverse
	uint64_t** QModp; // QModp[l][k] = Q_l mod p_k 
	uint64_t** QInvModp; // QModp[l][k] = Q_l^-1 mod p_k

    // qInvModq[i][j] = q_i^-1 mod q_j
	uint64_t** qInvModq; // [i][j] q_i^-1 mod q_j

	long* rotGroup; ///< precomputed rotation group indexes: 5^i

	complex<double>* ksiPows; ///< precomputed ksi powers (M + 1 complex array): e^(2 * pi * i / M)

	map<string, double*> taylorCoeffsMap; ///< precomputed taylor coefficients

    // size = L * N ex) 10 * 2^15
	// p = 2^55
	uint64_t* p2coeff; // index N * i + (0 ~ N-1): p^2 mod q_i
	uint64_t* pccoeff; // index N * i + (0 ~ N-1): (94.2372881 * p mod q_i) -> negateandequal
	uint64_t* p2hcoeff;// index N * i + (0 ~ N-1): 0.5 * p^2 mod q_i

	Context(long logN, long logp, long L, long K, long h = 64, double sigma = 3.2);

    // disable copy or move
	Context(const Context&) = delete;
	Context& operator=(const Context&) = delete;

	void arrayBitReverse(complex<double>* vals, const long size);
	void arrayBitReverse(uint64_t* vals, const long size);

    // plain Cooley-Tukey in-place DIT radix-2 FFT
	void fft(complex<double>* vals, const long size);
	
	// plain Cooley-Tukey in-place DIT radix-2 inverse FFT, but do not divide by size
	void fftInvLazy(complex<double>* vals, const long size);

	// plain Cooley-Tukey in-place DIT radix-2 inverse FFT, and divide by size
	void fftInv(complex<double>* vals, const long size);

    // special FFT used for decoding
	void fftSpecial(complex<double>* vals, const long size);
	
	// special inverse FFT used for encoding, lazy
	void fftSpecialInvLazy(complex<double>* vals, const long size);
	
	// special inverse FFT used for encoding
	void fftSpecialInv(complex<double>* vals, const long size);

    
	// ax = encode(vals)
	void encode(uint64_t* ax, complex<double>* vals, long slots, long l);

	// not used
	void encode(uint64_t* ax, double* vals, long slots, long l);
    
	// ax = encode(val) where val is a single complex number
	void encodeSingle(uint64_t* ax, complex<double>& val, long l);

	// not used
	void encodeSingle(uint64_t* ax, double val, long l);

    // vals = decode(ax)
	void decode(uint64_t* ax, complex<double>* vals, long slots, long l);

	// not used
	void decode(uint64_t* ax, double* vals, long slots, long l);

    // 
	void decodeSingle(uint64_t* ax, complex<double>& val, long l);
	
	// not used
	void decodeSingle(uint64_t* ax, double val, long l);

	void qiNTT(uint64_t* res, uint64_t* a, long index);
	void piNTT(uint64_t* res, uint64_t* a, long index);

	void NTT(uint64_t* res, uint64_t* a, long l, long k = 0);

	void qiNTTAndEqual(uint64_t* a, long index);
	void piNTTAndEqual(uint64_t* a, long index);

	void NTTAndEqual(uint64_t* a, long l, long k = 0);

	void qiINTT(uint64_t* res, uint64_t* a, long index);
	void piINTT(uint64_t* res, uint64_t* a, long index);

	void INTT(uint64_t* res, uint64_t* a, long l, long k = 0);

	void qiINTTAndEqual(uint64_t* a, long index);
	void piINTTAndEqual(uint64_t* a, long index);

	void INTTAndEqual(uint64_t* a, long l, long k = 0);

	// compute res = -a
	void qiNegate(uint64_t* res, uint64_t* a, long index);
	void piNegate(uint64_t* res, uint64_t* a, long index);

	void negate(uint64_t* res, uint64_t* a, long l, long k = 0);

    // compute a = -a
	void qiNegateAndEqual(uint64_t* a, long index);
	void piNegateAndEqual(uint64_t* a, long index);

	void negateAndEqual(uint64_t* a, long l, long k = 0);

	void qiAddConst(uint64_t* res, uint64_t* a, uint64_t c, long index);
	void piAddConst(uint64_t* res, uint64_t* a, uint64_t c, long index);

	void addConst(uint64_t* res, uint64_t* a, uint64_t c, long l, long k = 0);

	void qiAddConstAndEqual(uint64_t* a, uint64_t c, long index);
	void piAddConstAndEqual(uint64_t* a, uint64_t c, long index);

	void addConstAndEqual(uint64_t* a, uint64_t c, long l, long k = 0);

	void qiSubConst(uint64_t* res, uint64_t* a, uint64_t c, long index);
	void piSubConst(uint64_t* res, uint64_t* a, uint64_t c, long index);

	void subConst(uint64_t* res, uint64_t* a, uint64_t c, long l, long k = 0);

	void qiSubConstAndEqual(uint64_t* a, uint64_t c, long index);
	void piSubConstAndEqual(uint64_t* a, uint64_t c, long index);

	void subConstAndEqual(uint64_t* a, uint64_t c, long l, long k = 0);

	void qiAdd(uint64_t* res, uint64_t* a, uint64_t* b, long index);
	void piAdd(uint64_t* res, uint64_t* a, uint64_t* b, long index);

	void add(uint64_t* res, uint64_t* a, uint64_t* b, long l, long k = 0);

	void qiAddAndEqual(uint64_t* a, uint64_t* b, long index);
	void piAddAndEqual(uint64_t* a, uint64_t* b, long index);

	void addAndEqual(uint64_t* a, uint64_t* b, long l, long k = 0);

	void qiSub(uint64_t* res, uint64_t* a, uint64_t* b, long index);
	void piSub(uint64_t* res, uint64_t* a, uint64_t* b, long index);

	void sub(uint64_t* res, uint64_t* a, uint64_t* b, long l, long k = 0);

	void qiSubAndEqual(uint64_t* a, uint64_t* b, long index);
	void piSubAndEqual(uint64_t* a, uint64_t* b, long index);

	void subAndEqual(uint64_t* a, uint64_t* b, long l, long k = 0);

	void qiSub2AndEqual(uint64_t* a, uint64_t* b, long index);
	void piSub2AndEqual(uint64_t* a, uint64_t* b, long index);

	void sub2AndEqual(uint64_t* a, uint64_t* b, long l, long k = 0);

	void qiMulConst(uint64_t* res, uint64_t* a, uint64_t cnst, long index);
	void piMulConst(uint64_t* res, uint64_t* a, uint64_t cnst, long index);

	void mulConst(uint64_t* res, uint64_t* a, uint64_t cnst, long l, long k = 0);

	void qiMulConstAndEqual(uint64_t* res, uint64_t cnst, long index);
	void piMulConstAndEqual(uint64_t* res, uint64_t cnst, long index);

	void mulConstAndEqual(uint64_t* res, uint64_t cnst, long l, long k = 0);

	void qiMul(uint64_t* res, uint64_t* a, uint64_t* b, long index);
	void piMul(uint64_t* res, uint64_t* a, uint64_t* b, long index);

	void mul(uint64_t* res, uint64_t* a, uint64_t* b, long l, long k = 0);

	void mulKey(uint64_t* res, uint64_t* a, uint64_t* b, long l);

	void qiMulAndEqual(uint64_t* a, uint64_t* b, long index);
	void piMulAndEqual(uint64_t* a, uint64_t* b, long index);

	void mulAndEqual(uint64_t* a, uint64_t* b, long l, long k = 0);

	void qiSquare(uint64_t* res, uint64_t* a, long index);
	void piSquare(uint64_t* res, uint64_t* a, long index);

	void square(uint64_t* res, uint64_t* a, long l, long k = 0);

	void qiSquareAndEqual(uint64_t* a, long index);
	void piSquareAndEqual(uint64_t* a, long index);

	void squareAndEqual(uint64_t* a, long l, long k = 0);

	void evalAndEqual(uint64_t* a, long l);

    // ModUp
	void raiseAndEqual(uint64_t*& a, long l);
	void raise(uint64_t* res, uint64_t* a, long l);

    // ModDown
	void backAndEqual(uint64_t*& a, long l);
	void back(uint64_t* res, uint64_t* a, long l);

    // Rescale
	void reScaleAndEqual(uint64_t*& a, long l);
	void reScale(uint64_t* res, uint64_t* a, long l); // not implemented

    // Simply Drop higher dl mod q_i's
	void modDownAndEqual(uint64_t*& a, long l, long dl);
	uint64_t* modDown(uint64_t* a, long l, long dl);

	void leftRot(uint64_t* res, uint64_t* a, long l, long rotSlots);
	void leftRotAndEqual(uint64_t* a, long l, long rotSlots);

	void conjugate(uint64_t* res, uint64_t* a, long l);
	void conjugateAndEqual(uint64_t* a, long l);

	void mulByMonomial(uint64_t* res, uint64_t* a, long l, long mdeg);
	void mulByMonomialAndEqual(uint64_t* a, long l, long mdeg);

	void sampleGauss(uint64_t* res, long l, long k = 0);
	void sampleZO(uint64_t* res, long s, long l, long k = 0);
	void sampleUniform(uint64_t* res, long l, long k = 0);
	void sampleHWT(uint64_t* res, long l, long k = 0);
	
};

#endif
