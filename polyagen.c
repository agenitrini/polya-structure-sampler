/*
 * polyagen.c -- Uniform random sampler (up to floating-point precision)
 * of unlabelled Pólya structures of a given size:
 *
 *   polya    : unlabelled non-plane rooted trees (OEIS A000081),
 *              size = number of nodes;
 *   hier     : hierarchies (OEIS A000669), size = number of leaves,
 *              no unary node;
 *   identity : identity trees (OEIS A004111), the children of each node
 *              are pairwise non-isomorphic;
 *   binary   : non-plane binary trees (Wedderburn--Etherington, OEIS
 *              A001190), each node has 0 or 2 children, size = number of
 *              leaves.
 *
 * Method: recursive method with the "Cayley-first" order (a variant of
 * Nijenhuis--Wilf's RANRUT); tables of normalised weights w_k = rho^k T_k
 * computed by recurrence in double-double (~106 bits) for k <= K0, then by
 * the 16-term asymptotic expansion for k > K0.
 *
 * Code accompanying the paper by Antoine Genitrini, "Smooth disturbances of
 * tree functions: full asymptotic expansions and uniform random sampling of Pólya
 * structures".
 *
 * Copyright (C) 2026 Antoine Genitrini
 * SPDX-License-Identifier: GPL-3.0-or-later
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 * C99, no dependency other than libm.
 */
#define _POSIX_C_SOURCE 199309L
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <inttypes.h>
#include <math.h>
#include <time.h>

#define PI_D 3.14159265358979323846

enum { FAM_POLYA = 0, FAM_HIER = 1, FAM_IDENTITY = 2, FAM_BINARY = 3, NFAM = 4 };
static const char *const FAM_NAME[NFAM] = { "polya", "hier", "identity", "binary" };

/* ------------------------------------------------------------------------
 * Asymptotic constants (taken from asymptotic_constants.json).
 * w_k = rho^k T_k ~ (pi k^3)^{-1/2} sum_{l>=0} tau_l k^{-l}.
 * They are stored as strings (25 digits) and converted to double-double
 * at start-up, which preserves their full precision
 * (on Apple Silicon, long double == double).
 * ---------------------------------------------------------------------- */
#define NTAU_ALL 18
#define NTAU 16
static const char *const RHO_STR[NFAM] = {
    "0.3383218568992076951961126", /* polya */
    "0.2808326669842003553931876", /* hierarchies */
    "0.3972130968842400414856541", /* identity */
    "0.4026975036714412909690453"  /* non-plane binary */
};
static const char *const TAU_STR[NFAM][NTAU_ALL] = {
    { /* polya */
     "0.7797450101873204427711037", "0.07828911261061096206127536",
     "0.3929402676631860184743156", "1.537879315978838091102498",
     "8.200844090435596159219525", "57.2929147349434378773993",
     "503.044505026273581824774", "5359.600933884326033403143",
     "67342.06920114653048340563", "975425.4970695924732213026",
     "15996932.49293173312960557", "292822531.3353392231145688",
     "5914523441.293932519117405", "130599192789.897133196639",
     "3128498399789.518761197287", "80783054014688.84289107057",
     "2236301680891633.731798926", "66059608696991895.64090365" },
    { /* hierarchies */
     "0.3658015862381119383635091", "0.2409833212579280363367757",
     "0.3678657493849431874959922", "0.9991064877914853514605011",
     "4.137777553476907787441839", "23.4341024892157005930631",
     "170.118881151155534539984", "1514.745295656330157729138",
     "16007.82637588106896900866", "195812.3506172274832480797",
     "2719234.685827618773611401", "42224444.65223140029342266",
     "724386196.2702191531749218", "13597749264.15692501615747",
     "277090864449.8957298989293", "6089496262810.801393932868",
     "143526925489333.1073960124", "3610881990157578.415024002" },
    { /* identity */
     "0.6425790797442694726935542", "-0.1851197977766337042457091",
     "-0.4272427290060978715187072", "-2.255455568987212074472333",
     "-16.60970953335647840704865", "-157.900369337330271646374",
     "-1840.11051735935115513059", "-25387.34869954017828442983",
     "-404610.0663959841518232298", "-7313377.058487246535116277",
     "-147794913.8517813317904795", "-3301794456.762036714846278",
     "-80802296042.28356743359157", "-2149826267241.085225832466",
     "-61790758146990.61894442814", "-1908151484770832.689862097",
     "-63010632804365562.05915499", "-2215767775919040222.499927" },
    { /* non-plane binary (Wedderburn--Etherington) */
     "0.5650168581994860035720684", "0.3612823573677108197476931",
     "0.6527419901871130860619072", "2.617594173484722181427383",
     "16.97131750832575322949689", "148.2748599119699396127002",
     "1618.318101759760714693533", "21184.89172931456469398967",
     "323713.4928060808500714317", "5658853.653758183286108931",
     "111406245.710217386860471", "2439227653.862079683139129",
     "58793037040.75287651508711", "1546919156062.517561260404",
     "44116437573068.1411113219", "1355511125440257.626875044",
     "44639257999504804.87098245", "1568473134188784187.556593" }
};

/* ------------------------------------------------------------------------
 * Utilities
 * ---------------------------------------------------------------------- */
static void *xmalloc(size_t n)
{
    void *p = malloc(n ? n : 1);
    if (!p) { fprintf(stderr, "polyagen: out of memory (%zu bytes)\n", n); exit(2); }
    return p;
}
static void *xcalloc(size_t n, size_t s)
{
    void *p = calloc(n ? n : 1, s ? s : 1);
    if (!p) { fprintf(stderr, "polyagen: out of memory\n"); exit(2); }
    return p;
}
static void *xrealloc(void *q, size_t n)
{
    void *p = realloc(q, n ? n : 1);
    if (!p) { fprintf(stderr, "polyagen: out of memory (%zu bytes)\n", n); exit(2); }
    return p;
}
static double now_sec(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec;
}

/* ------------------------------------------------------------------------
 * Double-double arithmetic (~106-bit precision), sufficient for the
 * precomputation: the recurrence only adds positive terms (or mildly
 * alternating ones for identity trees).
 * ---------------------------------------------------------------------- */
typedef struct { double hi, lo; } dd;

static inline dd dd_qts(double a, double b) { double s = a + b; return (dd){ s, b - (s - a) }; }
static inline dd dd_ts(double a, double b)
{
    double s = a + b, bb = s - a;
    return (dd){ s, (a - (s - bb)) + (b - bb) };
}
static inline dd dd_d(double a) { return (dd){ a, 0.0 }; }
static inline dd dd_add(dd a, dd b)
{
    dd s = dd_ts(a.hi, b.hi), t = dd_ts(a.lo, b.lo);
    s.lo += t.hi;
    s = dd_qts(s.hi, s.lo);
    s.lo += t.lo;
    return dd_qts(s.hi, s.lo);
}
static inline dd dd_neg(dd a) { return (dd){ -a.hi, -a.lo }; }
static inline dd dd_sub(dd a, dd b) { return dd_add(a, dd_neg(b)); }
static inline dd dd_mul(dd a, dd b)
{
    double p = a.hi * b.hi, e = fma(a.hi, b.hi, -p);
    e += a.hi * b.lo + a.lo * b.hi;
    return dd_qts(p, e);
}
static inline dd dd_mul_d(dd a, double b)
{
    double p = a.hi * b, e = fma(a.hi, b, -p);
    e += a.lo * b;
    return dd_qts(p, e);
}
static dd dd_div(dd a, dd b)
{
    double q1 = a.hi / b.hi;
    dd r = dd_sub(a, dd_mul_d(b, q1));
    double q2 = r.hi / b.hi;
    r = dd_sub(r, dd_mul_d(b, q2));
    double q3 = r.hi / b.hi;
    return dd_add(dd_qts(q1, q2), dd_d(q3));
}
static dd dd_div_d(dd a, double b) { return dd_div(a, dd_d(b)); }
static dd dd_sqrt(dd a)
{
    double x = sqrt(a.hi);
    double p = x * x, e = fma(x, x, -p);
    dd diff = dd_sub(a, (dd){ p, e });
    return dd_add(dd_d(x), dd_d(diff.hi / (2.0 * x)));
}
/* Parse a decimal "[-]ddd.ddd" into a double-double. */
static dd dd_from_str(const char *s)
{
    int neg = 0, frac = 0, infrac = 0;
    dd acc = dd_d(0.0);
    if (*s == '-') { neg = 1; s++; } else if (*s == '+') s++;
    for (; *s; s++) {
        if (*s == '.') { infrac = 1; continue; }
        if (*s < '0' || *s > '9') break;
        acc = dd_add(dd_mul_d(acc, 10.0), dd_d((double)(*s - '0')));
        if (infrac) frac++;
    }
    dd p10 = dd_d(1.0);
    for (int i = 0; i < frac; i++) p10 = dd_mul_d(p10, 10.0);
    acc = dd_div(acc, p10);
    return neg ? dd_neg(acc) : acc;
}
static const dd DD_PI = { 3.141592653589793116e+00, 1.224646799147353207e-16 };

/* ------------------------------------------------------------------------
 * Pseudo-random generator xoshiro256** (Blackman--Vigna),
 * seeded by splitmix64.
 * ---------------------------------------------------------------------- */
typedef struct { uint64_t s[4]; } rng_t;
static inline uint64_t rotl64(uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }
static uint64_t splitmix64(uint64_t *x)
{
    uint64_t z = (*x += 0x9E3779B97F4A7C15ULL);
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}
static void rng_seed(rng_t *r, uint64_t seed)
{
    uint64_t x = seed;
    for (int i = 0; i < 4; i++) r->s[i] = splitmix64(&x);
}
static inline uint64_t rng_next(rng_t *r)
{
    uint64_t *s = r->s;
    uint64_t res = rotl64(s[1] * 5, 7) * 9, t = s[1] << 17;
    s[2] ^= s[0]; s[3] ^= s[1]; s[1] ^= s[2]; s[0] ^= s[3];
    s[2] ^= t; s[3] = rotl64(s[3], 45);
    return res;
}
/* Uniform float in [0,1) with 53 random bits. */
static inline double rng_double(rng_t *r) { return (double)(rng_next(r) >> 11) * 0x1.0p-53; }

/* ------------------------------------------------------------------------
 * Tables of normalised weights.
 *
 *   a_k = w_k = rho^k T_k   (class of the elements),
 *   f_k = rho^k F_k         (multisets of elements of total size k),
 *         polya : F_k = T_{k+1}, hence f_k = w_{k+1}/rho ... we use instead
 *         g_m = rho^{m+1} F_m = w_{m+1} (see below);
 *         hier  : f_0 = 1, f_1 = rho, f_k = 2 w_k (k >= 2).
 *
 * Probability of the pair (d,j) for a multiset of size m:
 *   P(d,j) = d T_d F_{m-jd} / (m F_m)
 *          = d a_d f_{m-jd} rho^{(j-1)d} / (m f_m)
 * (the powers of rho cancel out: rho^m = rho^d rho^{m-jd} rho^{(j-1)d}).
 * For polya with g_m = rho^{m+1} F_m the same formula holds with f = g.
 *
 * Recurrences (exact computation in double-double for k <= K0): with
 *   beta_k = sum_{d | k} eps(k/d) d w_d rho^{k-d}
 * (eps = 1 for polya/hier, eps(q) = (-1)^{q-1} for identity):
 *   polya, identity : (n-1) w_n = sum_{k=1}^{n-1} beta_k w_{n-k}
 *   hier            :  n f_n    = sum_{k=1}^{n}   beta_k f_{n-k}, f_n = 2 w_n
 * This is the "Euler" (convolution) form of the recurrences in the paper,
 * in O(K0^2) instead of O(K0^2 log K0).
 *
 * binary (T = z + T^2/2 + T(z^2)/2): w_1 = rho and, for n >= 2,
 *   2 w_n = sum_{i=1}^{n-1} w_i w_{n-i} + [n even] rho^{n/2} w_{n/2},
 * and the probability of a sampler step (pointed form, cf. gen_binary) is
 *   P(d)   = d w_d w_{n-d} / (n w_n)              (Cayley step, 1 <= d < n),
 *   P(dec) = (n/2) w_{n/2} rho^{n/2} / (n w_n)    (decoration, n even).
 * ---------------------------------------------------------------------- */
typedef struct {
    int fam;
    int64_t n, K;        /* target size, K = min(n, K0) */
    double rho;
    double tau[NTAU];
    double *w;           /* w[0..n+1] */
    double *f;           /* multiset weights (points into w for polya) */
    double *fbuf;        /* own buffer for f (hier) */
    double *rp;          /* rp[e] = rho^e for e < nrp, 0 beyond */
    int64_t nrp;
    double check_reldiff;  /* relative gap recurrence / asymptotic at K */
    int asym_used;
} tables_t;

/* Exact weights (double-double) w[0..K] by recurrence. */
static dd *build_dd(int fam, int64_t K)
{
    dd rho = dd_from_str(RHO_STR[fam]);
    dd *w = xcalloc((size_t)K + 2, sizeof(dd));
    dd *rp = xcalloc((size_t)K + 2, sizeof(dd));
    dd *beta = xcalloc((size_t)K + 2, sizeof(dd));
    dd *f = (fam == FAM_HIER) ? xcalloc((size_t)K + 2, sizeof(dd)) : NULL;
    rp[0] = dd_d(1.0);
    for (int64_t e = 1; e <= K; e++) rp[e] = dd_mul(rp[e - 1], rho);
    if (K >= 1) w[1] = rho;
    if (f) { f[0] = dd_d(1.0); if (K >= 1) f[1] = rho; }
    if (fam == FAM_BINARY) {
        for (int64_t n = 2; n <= K; n++) {
            dd s = dd_d(0.0);
            for (int64_t i = 1; 2 * i < n; i++) s = dd_add(s, dd_mul(w[i], w[n - i]));
            s = dd_mul_d(s, 2.0);
            if (n % 2 == 0) {
                dd h = w[n / 2];
                s = dd_add(s, dd_mul(h, dd_add(h, rp[n / 2])));
            }
            w[n] = dd_mul_d(s, 0.5);
        }
        free(rp); free(beta);
        return w;
    }
    for (int64_t n = 1; n <= K; n++) {
        if (n >= 2) {
            dd s = dd_d(0.0);
            if (fam == FAM_HIER) {
                s = beta[n]; /* contributions of the proper divisors of n */
                for (int64_t k = 1; k < n; k++) s = dd_add(s, dd_mul(beta[k], f[n - k]));
                w[n] = dd_div_d(s, (double)n);
                f[n] = dd_mul_d(w[n], 2.0);
            } else {
                for (int64_t k = 1; k < n; k++) s = dd_add(s, dd_mul(beta[k], w[n - k]));
                w[n] = dd_div_d(s, (double)(n - 1));
            }
        }
        /* sieve: contribution of d = n to beta_{qn} */
        dd c = dd_mul_d(w[n], (double)n);
        int64_t q = 1;
        for (int64_t k = n; k <= K; k += n, q++) {
            if (rp[k - n].hi == 0.0) break;
            dd t = dd_mul(c, rp[k - n]);
            if (fam == FAM_IDENTITY && (q % 2 == 0)) beta[k] = dd_sub(beta[k], t);
            else beta[k] = dd_add(beta[k], t);
        }
    }
    free(rp); free(beta); free(f);
    return w;
}

/* Asymptotic expansion (16 terms) in double-double. */
static dd asym_dd(int fam, int64_t k)
{
    dd x = dd_div(dd_d(1.0), dd_d((double)k)), s = dd_d(0.0);
    for (int l = NTAU - 1; l >= 0; l--) s = dd_add(dd_mul(s, x), dd_from_str(TAU_STR[fam][l]));
    dd k3 = dd_mul_d(dd_mul_d(dd_d((double)k), (double)k), (double)k);
    return dd_div(s, dd_sqrt(dd_mul(DD_PI, k3)));
}
/* Same in double (fast filling of the tables for k > K0). */
static inline double asym_d(const tables_t *T, int64_t k)
{
    double kk = (double)k, x = 1.0 / kk, s = 0.0;
    for (int l = NTAU - 1; l >= 0; l--) s = s * x + T->tau[l];
    return s / sqrt(PI_D * kk * kk * kk);
}

static void tables_build(tables_t *T, int fam, int64_t n, int64_t K0)
{
    memset(T, 0, sizeof *T);
    T->fam = fam; T->n = n;
    T->K = (K0 < n) ? K0 : n;
    if (T->K < 1) T->K = 1;
    T->rho = dd_from_str(RHO_STR[fam]).hi;
    for (int l = 0; l < NTAU; l++) T->tau[l] = dd_from_str(TAU_STR[fam][l]).hi;
    int64_t K = T->K;
    dd *wdd = build_dd(fam, K);
    T->w = xmalloc(sizeof(double) * (size_t)(n + 2));
    T->w[0] = 0.0;
    for (int64_t k = 1; k <= K; k++) T->w[k] = wdd[k].hi;
    for (int64_t k = K + 1; k <= n + 1; k++) T->w[k] = asym_d(T, k);
    T->asym_used = (n > K); /* w_{n+1} is computed but never used */
    /* consistency check at k = K */
    T->check_reldiff = -1.0;
    if (K >= 2) {
        dd a = asym_dd(fam, K), d = dd_sub(a, wdd[K]);
        T->check_reldiff = fabs(d.hi / wdd[K].hi);
    }
    free(wdd);
    if (fam == FAM_POLYA) {
        T->f = T->w + 1; /* g_m = w_{m+1} */
    } else if (fam == FAM_HIER) {
        T->fbuf = xmalloc(sizeof(double) * (size_t)(n + 2));
        T->fbuf[0] = 1.0;
        T->fbuf[1] = T->rho;
        for (int64_t k = 2; k <= n + 1; k++) T->fbuf[k] = 2.0 * T->w[k];
        T->f = T->fbuf;
    } else {
        T->f = NULL;
    }
    /* powers of rho in double (zero beyond, by underflow) */
    T->nrp = 1;
    {
        int64_t cap = n + 2;
        T->rp = xmalloc(sizeof(double) * (size_t)(cap < 4096 ? cap : 4096));
        dd r = dd_d(1.0), rho = dd_from_str(RHO_STR[fam]);
        T->rp[0] = 1.0;
        while (T->nrp < cap && T->nrp < 4096) {
            r = dd_mul(r, rho);
            if (r.hi == 0.0) break;
            T->rp[T->nrp++] = r.hi;
        }
    }
}
static void tables_free(tables_t *T)
{
    free(T->w); free(T->fbuf); free(T->rp);
    memset(T, 0, sizeof *T);
}
static inline double rho_pow(const tables_t *T, int64_t e) { return e < T->nrp ? T->rp[e] : 0.0; }

/* ------------------------------------------------------------------------
 * Hash-consing: canonical identifier of an isomorphism class
 * = identifier of the sorted list of the children's identifiers.
 * Exact (no collision possible: the lists themselves are compared).
 * ---------------------------------------------------------------------- */
typedef struct {
    uint32_t *pool; size_t plen, pcap;
    size_t *off; uint32_t *len; uint64_t *hash;
    uint32_t nent, ecap;
    uint32_t *tab; size_t tmask; /* stores id+1, 0 = empty */
} hcons_t;

static void hc_init(hcons_t *H)
{
    memset(H, 0, sizeof *H);
    H->tmask = 1023;
    H->tab = xcalloc(H->tmask + 1, sizeof(uint32_t));
}
static void hc_free(hcons_t *H)
{
    free(H->pool); free(H->off); free(H->len); free(H->hash); free(H->tab);
    memset(H, 0, sizeof *H);
}
static void hc_reset(hcons_t *H)
{
    H->plen = 0; H->nent = 0;
    memset(H->tab, 0, (H->tmask + 1) * sizeof(uint32_t));
}
static inline uint64_t hc_hashfn(const uint32_t *k, uint32_t n)
{
    uint64_t h = 0x9E3779B97F4A7C15ULL ^ ((uint64_t)n * 0xD6E8FEB86659FD93ULL);
    for (uint32_t i = 0; i < n; i++) {
        h ^= k[i];
        h *= 0xFF51AFD7ED558CCDULL;
        h ^= h >> 29;
    }
    h ^= h >> 32; h *= 0xC4CEB9FE1A85EC53ULL; h ^= h >> 29;
    return h;
}
static uint32_t hc_get(hcons_t *H, const uint32_t *k, uint32_t n)
{
    if ((size_t)(H->nent + 1) * 2 > H->tmask + 1) {
        size_t nm = (H->tmask + 1) * 2 - 1;
        uint32_t *nt = xcalloc(nm + 1, sizeof(uint32_t));
        for (uint32_t e = 0; e < H->nent; e++) {
            size_t i = H->hash[e] & nm;
            while (nt[i]) i = (i + 1) & nm;
            nt[i] = e + 1;
        }
        free(H->tab); H->tab = nt; H->tmask = nm;
    }
    uint64_t h = hc_hashfn(k, n);
    size_t i = h & H->tmask;
    while (H->tab[i]) {
        uint32_t e = H->tab[i] - 1;
        if (H->hash[e] == h && H->len[e] == n &&
            (n == 0 || memcmp(H->pool + H->off[e], k, n * sizeof(uint32_t)) == 0))
            return e;
        i = (i + 1) & H->tmask;
    }
    if (H->nent == H->ecap) {
        H->ecap = H->ecap ? 2 * H->ecap : 1024;
        H->off = xrealloc(H->off, H->ecap * sizeof(size_t));
        H->len = xrealloc(H->len, H->ecap * sizeof(uint32_t));
        H->hash = xrealloc(H->hash, H->ecap * sizeof(uint64_t));
    }
    if (H->plen + n > H->pcap) {
        while (H->plen + n > H->pcap) H->pcap = H->pcap ? 2 * H->pcap : 4096;
        H->pool = xrealloc(H->pool, H->pcap * sizeof(uint32_t));
    }
    uint32_t e = H->nent++;
    H->off[e] = H->plen; H->len[e] = n; H->hash[e] = h;
    if (n) memcpy(H->pool + H->plen, k, n * sizeof(uint32_t));
    H->plen += n;
    H->tab[i] = e + 1;
    return e;
}
static int cmp_u32(const void *a, const void *b)
{
    uint32_t x = *(const uint32_t *)a, y = *(const uint32_t *)b;
    return (x > y) - (x < y);
}
static void sort_u32(uint32_t *a, uint32_t n)
{
    if (n > 24) { qsort(a, n, sizeof(uint32_t), cmp_u32); return; }
    for (uint32_t i = 1; i < n; i++) {
        uint32_t x = a[i], j = i;
        while (j > 0 && a[j - 1] > x) { a[j] = a[j - 1]; j--; }
        a[j] = x;
    }
}

/* ------------------------------------------------------------------------
 * Tree representation: array sz[] in preorder, sz[p] = number of nodes
 * of the subtree rooted at p. The children of p are p+1,
 * p+1+sz[p+1], ... This representation is translation-invariant:
 * a subtree is a contiguous block, duplicated with memcpy.
 * ---------------------------------------------------------------------- */
typedef struct { int64_t r, m, s; int32_t j, restricted; } mframe;
typedef struct { int64_t r, m, i, s, t0; int32_t state, outer; } iframe;
typedef struct { int64_t r, m, s, m2; int32_t state; } bframe;

typedef struct {
    int32_t *sz;  int64_t cap, top;  /* tree buffer */
    uint32_t *cid;                   /* canonical identifiers (identity) */
    uint32_t *kids;                  /* work buffer */
    mframe *ms; int64_t mscap;
    iframe *is; int64_t iscap;
    bframe *bs; int64_t bscap;
    hcons_t hc;
    uint64_t steps, rejections;
    rng_t rng;
} sampler_t;

static void sampler_init(sampler_t *S, int fam, int64_t n, uint64_t seed)
{
    memset(S, 0, sizeof *S);
    S->cap = (fam == FAM_HIER || fam == FAM_BINARY) ? 2 * n : n + 1;
    S->sz = xmalloc(sizeof(int32_t) * (size_t)S->cap);
    if (fam == FAM_IDENTITY) {
        S->cid = xmalloc(sizeof(uint32_t) * (size_t)S->cap);
        S->kids = xmalloc(sizeof(uint32_t) * (size_t)S->cap);
    }
    S->mscap = 64; S->ms = xmalloc(sizeof(mframe) * (size_t)S->mscap);
    S->iscap = 64; S->is = xmalloc(sizeof(iframe) * (size_t)S->iscap);
    S->bscap = 64; S->bs = xmalloc(sizeof(bframe) * (size_t)S->bscap);
    hc_init(&S->hc);
    rng_seed(&S->rng, seed);
}
static void sampler_free(sampler_t *S)
{
    free(S->sz); free(S->cid); free(S->kids); free(S->ms); free(S->is); free(S->bs);
    hc_free(&S->hc);
}

/* Choice of the pair (d,j) for a multiset of size m, "Cayley-first"
 * order: first j = 1 with d = m, m-1, ..., 1, then j >= 2
 * (d increasing, j increasing). restricted: pair (m,1) forbidden
 * (root of a hierarchy: at least two children). */
static inline void choose_pair(const tables_t *T, rng_t *rng, int64_t m, int restricted,
                               uint64_t *steps, int64_t *pd, int64_t *pj)
{
    const double *a = T->w, *f = T->f;
    double denom = restricted ? (double)m * a[m] : (double)m * f[m];
    int64_t dmax = restricted ? m - 1 : m;
    for (;;) {
        double t = rng_double(rng) * denom;
        int64_t d;
        for (d = dmax; d >= 1; d--) {
            t -= (double)d * a[d] * f[m - d];
            if (t < 0.0) { *steps += (uint64_t)(dmax - d + 1); *pd = d; *pj = 1; return; }
        }
        *steps += (uint64_t)(dmax > 0 ? dmax : 0);
        for (d = 1; 2 * d <= m; d++) {
            double rd = rho_pow(T, d);
            if (rd == 0.0) break;
            double c = (double)d * a[d], p = 1.0;
            for (int64_t j = 2; j * d <= m; j++) {
                p *= rd;
                (*steps)++;
                t -= c * f[m - j * d] * p;
                if (t < 0.0) { *pd = d; *pj = j; return; }
                if (p == 0.0) break;
            }
        }
        /* remainder due to rounding (probability ~1e-16): start again */
    }
}

/* polya / hier sampler (recursive method, explicit stack). */
static void gen_multiset_family(sampler_t *S, const tables_t *T, int64_t n)
{
    const int hier = (T->fam == FAM_HIER);
    int32_t *sz = S->sz;
    int64_t top = 0, sp = 0;
    int64_t root = top++;
    if (n == 1) { sz[root] = 1; S->top = top; return; }
    S->ms[sp++] = (mframe){ root, hier ? n : n - 1, 0, 0, hier };
    while (sp > 0) {
        mframe *F = &S->ms[sp - 1];
        if (F->j) { /* sigma finished in [s, top): j-1 copies */
            int64_t L = top - F->s;
            for (int32_t c = 1; c < F->j; c++) {
                memcpy(sz + top, sz + F->s, (size_t)L * sizeof(int32_t));
                top += L;
            }
            F->j = 0;
        }
        if (F->m == 0) { sz[F->r] = (int32_t)(top - F->r); sp--; continue; }
        int64_t d, j;
        choose_pair(T, &S->rng, F->m, F->restricted, &S->steps, &d, &j);
        F->restricted = 0;
        F->m -= j * d;
        F->s = top;
        F->j = (int32_t)j;
        int64_t pos = top++;
        if (d == 1) { sz[pos] = 1; continue; }
        if (sp == S->mscap) {
            S->mscap *= 2;
            S->ms = xrealloc(S->ms, sizeof(mframe) * (size_t)S->mscap);
        }
        S->ms[sp++] = (mframe){ pos, hier ? d : d - 1, 0, 0, hier };
    }
    S->top = top;
}

/* v_i(m) = rho^m V_i(m) = sum_j (-1)^j rho^{ji} w_{m-ji} (indices >= 1). */
static inline double v_id(const tables_t *T, int64_t i, int64_t m)
{
    if (m <= i) return T->w[m];
    double ri = rho_pow(T, i), p = 1.0, v = 0.0, sg = 1.0;
    for (int64_t k = m; k >= 1; k -= i) {
        v += sg * p * T->w[k];
        p *= ri; sg = -sg;
        if (p == 0.0) break;
    }
    return v;
}
/* Choice of the size i of the subtree sigma (identity trees), i = m-1..1. */
static inline int64_t choose_identity(const tables_t *T, rng_t *rng, int64_t m, uint64_t *steps)
{
    const double *w = T->w;
    double denom = (double)(m - 1) * w[m];
    for (;;) {
        double t = rng_double(rng) * denom;
        for (int64_t i = m - 1; i >= 1; i--) {
            (*steps)++;
            t -= (double)i * w[i] * v_id(T, i, m - i);
            if (t < 0.0) return i;
        }
    }
}
/* Finalising a node (identity): size and canonical identifier. */
static void id_finalize(sampler_t *S, int64_t r, int64_t top)
{
    S->sz[r] = (int32_t)(top - r);
    uint32_t k = 0;
    for (int64_t c = r + 1; c < top; c += S->sz[c]) S->kids[k++] = S->cid[c];
    sort_u32(S->kids, k);
    S->cid[r] = hc_get(&S->hc, S->kids, k);
}
/* identity sampler: tree(m) = sigma(i) grafted onto the root of tau(m-i),
 * tau being redrawn as long as its root has a child isomorphic to sigma.
 * The children of the root are generated one after another (sigma_1,
 * sigma_2, ...), tau sharing the root r. */
static void gen_identity(sampler_t *S, const tables_t *T, int64_t n)
{
    int32_t *sz = S->sz;
    uint32_t *cid = S->cid;
    hc_reset(&S->hc);
    uint32_t leaf = hc_get(&S->hc, NULL, 0);
    int64_t top = 0, sp = 0;
    int64_t root = top++;
    if (n == 1) { sz[root] = 1; cid[root] = leaf; S->top = top; return; }
    S->is[sp++] = (iframe){ root, n, 0, 0, 0, 0, 1 };
    while (sp > 0) {
        if (sp + 1 >= S->iscap) {
            S->iscap *= 2;
            S->is = xrealloc(S->is, sizeof(iframe) * (size_t)S->iscap);
        }
        iframe *F = &S->is[sp - 1];
        if (F->state == 0) {
            int64_t i = choose_identity(T, &S->rng, F->m, &S->steps);
            int64_t pos = top++;
            F->i = i; F->s = pos; F->state = 1;
            if (i == 1) { sz[pos] = 1; cid[pos] = leaf; }
            else S->is[sp++] = (iframe){ pos, i, 0, 0, 0, 0, 1 };
            continue;
        }
        if (F->state == 1) {
            if (F->m - F->i == 1) { /* tau reduced to one node: nothing to test */
                if (F->outer) id_finalize(S, F->r, top);
                sp--; continue;
            }
            F->t0 = top; F->state = 2;
            S->is[sp++] = (iframe){ F->r, F->m - F->i, 0, 0, 0, 0, 0 };
            continue;
        }
        /* state 2: tau complete, children in [t0, top) */
        uint32_t cs = cid[F->s];
        int bad = 0;
        for (int64_t c = F->t0; c < top; c += sz[c])
            if (cid[c] == cs) { bad = 1; break; }
        if (bad) {
            S->rejections++;
            top = F->t0;
            S->is[sp++] = (iframe){ F->r, F->m - F->i, 0, 0, 0, 0, 0 };
            continue;
        }
        if (F->outer) id_finalize(S, F->r, top);
        sp--;
    }
    S->top = top;
}

/* Choice of the step for a non-plane binary tree of size m >= 2,
 * "Cayley-first" order: Cayley steps d = m-1, m-2, ..., 1 (probability
 * d w_d w_{m-d} / (m w_m)), then, if m is even, the decoration step
 * (probability (m/2) w_{m/2} rho^{m/2} / (m w_m)). Returns d, or 0 for the
 * decoration. */
static inline int64_t choose_binary(const tables_t *T, rng_t *rng, int64_t m, uint64_t *steps)
{
    const double *w = T->w;
    double denom = (double)m * w[m];
    for (;;) {
        double t = rng_double(rng) * denom;
        for (int64_t d = m - 1; d >= 1; d--) {
            (*steps)++;
            t -= (double)d * w[d] * w[m - d];
            if (t < 0.0) return d;
        }
        if (m % 2 == 0) {
            (*steps)++;
            t -= (double)(m / 2) * w[m / 2] * rho_pow(T, m / 2);
            if (t < 0.0) return 0;
        }
        /* remainder due to rounding (probability ~1e-16): start again */
    }
}
/* binary sampler (explicit stack). Cayley step d: root whose children
 * are sigma (size d) then phi (size m-d), drawn independently;
 * decoration: root whose two children are two copies of the same sigma
 * of size m/2 (the second one by memcpy). Uniformity: a tree {A,B}
 * (A != B) is reached via d = |A| and d = |B|, with total weight
 * |A| + |B| = m; a tree {A,A} via d = m/2 and via the decoration, with
 * total weight m; hence the probability m / (m T_m) = 1/T_m. */
static void gen_binary(sampler_t *S, const tables_t *T, int64_t n)
{
    int32_t *sz = S->sz;
    int64_t top = 0, sp = 0;
    int64_t root = top++;
    if (n == 1) { sz[root] = 1; S->top = top; return; }
    S->bs[sp++] = (bframe){ root, n, 0, 0, 0 };
    while (sp > 0) {
        if (sp == S->bscap) {
            S->bscap *= 2;
            S->bs = xrealloc(S->bs, sizeof(bframe) * (size_t)S->bscap);
        }
        bframe *F = &S->bs[sp - 1];
        int64_t c; /* size of the child to generate */
        if (F->state == 0) {
            int64_t d = choose_binary(T, &S->rng, F->m, &S->steps);
            F->s = top;
            F->m2 = d ? F->m - d : 0;   /* 0: second copy of sigma */
            F->state = 1;
            c = d ? d : F->m / 2;
        } else if (F->state == 1 && F->m2 == 0) {
            int64_t L = top - F->s;     /* sigma finished in [s, top) */
            memcpy(sz + top, sz + F->s, (size_t)L * sizeof(int32_t));
            top += L;
            sz[F->r] = (int32_t)(top - F->r); sp--;
            continue;
        } else if (F->state == 1) {
            F->state = 2;
            c = F->m2;
        } else {
            sz[F->r] = (int32_t)(top - F->r); sp--;
            continue;
        }
        int64_t pos = top++;
        if (c == 1) sz[pos] = 1;
        else S->bs[sp++] = (bframe){ pos, c, 0, 0, 0 };
    }
    S->top = top;
}

static void sample(sampler_t *S, const tables_t *T, int64_t n)
{
    if (T->fam == FAM_IDENTITY) gen_identity(S, T, n);
    else if (T->fam == FAM_BINARY) gen_binary(S, T, n);
    else gen_multiset_family(S, T, n);
}

/* ------------------------------------------------------------------------
 * Canonical forms and validation.
 * ---------------------------------------------------------------------- */
/* Canonical identifiers of all nodes (reverse preorder). */
static void canon_ids(const int32_t *sz, int64_t N, hcons_t *H, uint32_t *cid, uint32_t *kids)
{
    for (int64_t p = N - 1; p >= 0; p--) {
        uint32_t k = 0;
        for (int64_t c = p + 1; c < p + sz[p]; c += sz[c]) kids[k++] = cid[c];
        sort_u32(kids, k);
        cid[p] = hc_get(H, kids, k);
    }
}
/* Returns NULL if the tree is valid, an error message otherwise. */
static const char *validate(int fam, int64_t n, const int32_t *sz, int64_t N)
{
    if (N < 1 || sz[0] != N) return "inconsistent root size";
    int64_t leaves = 0;
    for (int64_t p = 0; p < N; p++) {
        if (sz[p] < 1 || p + sz[p] > N) return "sz out of bounds";
        int64_t deg = 0, c;
        for (c = p + 1; c < p + sz[p]; c += sz[c]) {
            if (sz[c] < 1) return "sz zero";
            deg++;
        }
        if (c != p + sz[p]) return "children do not cover the subtree";
        if (deg == 0) leaves++;
        if (fam == FAM_HIER && deg == 1) return "unary node in a hierarchy";
        if (fam == FAM_BINARY && deg != 0 && deg != 2) return "node of arity other than 0 and 2 (binary)";
    }
    if (fam == FAM_BINARY) {
        if (leaves != n) return "number of leaves != n";
        if (N != 2 * n - 1) return "number of nodes != 2n-1";
    } else if (fam == FAM_HIER) { if (leaves != n) return "number of leaves != n"; }
    else if (N != n) return "number of nodes != n";
    if (fam == FAM_IDENTITY) {
        hcons_t H; hc_init(&H);
        uint32_t *cid = xmalloc(sizeof(uint32_t) * (size_t)N);
        uint32_t *kids = xmalloc(sizeof(uint32_t) * (size_t)N);
        canon_ids(sz, N, &H, cid, kids);
        const char *err = NULL;
        for (int64_t p = 0; p < N && !err; p++) {
            uint32_t k = 0;
            for (int64_t c = p + 1; c < p + sz[p]; c += sz[c]) kids[k++] = cid[c];
            sort_u32(kids, k);
            for (uint32_t t = 1; t < k; t++)
                if (kids[t] == kids[t - 1]) { err = "two isomorphic children (identity tree)"; break; }
        }
        free(cid); free(kids); hc_free(&H);
        return err;
    }
    return NULL;
}

/* ------------------------------------------------------------------------
 * Buffered output.
 * ---------------------------------------------------------------------- */
static char obuf[1 << 16];
static size_t olen = 0;
static void oflush(void) { fwrite(obuf, 1, olen, stdout); olen = 0; }
static inline void oputc(char c) { if (olen == sizeof obuf) oflush(); obuf[olen++] = c; }
static void oputi(int64_t x)
{
    char t[24]; int k = 0;
    if (x < 0) { oputc('-'); x = -x; }
    do { t[k++] = (char)('0' + x % 10); x /= 10; } while (x);
    while (k) oputc(t[--k]);
}
static void oputs(const char *t) { while (*t) oputc(*t++); }
enum { FMT_PARENT, FMT_BRACKETS, FMT_DOT, FMT_NONE };
static void output_tree(const int32_t *sz, int64_t N, int fmt, int64_t *stk, int64_t idx)
{
    int64_t sp = 0; /* stack of the ends of open subtrees (and of the nodes) */
    if (fmt == FMT_NONE) return;
    if (fmt == FMT_DOT) {
        /* Graphviz format: one directed graph per structure, edges parent -> child,
         * nodes numbered in preorder (root 0) */
        oputs("digraph tree"); oputi(idx);
        oputs(" {\n  node [shape=circle, label=\"\", width=0.15, height=0.15];\n  edge [arrowhead=none];\n");
        if (N == 1) oputs("  0;\n");
        for (int64_t i = 0; i < N; i++) {
            while (sp > 0 && stk[2 * sp - 1] <= i) sp--;
            if (sp) { oputs("  "); oputi(stk[2 * sp - 2]); oputs(" -> "); oputi(i); oputs(";\n"); }
            stk[2 * sp] = i; stk[2 * sp + 1] = i + sz[i]; sp++;
        }
        oputs("}\n");
        return;
    }
    if (fmt == FMT_PARENT) {
        /* stk[2k] = node, stk[2k+1] = end */
        for (int64_t i = 0; i < N; i++) {
            while (sp > 0 && stk[2 * sp - 1] <= i) sp--;
            if (i) oputc(' ');
            oputi(sp ? stk[2 * sp - 2] : -1);
            stk[2 * sp] = i; stk[2 * sp + 1] = i + sz[i]; sp++;
        }
    } else {
        for (int64_t i = 0; i < N; i++) {
            while (sp > 0 && stk[sp - 1] <= i) { sp--; oputc(')'); }
            oputc('(');
            stk[sp++] = i + sz[i];
        }
        while (sp > 0) { sp--; oputc(')'); }
    }
    oputc('\n');
}

/* ------------------------------------------------------------------------
 * Self-tests.
 * ---------------------------------------------------------------------- */
typedef __int128 i128;
static const int64_t OEIS[NFAM][20] = {
    /* A000081, n = 1..20 */
    { 1, 1, 2, 4, 9, 20, 48, 115, 286, 719, 1842, 4766, 12486, 32973, 87811,
      235381, 634847, 1721159, 4688676, 12826228 },
    /* A000669, n = 1..20 */
    { 1, 1, 2, 5, 12, 33, 90, 261, 766, 2312, 7068, 21965, 68954, 218751,
      699534, 2253676, 7305788, 23816743, 78023602, 256738751 },
    /* A004111, n = 1..20 */
    { 1, 1, 1, 2, 3, 6, 12, 25, 52, 113, 247, 548, 1226, 2770, 6299, 14426,
      33209, 76851, 178618, 416848 },
    /* A001190, n = 1..20 */
    { 1, 1, 1, 2, 3, 6, 11, 23, 46, 98, 207, 451, 983, 2179, 4850, 10905,
      24631, 56011, 127912, 293547 }
};
/* Exact counts T[1..N] via the "(d,j)" recurrences of the paper. */
static void exact_counts(int fam, int N, i128 *T)
{
    i128 *F = xcalloc((size_t)N + 2, sizeof(i128));
    T[0] = 0;
    for (int n = 1; n <= N; n++) {
        if (fam == FAM_POLYA) {
            if (n == 1) { T[1] = 1; continue; }
            i128 s = 0;
            for (int d = 1; d <= n - 1; d++)
                for (int j = 1; j * d <= n - 1; j++) s += (i128)d * T[d] * T[n - j * d];
            T[n] = s / (n - 1);
        } else if (fam == FAM_HIER) {
            if (n == 1) { T[1] = 1; F[0] = 1; F[1] = 1; continue; }
            i128 s = 0;
            for (int d = 1; d <= n; d++)
                for (int j = 1; j * d <= n; j++) {
                    if (d == n && j == 1) continue;
                    s += (i128)d * T[d] * F[n - j * d];
                }
            T[n] = s / n; F[n] = 2 * T[n];
        } else if (fam == FAM_BINARY) {
            /* 2 T_n = sum_{i=1}^{n-1} T_i T_{n-i} + [n even] T_{n/2} */
            if (n == 1) { T[1] = 1; continue; }
            i128 s = 0;
            for (int i = 1; i <= n - 1; i++) s += T[i] * T[n - i];
            if (n % 2 == 0) s += T[n / 2];
            T[n] = s / 2;
        } else {
            if (n == 1) { T[1] = 1; continue; }
            i128 s = 0;
            for (int i = 1; i <= n - 1; i++) {
                int m = n - i; i128 V = 0, sg = 1;
                for (int k = m; k >= 1; k -= i) { V += sg * T[k]; sg = -sg; }
                s += (i128)i * T[i] * V;
            }
            T[n] = s / (n - 1);
        }
    }
    free(F);
}

/* Q(a,x): regularised upper incomplete gamma function. */
static double gammq(double a, double x)
{
    if (x <= 0) return 1.0;
    double gln = lgamma(a);
    if (x < a + 1.0) {
        double ap = a, sum = 1.0 / a, del = sum;
        for (int i = 0; i < 100000; i++) {
            ap += 1.0; del *= x / ap; sum += del;
            if (fabs(del) < fabs(sum) * 1e-15) break;
        }
        return 1.0 - sum * exp(-x + a * log(x) - gln);
    } else {
        double b = x + 1.0 - a, c = 1.0 / 1e-300, d = 1.0 / b, h = d;
        for (int i = 1; i < 100000; i++) {
            double an = -i * (i - a);
            b += 2.0;
            d = an * d + b; if (fabs(d) < 1e-300) d = 1e-300;
            c = b + an / c; if (fabs(c) < 1e-300) c = 1e-300;
            d = 1.0 / d;
            double del = d * c; h *= del;
            if (fabs(del - 1.0) < 1e-15) break;
        }
        return exp(-x + a * log(x) - gln) * h;
    }
}

static int nfail = 0;
#define CHECK(cond, ...) do { if (cond) { printf("  [OK]    "); } else { printf("  [FAIL]  "); nfail++; } \
                              printf(__VA_ARGS__); printf("\n"); } while (0)

static void test_tables(void)
{
    printf("(a) Tables: exact recurrences vs OEIS, double-double tables vs exact integers\n");
    for (int fam = 0; fam < NFAM; fam++) {
        const int N = 30;
        i128 T[32];
        exact_counts(fam, N, T);
        int ok = 1;
        for (int n = 1; n <= 20; n++) if (T[n] != OEIS[fam][n - 1]) ok = 0;
        CHECK(ok, "%-8s : T_1..T_20 = OEIS (T_30 = %.0f)", FAM_NAME[fam], (double)T[30]);
        if (fam == FAM_BINARY) {
            /* pointed form (Cayley/decoration reading):
             * n T_n = sum_{d=1}^{n-1} d T_d T_{n-d} + [n even] (n/2) T_{n/2} */
            int okp = 1;
            for (int n = 2; n <= N; n++) {
                i128 s = 0;
                for (int d = 1; d <= n - 1; d++) s += (i128)d * T[d] * T[n - d];
                if (n % 2 == 0) s += (i128)(n / 2) * T[n / 2];
                if (s != (i128)n * T[n]) okp = 0;
            }
            CHECK(okp, "%-8s : n T_n = sum_d d T_d T_{n-d} + [n even] (n/2) T_{n/2}, n = 2..30", FAM_NAME[fam]);
        }
        dd *w = build_dd(fam, N);
        dd rho = dd_from_str(RHO_STR[fam]), rp = dd_d(1.0);
        double maxrel = 0;
        for (int n = 1; n <= N; n++) {
            rp = dd_mul(rp, rho);
            double hi = (double)T[n];
            double lo = (double)(T[n] - (i128)hi);
            dd ex = dd_mul((dd){ hi, lo }, rp);
            double rel = fabs(dd_sub(w[n], ex).hi / ex.hi);
            if (rel > maxrel) maxrel = rel;
        }
        free(w);
        CHECK(maxrel < 1e-28, "%-8s : max_n<=30 |w_n - rho^n T_n|/w_n = %.2e (double-double)", FAM_NAME[fam], maxrel);
    }
    printf("(a') Recurrence (double-double) vs 16-term asymptotic expansion\n");
    for (int fam = 0; fam < NFAM; fam++) {
        const int64_t K = 3000;
        dd *w = build_dd(fam, K);
        int64_t ks[] = { 60, 100, 120, 150, 200, 500, 1000, 2000, 3000 };
        printf("  %-8s :", FAM_NAME[fam]);
        double at120 = 0, atK = 0;
        for (size_t t = 0; t < sizeof ks / sizeof ks[0]; t++) {
            dd a = asym_dd(fam, ks[t]);
            double rel = fabs(dd_sub(a, w[ks[t]]).hi / w[ks[t]].hi);
            printf(" k=%" PRId64 ":%.1e", ks[t], rel);
            if (ks[t] == 120) at120 = rel;
            if (ks[t] == K) atK = rel;
        }
        printf("\n");
        free(w);
        CHECK(atK < 1e-20 && at120 < 2e-16, "%-8s : relative gap %.1e at k=120, %.1e at k=3000",
              FAM_NAME[fam], at120, atK);
    }
}

/* (d) sum of probabilities = 1 (mixed recurrence + asymptotic tables). */
static void test_probsum(void)
{
    printf("(d) Sum of the pair probabilities (tables with K0=3000, n=2e5)\n");
    for (int fam = 0; fam < NFAM; fam++) {
        tables_t T;
        int64_t n = 200000;
        tables_build(&T, fam, n, 3000);
        int64_t ms[] = { 2, 3, 10, 100, 1000, 2999, 3000, 3001, 3002, 5000, 20000, 199999 };
        double maxerr = 0;
        for (size_t t = 0; t < sizeof ms / sizeof ms[0]; t++) {
            int64_t m = ms[t];
            for (int restr = 0; restr <= (fam == FAM_HIER); restr++) {
                double s = 0;
                if (fam == FAM_IDENTITY) {
                    for (int64_t i = m - 1; i >= 1; i--) s += (double)i * T.w[i] * v_id(&T, i, m - i);
                    s /= (double)(m - 1) * T.w[m];
                } else if (fam == FAM_BINARY) {
                    for (int64_t d = 1; d <= m - 1; d++) s += (double)d * T.w[d] * T.w[m - d];
                    if (m % 2 == 0) s += (double)(m / 2) * T.w[m / 2] * rho_pow(&T, m / 2);
                    s /= (double)m * T.w[m];
                } else {
                    const double *a = T.w, *f = T.f;
                    for (int64_t d = 1; d <= m; d++) {
                        double rd = rho_pow(&T, d), p = 1.0;
                        for (int64_t j = 1; j * d <= m; j++) {
                            if (restr && d == m && j == 1) continue;
                            s += (double)d * a[d] * f[m - j * d] * p;
                            p *= rd;
                            if (p == 0.0) break;
                        }
                    }
                    s /= restr ? (double)m * a[m] : (double)m * f[m];
                }
                if (fabs(s - 1) > maxerr) maxerr = fabs(s - 1);
            }
        }
        tables_free(&T);
        CHECK(maxerr < 1e-12, "%-8s : max |sum - 1| = %.2e", FAM_NAME[fam], maxerr);
    }
}

/* (b) uniformity: all structures appear, chi-squared test. */
static void test_uniform(void)
{
    printf("(b) Uniformity (2e5 samples per size, chi-squared test)\n");
    const int64_t NS = 200000;
    for (int fam = 0; fam < NFAM; fam++) {
        int nlo = 5, nhi = (fam == FAM_IDENTITY || fam == FAM_BINARY) ? 10 : 9;
        i128 Tex[16];
        exact_counts(fam, 12, Tex);
        for (int n = nlo; n <= nhi; n++) {
            tables_t T; tables_build(&T, fam, n, 3000);
            sampler_t S; sampler_init(&S, fam, n, 12345 + 1000 * fam + n);
            hcons_t H; hc_init(&H);
            uint32_t *cid = xmalloc(sizeof(uint32_t) * (size_t)S.cap);
            uint32_t *kids = xmalloc(sizeof(uint32_t) * (size_t)S.cap);
            int64_t *cnt = NULL; size_t ccap = 0;
            int64_t invalid = 0;
            for (int64_t s = 0; s < NS; s++) {
                sample(&S, &T, n);
                if (validate(fam, n, S.sz, S.top)) invalid++;
                canon_ids(S.sz, S.top, &H, cid, kids);
                uint32_t id = cid[0];
                if (id >= ccap) {
                    size_t nc = ccap ? ccap : 256;
                    while (nc <= id) nc *= 2;
                    cnt = xrealloc(cnt, nc * sizeof(int64_t));
                    memset(cnt + ccap, 0, (nc - ccap) * sizeof(int64_t));
                    ccap = nc;
                }
                cnt[id]++;
            }
            int64_t ncls = (int64_t)Tex[n], seen = 0;
            double E = (double)NS / (double)ncls, chi2 = 0;
            for (size_t id = 0; id < ccap; id++)
                if (cnt[id]) { seen++; chi2 += ((double)cnt[id] - E) * ((double)cnt[id] - E) / E; }
            chi2 += (double)(ncls - seen) * E; /* missing classes */
            double pval = (ncls > 1) ? gammq(0.5 * (double)(ncls - 1), 0.5 * chi2) : 1.0;
            CHECK(seen == ncls && invalid == 0 && pval > 1e-4,
                  "%-8s n=%2d : %4" PRId64 "/%4" PRId64 " classes seen, chi2=%8.2f (df=%" PRId64 "), p=%.3f, invalid=%" PRId64,
                  FAM_NAME[fam], n, seen, ncls, chi2, ncls - 1, pval, invalid);
            free(cnt); free(cid); free(kids); hc_free(&H);
            sampler_free(&S); tables_free(&T);
        }
    }
}

/* (c) valid structures for various sizes (including k > K0). */
static void test_valid(void)
{
    printf("(c) Validity of the structures (size, class constraints)\n");
    struct { int64_t n, K0, count; } cases[] = {
        { 1, 3000, 5 }, { 2, 3000, 20 }, { 3, 3000, 50 }, { 17, 3000, 1000 },
        { 100, 3000, 500 }, { 1000, 3000, 100 }, { 20000, 3000, 10 }, { 50000, 300, 5 }
    };
    for (int fam = 0; fam < NFAM; fam++) {
        for (size_t c = 0; c < sizeof cases / sizeof cases[0]; c++) {
            int64_t n = cases[c].n;
            tables_t T; tables_build(&T, fam, n, cases[c].K0);
            sampler_t S; sampler_init(&S, fam, n, 777 + c);
            const char *err = NULL;
            for (int64_t s = 0; s < cases[c].count && !err; s++) {
                sample(&S, &T, n);
                err = validate(fam, n, S.sz, S.top);
            }
            CHECK(!err, "%-8s n=%6" PRId64 " K0=%4" PRId64 " x%4" PRId64 " : %s", FAM_NAME[fam], n,
                  cases[c].K0, cases[c].count, err ? err : "valid");
            sampler_free(&S); tables_free(&T);
        }
    }
}

/* (e) statistical check at moderate size: P(root of degree 1).
 *   polya : T_{n-1}/T_n ; identity : T_{n-1}/T_n ;
 *   hier  : P(binary root) = number of 2-element multisets / T_n ;
 *   binary : P(root split into {1, n-1}) = T_{n-1}/T_n (n >= 3). */
static void test_rootstat(void)
{
    printf("(e) Root statistic at n=25 (2e5 samples, deviation in standard deviations)\n");
    const int n = 25; const int64_t NS = 200000;
    for (int fam = 0; fam < NFAM; fam++) {
        i128 Tex[32]; exact_counts(fam, 30, Tex);
        double pex;
        if (fam == FAM_HIER) {
            i128 c = 0;
            for (int a = 1; 2 * a <= n; a++) {
                int b = n - a;
                if (a < b) c += Tex[a] * Tex[b];
                else c += Tex[a] * (Tex[a] + 1) / 2;
            }
            pex = (double)c / (double)Tex[n];
        } else pex = (double)Tex[n - 1] / (double)Tex[n];
        tables_t T; tables_build(&T, fam, n, 3000);
        sampler_t S; sampler_init(&S, fam, n, 4242 + fam);
        int64_t hit = 0;
        for (int64_t s = 0; s < NS; s++) {
            sample(&S, &T, n);
            int64_t deg = 0;
            for (int64_t c = 1; c < S.top; c += S.sz[c]) deg++;
            if (fam == FAM_BINARY) {
                if (deg == 2 && (S.sz[1] == 1 || S.sz[1 + S.sz[1]] == 1)) hit++;
            } else if ((fam == FAM_HIER && deg == 2) || (fam != FAM_HIER && deg == 1)) hit++;
        }
        double p = (double)hit / NS, sd = sqrt(pex * (1 - pex) / NS), z = (p - pex) / sd;
        CHECK(fabs(z) < 4.5, "%-8s : P(%s) exact=%.5f, observed=%.5f, z=%+.2f", FAM_NAME[fam],
              fam == FAM_HIER ? "binary root" : fam == FAM_BINARY ? "root {1,n-1}" : "unary root",
              pex, p, z);
        sampler_free(&S); tables_free(&T);
    }
}

static int selftest(void)
{
    printf("=== polyagen: self-tests ===\n");
    test_tables();
    test_probsum();
    test_valid();
    test_uniform();
    test_rootstat();
    printf("=== %s (%d failure(s)) ===\n", nfail ? "FAILURE" : "SUCCESS", nfail);
    return nfail ? 1 : 0;
}

/* ------------------------------------------------------------------------
 * Main program.
 * ---------------------------------------------------------------------- */
static void usage(void)
{
    fprintf(stderr,
        "usage: polyagen --family polya|hier|identity|binary -n SIZE [--count K]\n"
        "                [--seed S] [--format parent|brackets|dot|none] [--stats] [--kexact K0]\n"
        "       polyagen --selftest\n"
        "  -n SIZE       size (nodes for polya/identity, leaves for hier/binary);\n"
        "                notation 1e6 accepted\n"
        "  --count K     number of structures (default 1)\n"
        "  --seed S      seed (default 1)\n"
        "  --format F    parent (default): parent array in preorder, root = -1;\n"
        "                brackets: Dyck word \"(()())\"; dot: Graphviz graph\n"
        "                (one digraph per structure); none: no output\n"
        "  --stats       per-structure statistics and means (on stderr)\n"
        "  --kexact K0   exact tables up to K0 (default 3000), asymptotic beyond\n");
    exit(1);
}
static int64_t parse_size(const char *s)
{
    char *end;
    double x = strtod(s, &end);
    if (*end || x < 0 || x > 1e9) { fprintf(stderr, "polyagen: invalid value '%s'\n", s); exit(1); }
    return (int64_t)llround(x);
}

int main(int argc, char **argv)
{
    int fam = -1, fmt = FMT_PARENT, stats = 0;
    int64_t n = -1, count = 1, K0 = 3000;
    uint64_t seed = 1;
    for (int i = 1; i < argc; i++) {
        const char *a = argv[i];
#define NEXT() (i + 1 < argc ? argv[++i] : (usage(), (char *)NULL))
        if (!strcmp(a, "--selftest")) return selftest();
        else if (!strcmp(a, "--family")) {
            const char *v = NEXT();
            if (!strcmp(v, "polya")) fam = FAM_POLYA;
            else if (!strcmp(v, "hier") || !strcmp(v, "hierarchies")) fam = FAM_HIER;
            else if (!strcmp(v, "identity")) fam = FAM_IDENTITY;
            else if (!strcmp(v, "binary")) fam = FAM_BINARY;
            else usage();
        } else if (!strcmp(a, "-n")) n = parse_size(NEXT());
        else if (!strcmp(a, "--count")) count = parse_size(NEXT());
        else if (!strcmp(a, "--seed")) seed = strtoull(NEXT(), NULL, 10);
        else if (!strcmp(a, "--kexact")) K0 = parse_size(NEXT());
        else if (!strcmp(a, "--stats")) stats = 1;
        else if (!strcmp(a, "--format")) {
            const char *v = NEXT();
            if (!strcmp(v, "parent")) fmt = FMT_PARENT;
            else if (!strcmp(v, "brackets")) fmt = FMT_BRACKETS;
            else if (!strcmp(v, "dot")) fmt = FMT_DOT;
            else if (!strcmp(v, "none")) fmt = FMT_NONE;
            else usage();
        } else usage();
#undef NEXT
    }
    if (fam < 0 || n < 1 || K0 < 1) usage();

    double t0 = now_sec();
    tables_t T;
    tables_build(&T, fam, n, K0);
    double tpre = now_sec() - t0;
    if (T.asym_used && (T.check_reldiff < 0 || T.check_reldiff > 1e-15))
        fprintf(stderr, "polyagen: WARNING: at k=K0=%" PRId64 ", the recurrence and the asymptotic "
                "expansion differ by %.2e (relative); increase --kexact.\n", T.K, T.check_reldiff);
    if (stats)
        fprintf(stderr, "# family=%s n=%" PRId64 " K0=%" PRId64 " precomputation=%.4fs gap(K0)=%.2e%s\n",
                FAM_NAME[fam], n, T.K, tpre, T.check_reldiff, T.asym_used ? "" : " (asymptotics not used)");

    sampler_t S;
    sampler_init(&S, fam, n, seed);
    int64_t *ostk = (fmt != FMT_NONE) ? xmalloc(sizeof(int64_t) * 2 * (size_t)(S.cap + 1)) : NULL;
    double nlogn = (double)n * log((double)n);
    double sum_t = 0, sum_ratio = 0, sum_rej = 0, sum_nodes = 0;
    int64_t nbad = 0;
    for (int64_t c = 0; c < count; c++) {
        S.steps = 0; S.rejections = 0;
        double ts = now_sec();
        sample(&S, &T, n);
        double dt = now_sec() - ts;
        output_tree(S.sz, S.top, fmt, ostk, c);
        if (stats) {
            const char *err = validate(fam, n, S.sz, S.top);
            if (err) nbad++;
            double ratio = nlogn > 0 ? (double)S.steps / nlogn : 0.0;
            sum_t += dt; sum_ratio += ratio; sum_rej += (double)S.rejections; sum_nodes += (double)S.top;
            fprintf(stderr, "# tree %" PRId64 ": nodes=%" PRId64 " %s steps=%" PRIu64 " steps/(n ln n)=%.4f",
                    c, S.top, err ? err : "size_ok", S.steps, ratio);
            if (fam == FAM_IDENTITY) fprintf(stderr, " rejections=%" PRIu64, S.rejections);
            fprintf(stderr, " time=%.6fs\n", dt);
        }
    }
    oflush();
    if (stats && count > 0)
        fprintf(stderr, "# MEAN family=%s n=%" PRId64 " count=%" PRId64 " nodes=%.1f steps/(n ln n)=%.4f"
                " time=%.6fs invalid=%" PRId64 " precomputation=%.4fs\n",
                FAM_NAME[fam], n, count, sum_nodes / count, sum_ratio / count, sum_t / count, nbad, tpre);
    if (stats && fam == FAM_IDENTITY && count > 0)
        fprintf(stderr, "# MEAN rejections=%.2f\n", sum_rej / count);
    free(ostk);
    sampler_free(&S);
    tables_free(&T);
    return 0;
}
