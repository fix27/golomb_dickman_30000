// Постоянная Голомба-Дикмана  lambda = integral_0^inf exp(-u - E1(u)) du = 0.62432998854355087099...
// Версия на C++ (FLINT/Arb, C API) для десятков тысяч знаков.  Повторяет golomb_dickman.py (v2).
//
// Сборка (FLINT >= 3.0, Arb внутри FLINT):
//     g++ -O3 -std=c++17 -pthread golomb_dickman.cpp -o golomb_dickman -lflint -lgmp -lmpfr
// FLINT 2.x + отдельный Arb: тот же файл (заголовки ищутся и как <flint/arb.h>, и как <arb.h>):
//     g++ -O3 -std=c++17 -pthread golomb_dickman.cpp -o golomb_dickman -larb -lflint -lgmp -lmpfr
// Если в вашей версии нет arb_hypgeom_expint, добавьте -DGD_USE_EI (E1(u) = -Ei(-u) через arb_hypgeom_ei).
//
// Метод (подробности и измерения - в golomb_dickman.py):
//  * подстановка Такахаси-Мори u = exp(t - e^{-t}), du = u (1 + e^{-t}) dt, шаг h = const/T,
//    T = (digits+guard) ln 10, const = 8 (измеренный показатель ошибки K ~ 9.8 не зависит от числа знаков);
//  * точность каждого узла своя: (tol + ln f)/ln 2 бит, tol = T + 30; E1 - с нужной относительной точностью
//    (tol - 2u)/ln 2 бит, а при 2u > tol отбрасывается;
//  * u < ~2^-500: F = e^gamma u exp(V(u)) вместо E1; u >= ~tol/2.2: собственный асимптотический ряд E1;
//  * чётные/нечётные узлы суммируются раздельно: Delta = 2h|S_even - S_odd| даёт оценку ошибки квадратуры ~ Delta^2;
//  * потоки (не процессы), куски узлов раздаются динамически, каждый готовый кусок пишется в чекпоинт;
//  * вся арифметика, включая сборку частичных сумм, - в шарах Arb, радиус результата проверяется.

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <memory>
#include <mutex>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

#if __has_include(<flint/arb.h>)
#include <flint/flint.h>
#include <flint/arb.h>
#include <flint/arb_hypgeom.h>
#else
#include <flint.h>
#include <arb.h>
#include <arb_hypgeom.h>
#endif

namespace fs = std::filesystem;

namespace {

constexpr double LN2 = 0.693147180559945309417232121458;
constexpr double LN10 = 2.302585092994045684017991454684;
constexpr double LOG2_10 = LN10 / LN2;
constexpr double EULER = 0.5772156649015329;

constexpr double DEFAULT_CONST = 8.0;  // h = const / T
constexpr int NODE_PAD_BITS = 96;      // запас бит на узле
constexpr int SERIES_KMAX = 400;       // макс. число членов ряда V(u)
constexpr int ASYM_KMAX = 3000;        // макс. число членов асимптотического ряда E1
constexpr double E1_SKIP_BITS = -8.0;  // нужная точность E1 хуже 2^8 -> E1 отбрасываем
// FLINT expint "висит" в окне u ~ [0.1*prec, 0.3*prec] бит: при 30000 знаках
// u=9990 - 0.5с, а u>=11000 - минуты и часы (проверено). Поэтому:
//   u <  0.08*prec  -> expint (быстро и точно);
//   0.08*prec <= u < 0.13*prec -> степенной ряд (см. e1_ps);
//   u >= 0.13*prec  -> цепная дробь (быстро при больших u, точна по построению).
constexpr double E1_EXPINT_MAX_FRAC = 0.08;  // граница expint, доля prec (бит)
constexpr double E1_CF_MIN_FRAC = 0.13;      // граница цепной дроби, доля prec (бит)
constexpr double E1_CF_MIN = 8.0;            // ниже этого u цепная дробь заведомо медленнее
constexpr long E1_CF_MAXK = 4000000;         // страховочный предел числа членов цепной дроби
constexpr long E1_PS_MAXK = 4000000;         // страховочный предел числа членов ряда
// expint теряет ~2u/ln2 бит относительной точности (абс. ошибка ~2^-pe e^{2u}),
// поэтому нужная pe = nb + 2u/ln2 = tol/ln2 + запас (постоянная; проверено u<=1e4)
constexpr int E1_EXPINT_PAD = 256;           // запас бит к tol/ln2 для expint

// Режим расчёта E1 для u>=1: auto_ - expint с правильной точностью (быстро везде),
// cf - цепная дробь (контрольный путь; на малых u медленная, но точна по построению).
enum class E1Mode { auto_, expint_, cf };
E1Mode g_e1_mode = E1Mode::auto_;
constexpr int RECOMP = 256;            // как часто e^{-t} пересчитывается заново

// 100 знаков lambda (OEIS A084945), только для --check
const char* KNOWN_100 =
    "0.6243299885435508709929363831008372441796426201805292869735519024"
    "956380888551132544624602761955398688";

// ---------------------------------------------------------------------------
// Обёртка над arb_t (RAII) - неявно приводится к arb_ptr / arb_srcptr
// ---------------------------------------------------------------------------
struct Arb {
    mutable arb_t v;
    Arb() { arb_init(v); }
    ~Arb() { arb_clear(v); }
    Arb(const Arb&) = delete;
    Arb& operator=(const Arb&) = delete;
    operator arb_ptr() const { return v; }
};

// ---------------------------------------------------------------------------
// Параметры квадратуры (всё в double - только для выбора границ и точностей)
// ---------------------------------------------------------------------------
struct Est { double lnf, xf, uf, emf; };

// Верхняя оценка ln f для узла t: F(u) <= min(e^gamma u, e^-u).  xf = ln u, uf = u, emf = e^-t.
Est est(double tf) {
    Est e;
    e.emf = std::exp(-tf);
    e.xf = tf - e.emf;
    e.uf = e.xf < 700.0 ? std::exp(e.xf) : std::numeric_limits<double>::infinity();
    double ln_j = e.xf + std::log1p(e.emf);
    e.lnf = ln_j + std::min(EULER + e.xf, -e.uf);
    return e;
}

void cutoffs(double thr, double& t_lo, double& t_hi) {
    double lo = -40.0, hi = 0.0;
    for (int k = 0; k < 200; ++k) {
        double m = 0.5 * (lo + hi);
        if (est(m).lnf < -thr) lo = m; else hi = m;
    }
    t_lo = hi;
    lo = 0.0; hi = 40.0;
    for (int k = 0; k < 200; ++k) {
        double m = 0.5 * (lo + hi);
        if (est(m).lnf > -thr) lo = m; else hi = m;
    }
    t_hi = lo;
}

double t_of_u(double u) {  // t: exp(t - e^{-t}) = u
    double lx = std::log(u), lo = -40.0, hi = 40.0;
    for (int k = 0; k < 200; ++k) {
        double m = 0.5 * (lo + hi);
        if (m - std::exp(-m) < lx) lo = m; else hi = m;
    }
    return lo;
}

struct Plan {
    long digits, guard;
    double cst, tgt, tol, h;
    slong prec;
    long k_lo, k_hi;

    Plan(long digits_, long guard_, double cst_) : digits(digits_), guard(guard_), cst(cst_) {
        if (digits < 10) throw std::invalid_argument("digits должно быть >= 10");
        if (cst <= 0) throw std::invalid_argument("const должно быть > 0");
        tgt = (double)(digits + guard) * LN10;     // ln 10^(digits+guard)
        tol = tgt + 30.0;                          // узел считаем до exp(-tol)
        h = cst / tgt;                             // шаг - точный double
        prec = (slong)((double)(digits + guard + 12) * LOG2_10) + 192;
        double t_lo, t_hi;
        cutoffs(tol, t_lo, t_hi);
        k_lo = (long)(-t_lo / h) + 1;
        k_hi = (long)(t_hi / h) + 1;
    }
    long n_nodes() const { return k_lo + k_hi + 1; }
    std::string meta(long chunk) const {
        char buf[600];
        std::snprintf(buf, sizeof buf,
                      "version=2\ndigits=%ld\nguard=%ld\nconst=%a\nh=%a\nprec=%ld\nk_lo=%ld\nk_hi=%ld\nchunk=%ld\n",
                      digits, guard, cst, h, (long)prec, k_lo, k_hi, chunk);
        return buf;
    }
};

// ---------------------------------------------------------------------------
// Один узел
// ---------------------------------------------------------------------------

// Наименьшее K с K!/u^K < 2^-bits (члены убывают при k < u); 0 - не нашлось.
int asym_terms(double u, double bits) {
    double thr = -bits * LN2, lu = std::log(u), lt = 0.0;
    long top = std::min<long>(ASYM_KMAX, (long)u - 1);
    for (long k = 1; k <= top; ++k) {
        lt += std::log((double)k) - lu;
        if (lt < thr) return (int)k;
    }
    return 0;
}

// Число членов K ряда V(u) при ln u = x < 0: u^K/(K K!) < 2^-bits; 0 - не нашлось.
int series_terms(double x, double bits) {
    double thr = -bits * LN2, lg = 0.0;
    for (int k = 2; k <= SERIES_KMAX; ++k) {
        lg += std::log((double)k);
        if (k * x - lg - std::log((double)k) < thr) return k;
    }
    return 0;
}

struct Work {  // временные arb на поток, чтобы не создавать их в цикле
    Arb x, u, em, j1, e1, w, S, acc, p, mu, t1, t2, res, H, one;
    Work() { arb_one(one.v); }
};

void expint1(arb_ptr res, arb_srcptr u, slong prec, Work& W) {
#ifdef GD_USE_EI
    arb_neg(W.t2.v, u);
    arb_hypgeom_ei(res, W.t2.v, prec);
    arb_neg(res, res);
#else
    arb_hypgeom_expint(res, W.one.v, u, prec);
#endif
}

// E1(z) = e^{-z} / f,  f = z + 1/(1 + 1/(z + 2/(1 + 2/(z + 3/(1 + ...))))).
// a_k = ceil(k/2); b_k = z (k чётное), 1 (k нечётное). Модифицированный Lentz.
// Для БОЛЬШИХ z сходится за ~O(bits^2/z) итераций (быстро там, где expint «висит»).
// Прогоняем до |delta-1| < 2^-(bits+8), не более maxK членов. Возвращает число членов.
long e1_cf(arb_ptr res, arb_srcptr z, slong bits, long maxK, slong prec) {
    Arb f, C, D, delta, a, b, t, tiny, one;
    arb_set_ui(tiny, 1); arb_mul_2exp_si(tiny, tiny, -prec - 20);
    arb_one(one);
    arb_set(f, z);
    if (arb_is_zero(f)) arb_set(f, tiny);
    arb_set(C, f);
    arb_zero(D);
    long used = maxK;
    for (long k = 1; k <= maxK; ++k) {
        ulong ak = (ulong)((k + 1) / 2);
        arb_set_ui(a, ak);
        if (k % 2 == 0) arb_set(b, z); else arb_set_ui(b, 1);
        arb_mul(D, D, a, prec); arb_add(D, D, b, prec);
        if (arb_is_zero(D)) arb_set(D, tiny);
        arb_div(C, a, C, prec); arb_add(C, C, b, prec);
        if (arb_is_zero(C)) arb_set(C, tiny);
        arb_inv(D, D, prec);
        arb_mul(delta, C, D, prec);
        arb_mul(f, f, delta, prec);
        arb_sub(t, delta, one, prec); arb_abs(t, t);
        if (k >= 8 && arf_cmpabs_2exp_si(arb_midref(t.v), -bits - 8) < 0) { used = k; break; }
    }
    arb_neg(t, z);
    arb_exp(t, t, prec);
    arb_div(res, t, f, prec);
    return used;
}

// E1(z) = -gamma - ln z - sum_{k>=1} (-z)^k/(k*k!) с АБСОЛЮТНОЙ точностью 2^-prec.
// Члены растут до ~e^z/z, поэтому считаем с запасом ~z/ln2 бит на сокращение.
// Число членов K ~ до z^K/(K*K!) < 2^-prec (при prec~1e5: z=1e3 -> 3e4, z=2e4 -> 1.2e5).
long e1_ps(arb_ptr res, arb_srcptr z, double zf, slong prec, long maxK) {
    slong pp = prec + (slong)(zf / LN2) + 32;
    Arb mz, p, s, term, gam, lz;
    arb_neg(mz, z);
    arb_set(p, mz);  // p = (-z)^k/k!
    arb_set(s, p);   // сумма, член k=1 равен -z
    long used = maxK;
    for (long k = 2; k <= maxK; ++k) {
        arb_mul(p, p, mz, pp);
        arb_div_ui(p, p, (ulong)k, pp);
        arb_div_ui(term, p, (ulong)k, pp);
        arb_add(s, s, term, pp);
        if (arf_cmpabs_2exp_si(arb_midref(term.v), -pp - 8) < 0) { used = k; break; }
    }
    arb_const_euler(gam, pp);
    arb_log(lz, z, pp);
    arb_add(res, gam, lz, pp);
    arb_add(res, res, s, pp);
    arb_neg(res, res);
    return used;
}

// E1(u_ball) для u >= 1 — запасной путь, когда асимптотический ряд не сходится.
// expint даёт относительную точность 2^-pe, но его АБСОЛЮТНАЯ ошибка ~2^-pe e^{2u},
// поэтому нужна pe = nb + 2u/ln2 = tol/ln2 (постоянная); при такой pe он быстр
// до u ~ 0.08*prec. Выше - степенной ряд (до 0.22*prec), затем цепная дробь.
// Режим --e1 cf принудительно считает цепной дробью (контрольный путь).
void e1_auto(arb_ptr res, arb_srcptr u_ball, double uf, const Plan& pl, slong nb, slong P, Work& W) {
    if (2.0 * uf >= pl.tol) {
        arb_zero(res);
        return;
    }
    double u_cf = std::max(E1_CF_MIN, E1_CF_MIN_FRAC * (double)pl.prec);
    double u_ps = E1_EXPINT_MAX_FRAC * (double)pl.prec;
    bool cf = (g_e1_mode == E1Mode::cf) ? (uf >= E1_CF_MIN) : (g_e1_mode == E1Mode::auto_ && uf >= u_cf);
    if (cf) {
        // Цепная дробь даёт относительную точность 2^-pe, pe = nb + 48 достаточно.
        slong pe = std::max<slong>(nb + 48, 64);
        e1_cf(res, u_ball, pe, E1_CF_MAXK, pe + 64);
    } else if (g_e1_mode == E1Mode::auto_ && uf >= u_ps) {
        e1_ps(res, u_ball, uf, P + 128, E1_PS_MAXK);
    } else {
        slong pe = (slong)(pl.tol / LN2) + E1_EXPINT_PAD;
        if (pe > pl.prec + 4096) pe = pl.prec + 4096;
        expint1(res, u_ball, pe, W);
    }
}

// res = F(u) * u (1 + e^{-t}),  t = i h.  em = e^{-t} (полная точность), EG = e^gamma.
// false - узел меньше exp(-tol), его не учитываем.
bool node(slong i, const Plan& pl, Work& W, arb_srcptr em, arb_srcptr EG) {
    double tf = (double)i * pl.h;
    Est e = est(tf);
    double need = pl.tol + e.lnf;  // нужная относительная точность узла, нат
    if (need < 0.0) return false;
    slong P = std::min<slong>(std::max<slong>((slong)(need / LN2) + NODE_PAD_BITS, 128), pl.prec);
    // Для u>=1 узел содержит exp(-(u + E1(u))); аргумент ~u, поэтому относительная
    // точность P = (tol - u + ln u)/ln2 уже обеспечивает абсолютную ~e^{-tol}. Этот
    // блок лишь сужает РАДИУС шара узла (и итоговой суммы), не меняя точность expint
    // (там точность ограничена min(nb+48, P+32), т.е. фактически nb+48).
    if (e.xf >= 0.0) {
        slong Pabs = (slong)((pl.tol + e.xf) / LN2) + NODE_PAD_BITS;
        if (Pabs > P) P = Pabs;
    }

    arb_mul_si(W.x, W.H, i, P);
    arb_sub(W.x, W.x, em, P);  // x = ln u
    arb_exp(W.u, W.x, P);
    arb_add_ui(W.j1, em, 1, P);  // du/dt = u * j1

    if (e.xf < 0.0) {  // ---- u < 1
        int K = series_terms(e.xf, (double)(P + 8));
        if (K) {  // F = e^gamma u exp(V(u)), V = -2u + sum_{k>=2} (-u)^k/(k k!)
            arb_neg(W.mu, W.u);
            arb_set(W.p, W.mu);
            arb_add(W.acc, W.mu, W.p, P);
            for (int k = 2; k <= K; ++k) {
                arb_mul(W.p, W.p, W.mu, P);
                arb_div_ui(W.p, W.p, (ulong)k, P);
                arb_div_ui(W.t1, W.p, (ulong)k, P);
                arb_add(W.acc, W.acc, W.t1, P);
            }
            arb_exp(W.acc, W.acc, P);
            arb_mul(W.res, EG, W.u, P);
            arb_mul(W.res, W.res, W.u, P);
            arb_mul(W.res, W.res, W.j1, P);
            arb_mul(W.res, W.res, W.acc, P);
            return true;
        }
        expint1(W.e1.v, W.u.v, P + 32, W);
        arb_add(W.t1, W.u, W.e1, P);
        arb_neg(W.t1, W.t1);
        arb_exp(W.t1, W.t1, P);
        arb_mul(W.res, W.t1, W.u, P);
        arb_mul(W.res, W.res, W.j1, P);
        return true;
    }

    // ---- u >= 1: нужная относительная точность E1 (в битах)
    double nb = (pl.tol - 2.0 * e.uf) / LN2;
    if (nb < E1_SKIP_BITS) {  // E1 ~ e^-u/u меньше допуска
        arb_neg(W.t1, W.u);
        arb_exp(W.t1, W.t1, P);
        arb_mul(W.res, W.t1, W.u, P);
        arb_mul(W.res, W.res, W.j1, P);
        return true;
    }
    int K = asym_terms(e.uf, nb + 20.0);
    if (K) {  // E1 = e^-u/u * sum_{k<K} (-1)^k k!/u^k на малой точности
        arb_neg(W.t1, W.u);
        arb_exp(W.w, W.t1, P);
        slong pa = std::max<slong>(64, (slong)nb + 48);
        arb_one(W.S);
        for (int k = K - 1; k >= 1; --k) {  // S = 1 - (k/u) S
            arb_mul_ui(W.S, W.S, (ulong)k, pa);
            arb_div(W.S, W.S, W.u, pa);
            arb_sub_ui(W.S, W.S, 1, pa);
            arb_neg(W.S, W.S);
        }
        arb_mul(W.e1, W.w, W.S, pa);
        arb_div(W.e1, W.e1, W.u, pa);
        arb_neg(W.t1, W.e1);
        arb_exp(W.t1, W.t1, P);
        arb_mul(W.res, W.w, W.t1, P);
        arb_mul(W.res, W.res, W.u, P);
        arb_mul(W.res, W.res, W.j1, P);
        return true;
    }
    // E1 нужна с АБСОЛЮТНОЙ точностью ~e^{-tol}; expint даёт относительную.
    // e1_auto сам выбирает точность и режим окна, чтобы обойти медленные пути FLINT.
    e1_auto(W.e1, W.u, e.uf, pl, (slong)nb, P, W);
    if (arb_is_zero(W.e1)) {  // E1 ~ 0: узел = e^{-u} * u * (1+e^{-t})
        arb_neg(W.t1, W.u);
        arb_exp(W.t1, W.t1, P);
        arb_mul(W.res, W.t1, W.u, P);
        arb_mul(W.res, W.res, W.j1, P);
        return true;
    }
    arb_add(W.t1, W.u, W.e1, P);
    arb_neg(W.t1, W.t1);
    arb_exp(W.t1, W.t1, P);
    arb_mul(W.res, W.t1, W.u, P);
    arb_mul(W.res, W.res, W.j1, P);
    return true;
}

// Эталонный узел: без сокращений точности и без рядов, E1 только через expint.
void node_ref(slong i, const Plan& pl, Work& W, Arb& out) {
    slong prec = pl.prec + 96;
    Arb t, em, u, e1, tmp;
    arb_set_d(W.H, pl.h);
    arb_mul_si(t, W.H, i, prec);
    arb_neg(tmp, t);
    arb_exp(em, tmp, prec);
    arb_sub(tmp, t, em, prec);
    arb_exp(u, tmp, prec);
    Est e = est((double)i * pl.h);
    if (e.xf >= 0.0 && 2.0 * e.uf > pl.tol + 60.0) {
        arb_zero(e1);
    } else if (e.xf >= 0.0) {
        // Как в node(): абсолютная точность E1 ~ e^{-tol}. Используем тот же e1_auto
        // (точный double-аргумент + обход медленных окон), но с запасом по точности,
        // чтобы эталон был не грубее рабочего пути.
        slong nbr = (slong)((pl.tol - 2.0 * e.uf) / LN2);
        e1_auto(e1, u, e.uf, pl, nbr, prec + 288, W);
    } else {
        expint1(e1, u, prec, W);
    }
    // Аргумент exp(-(u+e1)) ~ u, поэтому и здесь нужна абсолютная точность e^{-tol}.
    slong pout = prec;
    if (e.xf >= 0.0) pout = std::max<slong>(pout, (slong)((pl.tol + e.xf) / LN2) + NODE_PAD_BITS + 192);
    arb_add(tmp, u, e1, pout);
    arb_neg(tmp, tmp);
    arb_exp(tmp, tmp, pout);
    arb_mul(tmp, tmp, u, pout);
    arb_add_ui(em, em, 1, pout);
    arb_mul(out, tmp, em, pout);
}

// log10 |середины|, -inf для нуля
double lg10_mid(arb_srcptr x) {
    Arb m, t;
    arb_get_mid_arb(m, x);
    if (arb_is_zero(m)) return -std::numeric_limits<double>::infinity();
    arb_abs(m, m);
    arb_log(t, m, 64);
    return arf_get_d(arb_midref(t.v), ARF_RND_NEAR) / LN10;
}

// ---------------------------------------------------------------------------
// Куски узлов, чекпоинты
// ---------------------------------------------------------------------------
struct ChunkRes {
    slong lo = 0, hi = 0;
    Arb e, o;  // суммы по чётным / нечётным индексам
};

void run_chunk(slong lo, slong hi, const Plan& pl, const Arb& EG, ChunkRes& out) {
    Work W;
    arb_set_d(W.H, pl.h);
    slong prec = pl.prec;
    Arb R, em, tmp, s0, s1;
    arb_neg(tmp, W.H);
    arb_exp(R, tmp, prec);
    for (slong i = lo, n = 0; i < hi; ++i, ++n) {
        if (n % RECOMP == 0) {
            arb_mul_si(tmp, W.H, -i, prec);
            arb_exp(em, tmp, prec);
        } else {
            arb_mul(em, em, R, prec);
        }
        if (node(i, pl, W, em, EG)) {
            Arb& s = (i & 1) ? s1 : s0;
            arb_add(s, s, W.res, prec);
        }
    }
    out.lo = lo;
    out.hi = hi;
    arb_set(out.e, s0);
    arb_set(out.o, s1);
}

std::string dump(arb_srcptr x) {
    char* s = arb_dump_str(x);
    std::string r(s);
    flint_free(s);
    return r;
}

void load(arb_ptr x, const std::string& s) {
    if (arb_load_str(x, s.c_str())) throw std::runtime_error("arb_load_str: повреждённый чекпоинт");
}

void write_atomic(const fs::path& p, const std::string& text) {
    fs::path tmp = p;
    tmp += ".tmp";
    {
        std::ofstream f(tmp, std::ios::binary);
        f << text;
        if (!f) throw std::runtime_error("не удалось записать " + tmp.string());
    }
    fs::rename(tmp, p);
}

void save_chunk(const std::string& dir, const ChunkRes& c) {
    std::ostringstream os;
    os << c.lo << ' ' << c.hi << '\n' << dump(c.e) << '\n' << dump(c.o) << '\n';
    write_atomic(fs::path(dir) / ("c_" + std::to_string(c.lo) + ".txt"), os.str());
}

std::unique_ptr<ChunkRes> load_chunk(const std::string& dir, slong lo, slong hi) {
    fs::path p = fs::path(dir) / ("c_" + std::to_string(lo) + ".txt");
    std::ifstream f(p, std::ios::binary);
    if (!f) return nullptr;
    long a, b;
    std::string l1, l2, l3;
    if (!(f >> a >> b)) return nullptr;
    std::getline(f, l1);
    std::getline(f, l2);
    std::getline(f, l3);
    if (a != lo || b != hi || l2.empty() || l3.empty()) return nullptr;
    auto c = std::make_unique<ChunkRes>();
    c->lo = lo;
    c->hi = hi;
    try {
        load(c->e, l2);
        load(c->o, l3);
    } catch (const std::exception&) {
        return nullptr;
    }
    return c;
}

struct Bound { slong lo, hi; };

// Считает все куски (с чекпоинтами), возвращает их отсортированными по lo.
std::vector<std::unique_ptr<ChunkRes>> run_nodes(const Plan& pl, int jobs, long chunk, const std::string& ckpt,
                                                  const Arb& EG) {
    slong n0 = -pl.k_lo, n1 = pl.k_hi + 1;
    long total = n1 - n0;
    if (chunk <= 0) chunk = std::max<long>(50, std::min<long>(1500, total / ((long)std::max(1, jobs) * 24)));
    std::vector<Bound> bounds;
    for (slong lo = n0; lo < n1; lo += chunk) bounds.push_back({lo, std::min<slong>(lo + chunk, n1)});

    std::vector<std::unique_ptr<ChunkRes>> done;
    std::vector<Bound> todo;
    if (!ckpt.empty()) {
        fs::create_directories(ckpt);
        fs::path mp = fs::path(ckpt) / "meta.txt";
        std::string meta = pl.meta(chunk);
        if (fs::exists(mp)) {
            std::ifstream f(mp, std::ios::binary);
            std::stringstream ss;
            ss << f.rdbuf();
            if (ss.str() != meta)
                throw std::runtime_error(ckpt + ": чекпоинты другого прогона (другие digits/guard/const/chunk)");
        } else {
            write_atomic(mp, meta);
        }
    }
    for (const Bound& b : bounds) {
        std::unique_ptr<ChunkRes> c;
        if (!ckpt.empty()) c = load_chunk(ckpt, b.lo, b.hi);
        if (c) done.push_back(std::move(c)); else todo.push_back(b);
    }
    std::sort(todo.begin(), todo.end(), [](const Bound& a, const Bound& b) { return a.lo > b.lo; });  // дорогие - первыми

    if (!done.empty())
        std::fprintf(stderr, "  из чекпоинта взято кусков: %zu/%zu\n", done.size(), bounds.size());

    std::atomic<size_t> next{0};
    std::mutex mu;
    std::string error;
    size_t finished = 0;
    auto t0 = std::chrono::steady_clock::now();

    auto worker = [&]() {
        try {
            for (;;) {
                size_t k = next.fetch_add(1);
                if (k >= todo.size()) break;
                {
                    std::lock_guard<std::mutex> lk(mu);
                    if (!error.empty()) break;
                }
                auto c = std::make_unique<ChunkRes>();
                run_chunk(todo[k].lo, todo[k].hi, pl, EG, *c);
                std::lock_guard<std::mutex> lk(mu);
                if (!ckpt.empty()) save_chunk(ckpt, *c);
                done.push_back(std::move(c));
                ++finished;
                double el = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
                double eta = el / (double)finished * (double)(todo.size() - finished);
                std::fprintf(stderr, "\r  куски %zu/%zu  (новых %zu/%zu)  прошло %8.0fс  осталось ~%8.0fс   ",
                             done.size(), bounds.size(), finished, todo.size(), el, eta);
                std::fflush(stderr);
            }
        } catch (const std::exception& ex) {
            std::lock_guard<std::mutex> lk(mu);
            if (error.empty()) error = ex.what();
        }
        flint_cleanup();  // освободить кэши FLINT этого потока
    };

    int nt = (int)std::min<size_t>((size_t)std::max(1, jobs), std::max<size_t>(1, todo.size()));
    if (nt <= 1) {
        worker();
    } else {
        std::vector<std::thread> th;
        for (int k = 0; k < nt; ++k) th.emplace_back(worker);
        for (auto& t : th) t.join();
    }
    if (!todo.empty()) std::fprintf(stderr, "\n");
    if (!error.empty()) throw std::runtime_error(error);
    std::sort(done.begin(), done.end(), [](const auto& a, const auto& b) { return a->lo < b->lo; });
    return done;
}

// ---------------------------------------------------------------------------
// Сборка результата
// ---------------------------------------------------------------------------
struct Combined {
    Arb val;    // lambda (шар)
    Arb delta;  // 2h |S_even - S_odd|
};

void combine(const Plan& pl, const std::vector<std::unique_ptr<ChunkRes>>& chunks, Combined& out) {
    slong P2 = pl.prec + 256;
    Arb Se, So, H, t;
    for (const auto& c : chunks) {
        arb_add(Se, Se, c->e, P2);
        arb_add(So, So, c->o, P2);
    }
    arb_set_d(H, pl.h);
    arb_add(t, Se, So, P2);
    arb_mul(out.val, t, H, P2);
    arb_sub(t, Se, So, P2);
    arb_abs(t, t);
    arb_mul(out.delta, t, H, P2);
    arb_mul_ui(out.delta, out.delta, 2, P2);
}

// Печать с округлением до digits знаков. Как в Python-эталоне (to_decimal),
// печатается СЕРЕДИНА шара: радиус не является жёстким фильтром (он отдельно
// сообщается в отчёте и влияет на код возврата через q_rad), поэтому широкий
// шар сам по себе деление не срывает.
std::string to_decimal(arb_srcptr v, long digits) {
    Arb mid;
    arb_get_mid_arb(mid, v);
    char* s = arb_get_str(mid, digits + 25, ARB_STR_NO_RADIUS);
    std::string str(s);
    flint_free(s);
    if (str.compare(0, 2, "0.") != 0) throw std::runtime_error("неожиданный формат arb_get_str: " + str.substr(0, 40));
    std::string d = str.substr(2);
    if ((long)d.size() < digits + 5)
        throw std::runtime_error("arb_get_str напечатал слишком мало цифр");
    bool up = d[digits] >= '5';
    d.resize(digits);
    if (up) {
        long k = digits - 1;
        while (k >= 0) {
            if (d[k] == '9') { d[k] = '0'; --k; } else { ++d[k]; break; }
        }
        if (k < 0) throw std::runtime_error("перенос за старший разряд");
    }
    return "0." + d;
}

// ---------------------------------------------------------------------------
// Аудит: рабочий путь против эталонного
// ---------------------------------------------------------------------------
int audit(const Plan& pl, int n, const Arb& EG) {
    slong n0 = -pl.k_lo, n1 = pl.k_hi + 1;
    std::set<slong> idx = {n0, n1 - 1};
    for (int k = 0; k < n; ++k) idx.insert(n0 + (n1 - 1 - n0) * k / std::max(1, n - 1));
    const double us[] = {1e-300, 1e-30, 1e-3, 0.5, 1.0, 2.0, 50.0, 1e3, pl.tol / 6, pl.tol / 3,
                         0.45 * pl.tol, 0.55 * pl.tol, pl.tol};
    for (double u : us) {
        slong i = (slong)std::llround(t_of_u(u) / pl.h);
        if (i >= n0 && i < n1) idx.insert(i);
    }
    double tol_lg = -(pl.tol + 10.0) / LN10;  // допуск на узел, log10
    std::fprintf(stderr, "аудит: %zu узлов, допуск на узел 10^%.0f\n", idx.size(), tol_lg);
    Work W, W2;
    arb_set_d(W.H, pl.h);
    Arb em, tmp, ref, diff;
    int bad = 0;
    for (slong i : idx) {
        arb_mul_si(tmp, W.H, -i, pl.prec);
        arb_exp(em, tmp, pl.prec);
        if (!node(i, pl, W, em, EG)) {
            std::fprintf(stderr, "  i=%9ld пропущен (меньше допуска)\n", (long)i);
            continue;
        }
        node_ref(i, pl, W2, ref);
        arb_sub(diff, W.res, ref, pl.prec + 96);
        double d = lg10_mid(diff);
        bool ok = d <= tol_lg;
        bad += !ok;
        std::fprintf(stderr, "  i=%9ld ln u=%11.2f  |fast-ref| = 10^%9.1f  %s\n", (long)i, est((double)i * pl.h).xf, d,
                     ok ? "ok" : "ОШИБКА");
    }
    return bad;
}

}  // namespace

// ---------------------------------------------------------------------------
int main(int argc, char** argv) {
    long digits = 0, guard = 30, chunk = 0;
    double cst = DEFAULT_CONST;
    int jobs = 0, n_audit = 0;
    std::string ckpt, out;
    bool plan_only = false, check = false;

    auto usage = [&]() {
        std::fprintf(stderr,
                     "использование: %s DIGITS [--guard N] [--const X] [--jobs N] [--chunk N] [--ckpt DIR]\n"
                     "                [--out FILE] [--e1 auto|expint|cf] [--audit N] [--plan] [--check]\n",
                     argv[0]);
        return 1;
    };
    for (int k = 1; k < argc; ++k) {
        std::string a = argv[k];
        auto val = [&]() -> const char* {
            if (k + 1 >= argc) { std::fprintf(stderr, "для %s нужно значение\n", a.c_str()); std::exit(1); }
            return argv[++k];
        };
        if (a == "--guard") guard = std::atol(val());
        else if (a == "--const") cst = std::atof(val());
        else if (a == "--jobs") jobs = std::atoi(val());
        else if (a == "--chunk") chunk = std::atol(val());
        else if (a == "--ckpt") ckpt = val();
        else if (a == "--out") out = val();
        else if (a == "--e1") {
            std::string m = val();
            if (m == "cf") g_e1_mode = E1Mode::cf;
            else if (m == "expint") g_e1_mode = E1Mode::expint_;
            else if (m == "auto") g_e1_mode = E1Mode::auto_;
            else { std::fprintf(stderr, "--e1: ожидается auto|expint|cf\n"); return 1; }
        }
        else if (a == "--audit") n_audit = std::atoi(val());
        else if (a == "--plan") plan_only = true;
        else if (a == "--check") check = true;
        else if (a[0] != '-' && digits == 0) digits = std::atol(a.c_str());
        else return usage();
    }
    if (digits == 0) return usage();
    if (jobs <= 0) jobs = (int)std::max(1u, std::thread::hardware_concurrency());

    try {
        Plan pl(digits, guard, cst);
        std::fprintf(stderr,
                     "lambda: %ld знаков (+%ld запас), узлов %ld (u<1: %ld, u>1: %ld), h=%.4e (const=%g), "
                     "полная точность %ld бит, потоков %d\n",
                     pl.digits, pl.guard, pl.n_nodes(), pl.k_lo, pl.k_hi + 1, pl.h, pl.cst, (long)pl.prec, jobs);
        if (plan_only) return 0;

        Arb EG, g;
        arb_const_euler(g, pl.prec);
        arb_exp(EG, g, pl.prec);

        if (check && n_audit == 0) n_audit = 12;
        if (n_audit > 0 && audit(pl, n_audit, EG) != 0) {
            std::fprintf(stderr, "аудит НЕ пройден - результат считать нельзя\n");
            return 3;
        }

        auto t0 = std::chrono::steady_clock::now();
        auto chunks = run_nodes(pl, jobs, chunk, ckpt, EG);
        Combined cb;
        combine(pl, chunks, cb);
        std::string dec = to_decimal(cb.val, pl.digits);
        double dt = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();

        Arb rad;
        arb_get_rad_arb(rad, cb.val);
        double q_est = 2.0 * lg10_mid(cb.delta);  // log10 оценки ошибки квадратуры (~Delta^2)
        double q_rad = lg10_mid(rad);             // log10 радиуса шара результата
        double need = -(double)(pl.digits + 3);
        bool ok = q_est < need && q_rad < need;
        std::fprintf(stderr, "готово за %.0fс;  оценка ошибки квадратуры ~10^%.0f,  радиус шара ~10^%.0f  "
                             "(нужно < 10^%.0f)  -> %s\n",
                     dt, q_est, q_rad, need, ok ? "OK" : "НЕДОСТАТОЧНО");

        if (check && pl.digits >= 100) {
            bool same = dec.compare(0, 99, std::string(KNOWN_100).substr(0, 99)) == 0;
            std::fprintf(stderr, "сверка с OEIS A084945 (первые 99 символов): %s\n", same ? "OK" : "НЕ СОВПАЛО");
            ok = ok && same;
        }

        char hdr[400];
        std::snprintf(hdr, sizeof hdr,
                      "# Golomb-Dickman constant lambda, %ld digits; Takahasi-Mori, nodes=%ld, const=%g; "
                      "quad.err.est=10^%.0f, ball.radius=10^%.0f%s\n",
                      pl.digits, pl.n_nodes(), pl.cst, q_est, q_rad, ok ? "" : "  *** НЕ ПОДТВЕРЖДЕНО ***");
        if (!out.empty()) {
            std::ofstream f(out, std::ios::binary);
            f << hdr << dec << '\n';
            if (!f) throw std::runtime_error("не удалось записать " + out);
            std::fprintf(stderr, "записано %s\n", out.c_str());
        } else {
            std::cout << hdr << dec << '\n';
        }
        return ok ? 0 : 2;
    } catch (const std::exception& ex) {
        std::fprintf(stderr, "ошибка: %s\n", ex.what());
        return 1;
    }
}
