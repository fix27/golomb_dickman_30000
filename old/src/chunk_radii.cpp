// Проверка сертификации по кускам: читает чекпоинты golomb_dickman (c_<lo>.txt)
// и печатает радиус шара каждого куска (log10) и суммарный радиус.
//
// Формат куска (см. save_chunk в golomb_dickman.cpp):
//     <lo> <hi>\n
//     <arb_dump_str(e)>\n
//     <arb_dump_str(o)>\n
//
// Сборка:
//     g++ -O3 -std=c++17 chunk_radii.cpp -o chunk_radii -lflint -lgmp -lmpfr
// Запуск:
//     chunk_radii run_v3 30000
// Второй аргумент (число знаков) включает проверку порога 10^-(digits+3);
// код возврата 2, если хоть один кусок шире.

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

#if __has_include(<flint/arb.h>)
#include <flint/flint.h>
#include <flint/arb.h>
#else
#include <flint.h>
#include <arb.h>
#endif

namespace fs = std::filesystem;

static const double LN10 = 2.302585092994045684017991454684;

// log10 |x| для середины шара
static double lg10_mid(arb_srcptr x) {
    arb_t t;
    arb_init(t);
    double r;
    if (arb_is_zero(x)) {
        r = -INFINITY;
    } else {
        arb_abs(t, x);
        arb_log(t, t, 64);
        r = arf_get_d(arb_midref(t), ARF_RND_NEAR) / LN10;
    }
    arb_clear(t);
    return r;
}

int main(int argc, char** argv) {
    if (argc < 2) {
        std::fprintf(stderr, "использование: %s DIR [DIGITS]\n", argv[0]);
        return 1;
    }
    const fs::path dir = argv[1];
    const long digits = argc > 2 ? std::atol(argv[2]) : 0;

    struct Rec {
        long lo;
        double rad;
    };
    std::vector<Rec> recs;
    arb_t x, rad, radsum;
    arb_init(x);
    arb_init(rad);
    arb_init(radsum);
    arb_zero(radsum);

    long bad = 0;
    for (const auto& ent : fs::directory_iterator(dir)) {
        const std::string name = ent.path().filename().string();
        if (name.size() < 6 || name.compare(0, 2, "c_") != 0) continue;
        if (name.compare(name.size() - 4, 4, ".txt") != 0) continue;
        std::ifstream f(ent.path(), std::ios::binary);
        if (!f) continue;
        long lo = 0, hi = 0;
        std::string dummy, s1;
        if (!(f >> lo >> hi)) continue;
        std::getline(f, dummy);   // остаток строки "lo hi"
        std::getline(f, s1);      // шар частичной суммы
        if (s1.empty()) continue;
        if (arb_load_str(x, s1.c_str())) {
            std::fprintf(stderr, "%s: не разобран шар\n", name.c_str());
            continue;
        }
        arb_get_rad_arb(rad, x);
        arb_add(radsum, radsum, rad, 64);
        const double r = lg10_mid(rad);
        if (digits > 0 && r > -(double)digits - 3.0) ++bad;
        recs.push_back({lo, r});
    }
    arb_clear(x);

    std::sort(recs.begin(), recs.end(), [](const Rec& a, const Rec& b) { return a.rad < b.rad; });
    std::printf("chunks=%zu\n", recs.size());
    if (digits > 0) std::printf("worse than 1e-%ld: %ld\n", digits + 3, bad);
    for (const Rec& r : recs) std::printf("  lo=%9ld  rad=1e%.0f\n", r.lo, r.rad);
    std::printf("sum radius: 1e%.0f\n", lg10_mid(radsum));
    arb_clear(rad);
    arb_clear(radsum);
    return bad == 0 ? 0 : 2;
}
