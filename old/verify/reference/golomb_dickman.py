#!/usr/bin/env python3
"""Постоянная Голомба-Дикмана  lambda = integral_0^inf exp(-u - E1(u)) du
= 0.62432998854355087099293638310083724417964262018...   (версия для десятков тысяч знаков)

Что здесь устроено иначе, чем в первой версии, и почему
-------------------------------------------------------
1. Подстановка Такахаси-Мори  u = exp(t - exp(-t)),  du = u (1 + exp(-t)) dt.
   Прежняя u = exp(a sinh t) на правом конце имеет фазу a cosh(t) sin(y); при сдвиге
   t -> t + iy она выходит за pi/2, Re u < 0 и e^{-u} начинает расти.  Поэтому полоса
   аналитичности сужалась как d ~ pi / (2 sqrt(a^2 + ln^2 T)), и шаг h = const/T должен был
   уменьшаться с ростом числа знаков: const = 1.2 давал только ~11 700 верных знаков из
   15 000.  У подстановки Такахаси-Мори u ~ e^t, фаза -> y, полоса d -> pi/2 и показатель
   ошибки K = -log10(err) * const / (digits + guard) остаётся ~9.8 при любом числе знаков
   (измерено при 100, 300 и 600 знаках).  Шаг h = const/T, const = 8 даёт запас ~1.2 по
   показателю и в ~6 раз меньше узлов, чем прежняя схема при безопасном шаге.

2. Оценка ошибки квадратуры по ходу дела.  Чётные и нечётные узлы - это два грубых
   трапециевидных правила с шагом 2h, сдвинутых на полшага.  Их ошибки почти равны и
   противоположны, а ошибка полного правила ~ (ошибка грубого)^2.  Поэтому
   Delta = 2h |S_even - S_odd|  даёт оценку ошибки h*(S_even + S_odd) порядка Delta^2
   (на практике консервативную).  Программа печатает её и падает с кодом 2, если
   Delta^2 > 10^-(digits+3).  Стоит это ровно ноль дополнительных вычислений.

3. Точность на узле своя.  Значение узла f = F(u) * du/dt нужно с абсолютной погрешностью
   exp(-tol), tol = (digits+guard) ln 10 + 30, поэтому относительная точность узла
   (tol + ln f)/ln 2 бит (+ запас), а не полные digits*3.32 бит.  Для E1 при u > 1 нужна
   относительная точность всего (tol - 2u)/ln 2 бит: при 2u > tol E1 просто отбрасывается.

   ВАЖНАЯ ПОПРАВКА (найдена по шару, а не по середине): flint expint(1) возвращает
   шар с относительным радиусом ~2^(-P/2) (измерено: при P = 100000 и u ~ 11400
   радиус ~2^-50600).  Средняя точка при этом точнее, но сертифицировать надо шар.
   Поэтому на P = plan.prec шар E1 после умножения на f = e^{-u} u даёт лишь
   ~10^-25 вместо 10^-30, и понижать точность E1 нельзя.  Решение: E1 всегда
   считается от u.mid() на 200000 битах (быстрая зона при u < 34650; аргумент
   нельзя брать как float - binary64 даёт только 53 бита).  Наследуемая ошибка u
   (u * 2^-P, P - точность узла) даёт <= 2^-83000, а собственный шар E1 на 200k
   битах < 2^-100000, т.е. вклад E1 в радиус узла пренебрежим.  При u >= 34650
   E1 = 0 (вклад < 10^-(digits+guard)).

4. E1 при больших u считается собственным асимптотическим рядом
   E1(u) = e^{-u}/u * sum_{k<K} (-1)^k k!/u^k  (для u > 0 остаток по модулю не больше
   первого отброшенного члена), на малой точности.  При малых u (u < ~2^-500) вместо
   E1 (log + ряд) используется  F(u) = e^gamma u exp(V(u)), V(u) = T(u) - u  - без log.
   Остальное считает flint (expint(1)) на 200k битах, см. п. 3.

5. Параллелизм и чекпоинты.  Узлы режутся на много мелких кусков (не "по одному на
   процесс": стоимость узла растёт с u), каждый готовый кусок сразу пишется в --ckpt.
   Прерванный прогон продолжается с того же места.

6. Сумма частичных сумм передаётся ТОЧНО (мантисса, порядок), а не через str(arb):
   десятичная печать у arb обрезается до видимых цифр, и при неудачном ctx.dps хвост
   молча терялся бы.  Радиусы шаров складываются и проверяются.

Использование
    python golomb_dickman.py 45000 --plan                      # только параметры и число узлов
    python golomb_dickman.py 300 --check                       # быстрый тест (сверка с OEIS + аудит узлов)
    python golomb_dickman.py 45000 --jobs 16 --ckpt run45k --out lambda_45000.txt
    python golomb_dickman.py 45000 --audit 40 --jobs 1 ...     # выборочная перепроверка узлов эталонным путём

Требует: python-flint (arb, expint) и mpmath (только для сборки и вывода).
"""

from __future__ import annotations

import argparse
import concurrent.futures as cf
import json
import math
import os
import sys
import time

import mpmath as mp
import flint
from flint import arb

# Python >= 3.11 по умолчанию запрещает int <-> str длиннее 4300 цифр; у нас десятки тысяч.
if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(0)

LN2 = math.log(2.0)
LN10 = math.log(10.0)
LOG2_10 = LN10 / LN2
EULER = 0.5772156649015329

DEFAULT_CONST = 8.0      # h = const / T; измерено K ~ 9.8, т.е. запас ~1.2 по показателю
NODE_PAD_BITS = 160      # запас бит на узле: f*2^-P = e^-tol*2^-160 ~= 10^-(digits+guard+48)
SERIES_KMAX = 400        # макс. число членов ряда V(u); иначе берём flint-овский E1
ASYM_KMAX = 3000         # макс. число членов асимптотического ряда E1
E1_SKIP_BITS = -8.0      # если нужная относительная точность E1 хуже 2^8 - E1 отбрасываем
E1_HIGH_BITS = 200000    # expint всегда на 200k битах: при plan.prec шар E1 теряет ~P/2
E1_ZERO_U = 34650.0      # выше - E1 отбрасывается: e^{-2u}/u < 10^-(digits+guard)
RECOMP = 256             # как часто e^{-t} пересчитывается заново (между - рекуррентно)

# 100 знаков lambda (OEIS A084945), только для --check.
KNOWN_100 = (
    "0.6243299885435508709929363831008372441796426201805292869735519024"
    "956380888551132544624602761955398688"
)


# ---------------------------------------------------------------------------
# Параметры квадратуры
# ---------------------------------------------------------------------------

def _est(tf):
    """Грубая (float) оценка ln f для узла t = tf, где f = F(u) du/dt.

    F(u) <= min(e^gamma u, e^-u), поэтому оценка верхняя - на точность узла она влияет
    только в безопасную сторону.  Возвращает (ln f, ln u, u, e^-t)."""
    emf = math.exp(-tf)
    xf = tf - emf
    uf = math.exp(xf) if xf < 700.0 else math.inf
    ln_j = xf + math.log1p(emf)
    return ln_j + min(EULER + xf, -uf), xf, uf, emf


def _cutoffs(thr):
    """Границы по t, за которыми ln f < -thr (там узлы отбрасываются)."""
    lo, hi = -40.0, 0.0
    for _ in range(200):
        m = 0.5 * (lo + hi)
        if _est(m)[0] < -thr:
            lo = m
        else:
            hi = m
    t_lo = hi
    lo, hi = 0.0, 40.0
    for _ in range(200):
        m = 0.5 * (lo + hi)
        if _est(m)[0] > -thr:
            lo = m
        else:
            hi = m
    return t_lo, lo


def _t_of_u(u):
    """t, для которого exp(t - exp(-t)) = u."""
    lx = math.log(u)
    lo, hi = -40.0, 40.0
    for _ in range(200):
        m = 0.5 * (lo + hi)
        if m - math.exp(-m) < lx:
            lo = m
        else:
            hi = m
    return lo


class Plan:
    """Все параметры прогона, зависящие только от (digits, guard, const)."""

    def __init__(self, digits, guard=30, const=DEFAULT_CONST):
        if digits < 10:
            raise ValueError("digits должно быть >= 10")
        if const <= 0:
            raise ValueError("const должно быть > 0")
        self.digits = int(digits)
        self.guard = int(guard)
        self.const = float(const)
        self.tgt = (self.digits + self.guard) * LN10      # ln 10^(digits+guard)
        self.tol = self.tgt + 30.0                        # узел считаем до exp(-tol)
        self.h = self.const / self.tgt                    # шаг - точный double
        self.prec = int((self.digits + self.guard + 12) * LOG2_10) + 192   # бит, полная точность
        t_lo, t_hi = _cutoffs(self.tol)
        self.k_lo = int(-t_lo / self.h) + 1
        self.k_hi = int(t_hi / self.h) + 1

    @property
    def n_nodes(self):
        return self.k_lo + self.k_hi + 1

    def meta(self, chunk):
        return {
            "digits": self.digits, "guard": self.guard, "const": float.hex(self.const),
            "h": float.hex(self.h), "prec": self.prec, "k_lo": self.k_lo,
            "k_hi": self.k_hi, "chunk": chunk, "version": 2,
        }


# ---------------------------------------------------------------------------
# Один узел
# ---------------------------------------------------------------------------

def _asym_terms(u, bits, kmax=ASYM_KMAX):
    """Наименьшее K с K!/u^K < 2^-bits (члены убывают при k < u); иначе None."""
    thr = -bits * LN2
    lu = math.log(u)
    lt = 0.0
    for k in range(1, min(kmax, int(u) - 1) + 1):
        lt += math.log(k) - lu
        if lt < thr:
            return k
    return None


def _series_terms(x, bits, kmax=SERIES_KMAX):
    """Число членов K ряда V(u) при ln u = x < 0: u^K/(K K!) < 2^-bits; иначе None."""
    thr = -bits * LN2
    lg = 0.0
    for k in range(2, kmax + 1):
        lg += math.log(k)
        if k * x - lg - math.log(k) < thr:
            return k
    return None


def _exp1_full(u, uf, prec):
    """E1(u) на полной точности с обходом медленных окон flint.expint.

    flint.expint возвращает шар с относительным радиусом ~2^(-P/2) (измерено:
    при P=100k и u~11400 радиус ~2^-50600).  Поэтому на plan.prec шар E1 после
    умножения на f=e^{-u}u даёт ~10^-25 вместо нужных 10^-30.  Считаем E1 на
    200k бит от midpoint шара u (не float: binary64 даёт лишь 53 бита).  Для
    u < E1_ZERO_U это быстрая зона, а вклад шара E1 становится < 2^-prec
    (наследуемая ошибка u ~ u*2^-P <= 2^-83000).  При u >= E1_ZERO_U
    возвращает 0: вклад E1 в узел там < 2^-prec."""
    if uf >= E1_ZERO_U:
        return arb(0)
    old = flint.ctx.prec
    flint.ctx.prec = E1_HIGH_BITS
    e1 = u.mid().expint(1)
    flint.ctx.prec = old
    return e1


def _node(i, plan, H, em, EG):
    """Значение узла f_i = F(u_i) * u_i (1 + e^{-t_i}), t_i = i h.

    em = e^{-t_i} (arb полной точности), H = arb(h), EG = e^gamma.
    Возвращает arb либо None, если узел меньше exp(-tol)."""
    tf = i * plan.h
    lnf, xf, uf, _ = _est(tf)
    need = plan.tol + lnf                     # нужная относительная точность узла, нат
    if need < 0.0:
        return None
    P = min(max(int(need / LN2) + NODE_PAD_BITS, 128), plan.prec)

    flint.ctx.prec = P
    x = arb(i) * H - em                       # ln u
    u = x.exp()
    j1 = 1 + em                               # du/dt = u * j1

    if xf < 0.0:                              # ---- u < 1
        K = _series_terms(xf, P + 8)
        if K is not None:                     # F = e^gamma u exp(V(u)), без log и E1
            mu = -u
            p = mu
            acc = mu + p                      # -2u  (= -u + (-u)^1/1!)
            for k in range(2, K + 1):
                p = p * mu / k
                acc = acc + p / k
            return EG * u * u * j1 * acc.exp()
        flint.ctx.prec = P + 32
        e1 = u.expint(1)
        flint.ctx.prec = P
        return (-u - e1).exp() * u * j1

    nb = (plan.tol - 2.0 * uf) / LN2          # ---- u >= 1: E1 на полной точности
    if uf >= E1_ZERO_U:
        return (-u).exp() * u * j1            # E1 ~ e^-u/u < 2^-prec
    K = _asym_terms(uf, nb + 20.0)
    if K is not None:                         # асимптотический ряд на малой точности
        w = (-u).exp()
        flint.ctx.prec = max(64, int(nb) + 48)
        S = arb(1)
        for k in range(K - 1, 0, -1):
            S = 1 - (S * k) / u
        e1 = w * S / u
        flint.ctx.prec = P
        return w * (-e1).exp() * u * j1
    e1 = _exp1_full(u, uf, plan.prec)
    flint.ctx.prec = P
    return (-u - e1).exp() * u * j1


def _node_ref(i, plan):
    """Эталонное значение узла: без сокращений точности, без рядов, E1 только через expint."""
    prec = plan.prec + 96
    flint.ctx.prec = prec
    t = arb(i) * arb(plan.h)
    em = (-t).exp()
    u = (t - em).exp()
    _, xf, uf, _ = _est(i * plan.h)
    if xf >= 0.0:
        e1 = _exp1_full(u, uf, prec)
    else:
        flint.ctx.prec = prec
        e1 = u.expint(1)
    flint.ctx.prec = prec
    return (-u - e1).exp() * u * (1 + em)


# ---------------------------------------------------------------------------
# Сериализация arb (точная) и разбор
# ---------------------------------------------------------------------------

def _parse_num(s):
    s = s.strip().lstrip("[").strip()
    tok = s.split()[0] if s else "0"
    tok = tok.split("+/-")[0]
    return mp.mpf(tok) if tok else mp.mpf(0)


def _dump(x):
    """Середина шара точно: (мантисса, порядок) + радиус (текстом)."""
    m = x.mid()
    try:
        man, ex = m.man_exp()
        rec = {"man": str(int(man)), "exp": int(ex)}
    except Exception:                          # запасной путь: десятичная печать
        rec = {"dec": str(m)}
    try:
        rec["rad"] = str(x.rad())
    except Exception:
        rec["rad"] = "0"
    return rec


def _load(rec, need_digits):
    if "man" in rec:
        return mp.ldexp(mp.mpf(int(rec["man"])), int(rec["exp"]))
    tok = rec["dec"].strip().lstrip("[").split()[0]
    mant = tok.lower().split("e")[0]
    if sum(ch.isdigit() for ch in mant) < need_digits:
        raise RuntimeError(
            f"str(arb) напечатал меньше {need_digits} цифр - точность потеряна; "
            "нужна версия python-flint с arb.man_exp()")
    return mp.mpf(tok)


# ---------------------------------------------------------------------------
# Работа над куском узлов
# ---------------------------------------------------------------------------

_EG_CACHE = {}


def _euler_exp(prec):
    g = _EG_CACHE.get(prec)
    if g is None:
        old = flint.ctx.prec
        flint.ctx.prec = prec
        g = arb.const_euler().exp()
        flint.ctx.prec = old
        _EG_CACHE[prec] = g
    return g


def _chunk(job):
    """Сумма узлов [lo, hi) раздельно по чётным и нечётным индексам."""
    lo, hi, digits, guard, const = job
    plan = Plan(digits, guard, const)
    flint.ctx.prec = plan.prec
    EG = _euler_exp(plan.prec)
    H = arb(plan.h)
    R = (-H).exp()
    s = [arb(0), arb(0)]
    em = None
    for n, i in enumerate(range(lo, hi)):
        flint.ctx.prec = plan.prec
        if em is None or n % RECOMP == 0:
            em = (-(arb(i) * H)).exp()
        else:
            em = em * R
        v = _node(i, plan, H, em, EG)
        if v is not None:
            flint.ctx.prec = plan.prec
            s[i & 1] = s[i & 1] + v
    flint.ctx.prec = plan.prec
    return {"lo": lo, "hi": hi, "e": _dump(s[0]), "o": _dump(s[1])}


def _atomic_write(path, text):
    tmp = path + ".tmp"
    with open(tmp, "w", encoding="utf-8") as fh:
        fh.write(text)
    os.replace(tmp, path)


def run_nodes(plan, jobs=1, chunk=None, ckpt=None, progress=True):
    """Считает все куски (с чекпоинтами).  Возвращает {lo: запись}."""
    n0, n1 = -plan.k_lo, plan.k_hi + 1
    total = n1 - n0
    if chunk is None:
        chunk = max(50, min(1500, total // (max(1, jobs) * 24)))
    bounds = [(lo, min(lo + chunk, n1)) for lo in range(n0, n1, chunk)]

    done = {}
    if ckpt:
        os.makedirs(ckpt, exist_ok=True)
        mpath = os.path.join(ckpt, "meta.json")
        meta = plan.meta(chunk)
        if os.path.exists(mpath):
            with open(mpath, encoding="utf-8") as fh:
                if json.load(fh) != meta:
                    raise SystemExit(f"{ckpt}: чекпоинты другого прогона (другие digits/guard/const/chunk)")
        else:
            _atomic_write(mpath, json.dumps(meta))
        for lo, hi in bounds:
            p = os.path.join(ckpt, f"c_{lo}.json")
            if os.path.exists(p):
                try:
                    with open(p, encoding="utf-8") as fh:
                        rec = json.load(fh)
                    if rec["lo"] == lo and rec["hi"] == hi:
                        done[lo] = rec
                except (ValueError, KeyError):
                    pass

    todo = sorted((b for b in bounds if b[0] not in done), reverse=True)   # дорогие (большие u) - первыми
    t0 = time.time()
    finished = 0

    def store(rec):
        nonlocal finished
        done[rec["lo"]] = rec
        if ckpt:
            _atomic_write(os.path.join(ckpt, f"c_{rec['lo']}.json"), json.dumps(rec))
        finished += 1
        if progress:
            el = time.time() - t0
            eta = el / finished * (len(todo) - finished)
            print(f"\r  куски {len(done)}/{len(bounds)}  (новых {finished}/{len(todo)})  "
                  f"прошло {el:8.0f}с  осталось ~{eta:8.0f}с   ", end="", file=sys.stderr, flush=True)

    if progress and done:
        print(f"  из чекпоинта взято кусков: {len(done)}/{len(bounds)}", file=sys.stderr, flush=True)
    if jobs <= 1:
        for lo, hi in todo:
            store(_chunk((lo, hi, plan.digits, plan.guard, plan.const)))
    else:
        with cf.ProcessPoolExecutor(max_workers=jobs) as ex:
            futs = [ex.submit(_chunk, (lo, hi, plan.digits, plan.guard, plan.const)) for lo, hi in todo]
            for fut in cf.as_completed(futs):
                store(fut.result())
    if progress and todo:
        print(file=sys.stderr, flush=True)
    return done


def combine(plan, done):
    """Собирает частичные суммы; возвращает (lambda, Delta, радиус) как mpf."""
    mp.mp.prec = plan.prec + 256
    need = plan.digits + plan.guard
    Se = mp.mpf(0)
    So = mp.mpf(0)
    rad = mp.mpf(0)
    for lo in sorted(done):
        rec = done[lo]
        Se += _load(rec["e"], need)
        So += _load(rec["o"], need)
        rad += _parse_num(rec["e"]["rad"]) + _parse_num(rec["o"]["rad"])
    H = mp.mpf(plan.h)
    return H * (Se + So), 2 * H * abs(Se - So), H * rad


def _lg(x):
    return float(mp.log10(x)) if x > 0 else -math.inf


def to_decimal(x, digits):
    mp.mp.dps = digits + 30
    return mp.nstr(x, digits, strip_zeros=False)


# ---------------------------------------------------------------------------
# Аудит узлов
# ---------------------------------------------------------------------------

def audit(plan, n):
    """Сверяет рабочий путь с эталонным на n узлах (равномерно + на границах режимов)."""
    n0, n1 = -plan.k_lo, plan.k_hi + 1
    idx = {n0, n1 - 1}
    for k in range(n):
        idx.add(n0 + (n1 - 1 - n0) * k // max(1, n - 1))
    for u in (1e-300, 1e-30, 1e-3, 0.5, 1.0, 2.0, 50.0, 1e3, plan.tol / 6, plan.tol / 3,
              0.45 * plan.tol, 0.55 * plan.tol, plan.tol):
        i = int(round(_t_of_u(u) / plan.h))
        if n0 <= i < n1:
            idx.add(i)
    mp.mp.prec = plan.prec + 512
    tol_abs = mp.exp(-(plan.tol + 10.0))
    flint.ctx.prec = plan.prec
    EG = _euler_exp(plan.prec)
    H = arb(plan.h)
    bad = 0
    print(f"аудит: {len(idx)} узлов, допуск на узел 10^{_lg(tol_abs):.0f}", file=sys.stderr, flush=True)
    for i in sorted(idx):
        flint.ctx.prec = plan.prec
        em = (-(arb(i) * H)).exp()
        a = _node(i, plan, H, em, EG)
        b = _node_ref(i, plan)
        if a is None:
            print(f"  i={i:9d} пропущен (меньше допуска)", file=sys.stderr, flush=True)
            continue
        da = _load(_dump(a), 0)
        db = _load(_dump(b), 0)
        d = abs(da - db)
        ok = d <= tol_abs
        bad += not ok
        xf = _est(i * plan.h)[1]
        print(f"  i={i:9d} ln u={xf:11.2f}  |fast-ref| = 10^{_lg(d):9.1f}  {'ok' if ok else 'ОШИБКА'}",
              file=sys.stderr, flush=True)
    return bad


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(argv=None):
    ap = argparse.ArgumentParser(description="Постоянная Голомба-Дикмана (lambda), Takahasi-Mori + flint")
    ap.add_argument("digits", type=int, help="число знаков после запятой")
    ap.add_argument("--guard", type=int, default=30, help="запас знаков (по умолч. 30)")
    ap.add_argument("--const", type=float, default=None,
                    help=f"шаг h = const/T (по умолч. {DEFAULT_CONST}); больше 9 - нельзя")
    ap.add_argument("--jobs", type=int, default=1, help="число процессов")
    ap.add_argument("--chunk", type=int, default=None, help="узлов в куске (по умолч. автоматически)")
    ap.add_argument("--ckpt", type=str, default=None, help="каталог чекпоинтов (продолжение прерванного прогона)")
    ap.add_argument("--out", type=str, default=None, help="файл для результата")
    ap.add_argument("--audit", type=int, default=0, help="сверить N узлов с эталонным путём перед прогоном")
    ap.add_argument("--plan", action="store_true", help="только показать параметры и выйти")
    ap.add_argument("--check", action="store_true", help="сверить с 100 знаками OEIS (и аудит 12 узлов)")
    args = ap.parse_args(argv)

    const = DEFAULT_CONST if args.const is None else args.const
    plan = Plan(args.digits, args.guard, const)
    jobs = max(1, min(args.jobs, os.cpu_count() or 1))

    n_left = plan.k_lo
    print(f"lambda: {plan.digits} знаков (+{plan.guard} запас), узлов {plan.n_nodes} "
          f"(u<1: {n_left}, u>1: {plan.k_hi + 1}), h={plan.h:.4e} (const={plan.const}), "
          f"полная точность {plan.prec} бит, процессов {jobs}", file=sys.stderr, flush=True)
    if args.plan:
        return 0

    n_audit = args.audit if args.audit else (12 if args.check else 0)
    if n_audit:
        if audit(plan, n_audit):
            print("аудит НЕ пройден - результат считать нельзя", file=sys.stderr)
            return 3

    t0 = time.time()
    done = run_nodes(plan, jobs=jobs, chunk=args.chunk, ckpt=args.ckpt)
    val, delta, rad = combine(plan, done)
    dec = to_decimal(val, plan.digits)
    dt = time.time() - t0

    q_est = 2.0 * _lg(delta)        # log10 оценки ошибки квадратуры  (~Delta^2)
    q_rad = _lg(rad)                # log10 суммарного радиуса шаров
    need = -(plan.digits + 3)
    ok = q_est < need and q_rad < need
    print(f"готово за {dt:.0f}с;  оценка ошибки квадратуры ~10^{q_est:.0f},  радиус шаров ~10^{q_rad:.0f}  "
          f"(нужно < 10^{need})  -> {'OK' if ok else 'НЕДОСТАТОЧНО'}", file=sys.stderr, flush=True)

    if args.check and plan.digits >= 100:
        n = 99
        same = dec[:n] == KNOWN_100[:n]
        print(f"сверка с OEIS A084945 (первые {n} символов): {'OK' if same else 'НЕ СОВПАЛО'}",
              file=sys.stderr, flush=True)
        ok = ok and same

    header = (f"# Golomb-Dickman constant lambda, {plan.digits} digits; Takahasi-Mori, "
              f"nodes={plan.n_nodes}, const={plan.const}; quad.err.est=10^{q_est:.0f}, "
              f"ball.radius=10^{q_rad:.0f}{'' if ok else '  *** НЕ ПОДТВЕРЖДЕНО ***'}\n")
    if args.out:
        with open(args.out, "w", encoding="utf-8") as fh:
            fh.write(header + dec + "\n")
        print(f"записано {args.out}", file=sys.stderr, flush=True)
    else:
        sys.stdout.write(header + dec + "\n")
    return 0 if ok else 2


if __name__ == "__main__":
    sys.exit(main())
