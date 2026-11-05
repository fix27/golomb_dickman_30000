#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Convert golomb_dickman_30000_certified.txt -> golomb_dickman_30000_digits.txt

Запуск из Windows:
    python certified_to_digits.py
или двойным щелчком (если .py ассоциирован с Python).

Опции (все необязательны):
    -i, --input   FILE    входной файл   (по умолчанию golomb_dickman_30000_certified.txt)
    -o, --output  FILE    выходной файл  (по умолчанию golomb_dickman_30000_digits.txt)
    -g, --group   N       разбить цифры на строки по N знаков (0 = одна строка)
    --prefix              добавить "0." перед цифрами
    --expect D            ожидаемая длина (по умолчанию 30000); 0 = не проверять

Формат выходного файла по умолчанию:
    30000 цифр после запятой, одна непрерывная строка, без "0.", без переводов строк.
"""

import argparse
import os
import sys


def read_digits(path):
    """Извлекает цифры из сертификата. Возвращает строку только из 0-9."""
    # Python >= 3.11: поднимаем лимит на str<->int для больших чисел.
    if hasattr(sys, "set_int_max_str_digits"):
        sys.set_int_max_str_digits(200_000)

    digits = None
    with open(path, "r", encoding="utf-8") as f:
        for lineno, raw in enumerate(f, 1):
            s = raw.strip()
            if not s or s.startswith("#"):
                continue
            if s.startswith("0."):
                digits = s[2:].strip()
            elif s.isdigit():
                digits = s
            else:
                raise ValueError(f"{path}:{lineno}: unexpected line: {s[:60]!r}...")
            break

    if digits is None:
        raise ValueError(f"{path}: no data line (starting with '0.' or all digits) found")
    if not digits.isdigit():
        bad = next((i for i, c in enumerate(digits) if not c.isdigit()), None)
        raise ValueError(f"{path}: non-digit character at position {bad}: {digits[bad]!r}")
    return digits


def write_digits(path, digits, group=0, prefix=False):
    if group > 0:
        body = "\n".join(digits[i:i + group] for i in range(0, len(digits), group))
    else:
        body = digits
    if prefix:
        body = "0." + body
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        f.write(body)
        if not body.endswith("\n"):
            f.write("\n")


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Extract digits from golomb_dickman_30000_certified.txt",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    ap.add_argument("-i", "--input", default="golomb_dickman_30000_certified.txt")
    ap.add_argument("-o", "--output", default="golomb_dickman_30000_digits.txt")
    ap.add_argument("-g", "--group", type=int, default=0,
                    help="split into lines of N digits (0 = single line, default)")
    ap.add_argument("--prefix", action="store_true",
                    help="prefix with '0.'")
    ap.add_argument("--expect", type=int, default=30000,
                    help="expected number of digits (0 = don't check; default 30000)")
    args = ap.parse_args(argv)

    if not os.path.isfile(args.input):
        sys.exit(f"error: input file not found: {os.path.abspath(args.input)}")

    digits = read_digits(args.input)
    if args.expect and len(digits) != args.expect:
        sys.exit(f"error: got {len(digits)} digits, expected {args.expect} "
                 f"(pass --expect 0 to skip this check)")

    # независимая сверка с OEIS A084945 (первые 94 знака после запятой)
    KNOWN = ("6243299885435508709929363831008372441796426201805292869735519024"
             "956380888551132544624602761955398688")
    k = min(len(digits), len(KNOWN))
    if digits[:k] != KNOWN[:k]:
        print(f"WARNING: first {k} digits do not match OEIS A084945", file=sys.stderr)

    write_digits(args.output, digits, group=args.group, prefix=args.prefix)
    print(f"OK: {args.output}  ({len(digits)} digits, "
          f"{'one line' if args.group <= 0 else f'grouped by {args.group}'}"
          f"{', prefix 0.' if args.prefix else ''})")


if __name__ == "__main__":
    main()