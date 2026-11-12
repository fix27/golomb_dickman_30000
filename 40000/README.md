# Постоянная Голомба–Дикмана: 40 000 доказанных знаков

$$\lambda=\int_0^\infty F(u)\,du,\qquad F(u)=e^{-u-E_1(u)}
=0.62432998854355087099\ldots$$

Здесь $E_1$ — интегральная показательная функция.

## Результат

В файле [`golomb_dickman_40000_certified.txt`](golomb_dickman_40000_certified.txt)
лежит целочисленный интервал $[\Lambda_-,\Lambda_+]\ni\lambda$, для которого

$$\left\lfloor 10^{40000}\Lambda_-\right\rfloor
=\left\lfloor 10^{40000}\Lambda_+\right\rfloor .$$

Значит, первые **40 000** десятичных знаков доказаны (усечение, без округления),
а не просто вычислены с плавающей точкой. Ширина интервала — примерно
$10^{-40042.4}$, то есть остаётся около 42 запасных знаков.

- [`golomb_dickman_40000_digits.txt`](golomb_dickman_40000_digits.txt) — 40 000 цифр после запятой одной строкой.
- [`b084945.txt`](b084945.txt) — готовый b-file для OEIS A084945, термы `n = 0..39999`.

## Файлы

| Файл | Размер | Назначение |
|---|---:|---|
| `golomb_dickman.py` | ~57 KB | Вычислитель (экспорт ноутбука); `main()` делает пробу на 60 знаках, затем расчёт на 40 000 |
| `golomb_dickman_40000_certified.txt` | ~117 KB | Сертификат: границы интервала и 40 000 знаков |
| `golomb_dickman_40000_digits.txt` | 40 001 B | Только 40 000 цифр после запятой |
| `b084945.txt` | ~302 KB | b-file для OEIS A084945: 40 000 термов, только данные |
| `verify_40000.py` | ~5 KB | Автономная проверка сертификата, цифр и b-file |
| `smoke_60.txt` | 435 B | Сертификат пробного прогона на 60 знаках |
| `gd_rigorous_ckpt_40000/` | ~70 MB | Чекпоинты кусков для D = 40 000 (1210 файлов) |

## Быстрая проверка

Нужен только Python 3.11+; сторонние пакеты не требуются.

```bash
python verify_40000.py
```

Ожидаемый вывод:

```text
OK: 40000 digits certified
  interval: lam_lo/2^133074 .. lam_hi/2^133074
  interval width ~ 10^-40042.4
  b-file: 40000 terms, indices 0..39999
  first 102 digits match OEIS A084945
  SHA-256:
    ...
```

Скрипт проверяет:

1. что сертифицированный интервал не пересекает границу десятичной ячейки;
2. что цифры в сертификате равны $\lfloor 10^{40000}\lambda\rfloor$;
3. что отдельный файл цифр совпадает с сертификатом;
4. что b-file содержит ровно 40 000 термов `n = 0..39999` и каждая цифра совпадает;
5. что первые 102 цифры совпадают с OEIS A084945.

## Как воспроизвести

### Требования

- Python 3.11+;
- [`python-flint`](https://github.com/flintlib/python-flint) (Arb);
- желательно многоядерный CPU; расчёт распараллеливается по процессам.

```bash
python -m pip install python-flint
```

### Полный прогон

```bash
python golomb_dickman.py
```

`main()` сначала выполняет пробный прогон на 60 знаках, затем основной расчёт
на 40 000 знаков. Основной прогон занимает **часы CPU**; его можно прерывать —
прогресс сохраняется в `gd_rigorous_ckpt_40000/`. Чекпоинты привязаны к набору
параметров `(D, G, h, Ebits, K_L, K_R, chunk)`, поэтому при смене `G` или `D`
нужен новый каталог.

### База доверия

Доказательство условно только по отношению к явно перечисленной базе доверия:

1. корректность интервальной арифметики [Arb](https://arblib.org/) (`python-flint`);
2. теорема о трапециях для аналитических функций в полосе
   (Trefethen & Weideman, *SIAM Review* 56 (2014), Thm 5.1);
3. классическая теорема Стилтьеса о двусторонних оценках подходящими дробями;
4. леммы из docstring `golomb_dickman.py` (§1–§7);
5. точная целочисленная арифметика Python.

Это не проверка в proof assistant. Подробности метода — в docstring
`golomb_dickman.py` и в [корневом README](../README.md) для версии на 30 000 знаков.

## B-file для OEIS

Файл `b084945.txt`:

- 40 000 термов, `n = 0..39999`;
- только данные, без комментариев;
- ASCII, LF, без BOM и табов;
- завершающая пустая строка;
- 40 001 строка, 308 891 байт;
- SHA-256: `c73af643198da299d2e624e3c2f69d05439efcddbbefa7cec176ac63f993bed3`.

Текущий b-file в OEIS A084945 содержит только `n = 0..101` (синтезирован из
записи), поэтому этот файл его заменяет. Текст ссылки:

```text
YOUR NAME, <a href="b084945.txt">Table of n, a(n) for n=0..39999</a>
```

Пошаговая инструкция и письмо для редакторов — в
[`OEIS_SUBMISSION.md`](OEIS_SUBMISSION.md). Публикация в Zenodo/arXiv и анонс
в SeqFan — в [`PUBLICATION_CHECKLIST.md`](PUBLICATION_CHECKLIST.md).

## Цитирование

```bibtex
@misc{pozdnyakov2026golombdickman40000,
  title        = {Golomb--Dickman constant: 40{,}000 certified decimal digits},
  author       = {Pozdnyakov, Konstantin},
  year         = {2026},
  howpublished = {\url{https://github.com/fix27/golomb_dickman_30000/tree/main/40000}},
  note         = {Rigorous interval-arithmetic certificate using Arb;
                  OEIS A084945 b-file with 40{,}000 terms}
}
```

## Лицензия

[MIT](../LICENSE) © 2026 Поздняков Константин.

---

<details>
<summary><b>English summary</b></summary>

This directory contains a rigorous interval-arithmetic certificate for the
first **40,000 decimal digits** of the Golomb–Dickman constant
$\lambda = \int_0^\infty e^{-u-E_1(u)}\,du = 0.6243299885\ldots$.
The certificate gives integer bounds $[\Lambda_-,\Lambda_+]\ni\lambda$ with
$\lfloor 10^{40000}\Lambda_-\rfloor=\lfloor 10^{40000}\Lambda_+\rfloor$, so the
digits are proven by truncation. The interval width is about $10^{-40042.4}$.

The computation uses the trapezoidal rule after the Takahashi–Mori substitution,
explicit tail and discretization bounds (Trefethen–Weideman), Stieltjes continued
fractions, and Arb balls via `python-flint`. Trust base: Arb correctness and the
cited theorems. See `golomb_dickman.py` for the full proof and
`golomb_dickman_40000_certified.txt` for the result. `b084945.txt` is a ready
b-file for OEIS A084945 with terms $n=0\ldots39999$.

Run `python verify_40000.py` for a dependency-free verification.

</details>
