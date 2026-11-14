# Постоянная Голомба–Дикмана: сертифицированные вычисления

$$\lambda=\int_0^\infty F(u)\,du,\qquad F(u)=e^{-u-E_1(u)}
=0.624329988543550870992936383100837244\ldots$$

Здесь $E_1$ — интегральная показательная функция. В репозитории лежат
вычисления этой константы в интервальной арифметике с явным сертификатом
точности: 30 000 знаков (первая версия) и **40 000 знаков** (текущий результат).

## Текущий результат: 40 000 знаков

Сертификат
[`40000/golomb_dickman_40000_certified.txt`](40000/golomb_dickman_40000_certified.txt)
содержит целочисленный интервал $[\Lambda_-,\Lambda_+]\ni\lambda$ такой, что

$$\left\lfloor 10^{40000}\Lambda_-\right\rfloor
=\left\lfloor 10^{40000}\Lambda_+\right\rfloor .$$

Поэтому первые **40 000** десятичных знаков (усечение, без округления)
доказаны, а не просто вычислены с плавающей точкой. Ширина интервала —
около $10^{-40042.26}$, то есть в запасе остаётся ещё примерно 42 знака.

| Артефакт | Файл |
|---|---|
| Сертификат: границы интервала и 40 000 знаков | [`40000/golomb_dickman_40000_certified.txt`](40000/golomb_dickman_40000_certified.txt) |
| Только цифры, одной строкой | [`40000/golomb_dickman_40000_digits.txt`](40000/golomb_dickman_40000_digits.txt) |
| b-file для OEIS A084945, термы `n = 0..39999` | [`40000/b084945.txt`](40000/b084945.txt) |
| Автономная проверка без сторонних пакетов | [`40000/verify_40000.py`](40000/verify_40000.py) |

Контрольные суммы:

```text
8892c0715f1afbb271c9f07cc9658cea0f3ca83c2dda3ec7d23b447dc908de9c  golomb_dickman_40000_certified.txt
fbae757900983bd4fb6421c5f27bf9602af1a1bfc013e4c501bb9da2b5a65327  golomb_dickman_40000_digits.txt
c73af643198da299d2e624e3c2f69d05439efcddbbefa7cec176ac63f993bed3  b084945.txt
```

Полный список — в [`40000/SHA256SUMS.txt`](40000/SHA256SUMS.txt).

## Быстрая проверка

Нужен только Python 3.11+; сторонние пакеты не требуются.

```bash
git clone https://github.com/fix27/golomb-dickman-constant.git
cd golomb-dickman-constant/40000
python verify_40000.py
```

Ожидаемый вывод:

```text
OK: 40000 digits certified
  interval: lam_lo/2^133074 .. lam_hi/2^133074
  interval width ~ 10^-40042.3
  b-file: 40000 terms, indices 0..39999
  first 102 digits match OEIS A084945
  SHA-256:
    ...
```

Скрипт проверяет, что сертифицированный интервал не пересекает границу
десятичной ячейки, что цифры равны $\lfloor 10^{40000}\lambda\rfloor$, что
файл цифр и b-file согласованы с сертификатом и что первые 102 цифры
совпадают с OEIS A084945. Более глубокий аудит чекпоинтов описан в
[`40000/VERIFICATION.md`](40000/VERIFICATION.md).

## Структура репозитория

| Путь | Назначение |
|---|---|
| [`40000/`](40000/) | Текущий результат: 40 000 знаков, b-file, проверки, документация |
| [`30000/`](30000/) | Первый сертифицированный результат: 30 000 знаков |
| [`30000-v2/`](30000-v2/) | Вторая итерация пайплайна для 30 000 знаков (v23); сертификаты на 2 000 и 30 000 |
| [`alt2-50000/`](alt2-50000/) | Незавершённая попытка на 50 000 знаков: чекпоинты и текущий CLI; сертификата пока нет |
| [`old_version/`](old_version/) | Архив ранних версий программы (v22–v25) и `archive/` |
| [`coq/`](coq/), [`lean/`](lean/) | Экспериментальная формализация отдельных лемм (Coq/mathcomp, Lean/Mathlib) |
| [`golomb_dickman_26.py`](golomb_dickman_26.py) | Текущий самодостаточный CLI (v26) |

В корне также лежат копии трёх артефактов 40 000 знаков: `b084945.txt`,
`golomb_dickman_40000_certified.txt`, `golomb_dickman_40000_digits.txt`.
Это удобные копии в CRLF; **канонические файлы** находятся в `40000/` —
они в LF, и именно для них приведены SHA-256. Для OEIS и проверок
используйте файлы из `40000/`.

## Как воспроизвести

### Требования

- Python 3.11+;
- [`python-flint`](https://github.com/flintlib/python-flint) (Arb) — обязателен;
- [`mpmath`](https://mpmath.org/) — только для `crosscheck`;
- желательно многоядерный CPU: расчёт распараллеливается по процессам.

```bash
python -m pip install python-flint mpmath
```

### Текущий CLI

```bash
python golomb_dickman_26.py test
python golomb_dickman_26.py smoke
python golomb_dickman_26.py certify --digits 40000
python golomb_dickman_26.py crosscheck --digits 200
python golomb_dickman_26.py boundary-test
```

- `test` — самопроверка окружения и точные рациональные проверки леммы 7;
- `smoke` — сквозной прогон на 60 знаках;
- `certify --digits 40000` — полный сертификат (значение по умолчанию для
  `--digits` — 40 000); пишет `golomb_dickman_40000_certified.txt` в текущий
  каталог и сохраняет чекпоинты в `gd_rigorous_ckpt_40000/`;
- `crosscheck` — независимые сверки через `mpmath`/`mpmath.iv`, цепные дроби,
  форсированные ветви, альтернативный шаг и лестницу точности;
- `boundary-test` — граничные тесты леммы 7.

Полный прогон на 40 000 знаков занимает **часы CPU**; его можно прерывать,
прогресс сохраняется в чекпоинтах. Режим `--resume` — только аудиторский:
он доверяет готовым чекпоинтам, не пересчитывает их и **не выдаёт
сертификат**. По умолчанию работает безопасный режим `--no-resume`:
каждый кусок пересчитывается в Arb.

Текущий CLI использует формат чекпоинтов v12 (`meta_v12.json`). Опубликованные
чекпоинты `40000/gd_rigorous_ckpt_40000/` имеют старый формат `meta.pkl` и
относятся к [`40000/golomb_dickman.py`](40000/golomb_dickman.py) — пайплайну,
которым получен опубликованный сертификат. Его можно запустить так:

```bash
cd 40000
python golomb_dickman.py
```

`main()` сначала выполняет пробный прогон на 60 знаках, затем основной расчёт
на 40 000 знаков, используя существующие чекпоинты.

### База доверия

Доказательство условно только по отношению к явно перечисленной базе доверия:

1. корректность интервальной арифметики [Arb](https://arblib.org/)
   (`python-flint`);
2. теорема о трапециях для аналитических функций в полосе
   (Trefethen & Weideman, *SIAM Review* 56 (2014), Thm 5.1);
3. классическая теорема Стилтьеса о двусторонних оценках подходящими дробями;
4. леммы из docstring вычислителя (§1–§7);
5. точная целочисленная арифметика Python.

Это не проверка в proof assistant. Каталоги `coq/` и `lean/` содержат лишь
экспериментальную формализацию отдельных лемм и не заменяют сертификат.

## B-file для OEIS

Файл [`40000/b084945.txt`](40000/b084945.txt):

- 40 000 термов, `n = 0..39999`;
- только данные, без комментариев;
- ASCII, LF, без BOM и табов;
- завершающая пустая строка;
- 40 001 строка, 308 891 байт;
- SHA-256: `c73af643198da299d2e624e3c2f69d05439efcddbbefa7cec176ac63f993bed3`.

На момент подготовки файла b-file в OEIS A084945 содержал только `n = 0..101`
(он синтезирован из записи), поэтому новый файл его заменяет. Текст ссылки:

```text
YOUR NAME, <a href="b084945.txt">Table of n, a(n) for n=0..39999</a>
```

Пошаговая инструкция и письмо для редакторов — в
[`40000/OEIS_SUBMISSION.md`](40000/OEIS_SUBMISSION.md). Публикация в
Zenodo/arXiv и анонс в SeqFan — в
[`40000/PUBLICATION_CHECKLIST.md`](40000/PUBLICATION_CHECKLIST.md).

## Документация

| Документ | Содержание |
|---|---|
| [`40000/README.md`](40000/README.md) | Подробное описание результата на 40 000 знаков |
| [`40000/VERIFICATION.md`](40000/VERIFICATION.md) | Что и как проверяется, найденные дефекты, база доверия |
| [`40000/OEIS_SUBMISSION.md`](40000/OEIS_SUBMISSION.md) | Отправка b-file в OEIS |
| [`40000/PUBLICATION_CHECKLIST.md`](40000/PUBLICATION_CHECKLIST.md) | Zenodo, arXiv, SeqFan |
| [`30000/README.md`](30000/README.md) | Детали результата на 30 000 знаков |

## Цитирование

```bibtex
@misc{pozdnyakov2026golombdickman,
  title        = {Golomb--Dickman constant: 40{,}000 certified decimal digits},
  author       = {Pozdnyakov, Konstantin},
  year         = {2026},
  howpublished = {\url{https://github.com/fix27/golomb-dickman-constant/tree/main/40000}},
  note         = {Rigorous interval-arithmetic certificate using Arb;
                  OEIS A084945 b-file with 40{,}000 terms}
}
```

## Лицензия

[MIT](LICENSE) © 2026 Поздняков Константин.

---

<details>
<summary><b>English summary</b></summary>

This repository contains rigorous interval-arithmetic certificates for decimal
digits of the Golomb–Dickman constant
$\lambda = \int_0^\infty e^{-u-E_1(u)}\,du = 0.6243299885\ldots$.
The current result is **40,000 certified decimal digits**; an earlier 30,000-digit
result is kept in `30000/` and `30000-v2/`.

The certificate in `40000/golomb_dickman_40000_certified.txt` gives integer
bounds $[\Lambda_-,\Lambda_+]\ni\lambda$ with
$\lfloor 10^{40000}\Lambda_-\rfloor=\lfloor 10^{40000}\Lambda_+\rfloor$,
so the digits are proven by truncation. The interval width is about
$10^{-40042.26}$.

The computation uses the trapezoidal rule after the Takahashi–Mori substitution,
explicit tail and discretization bounds (Trefethen–Weideman), Stieltjes continued
fractions, and Arb balls via `python-flint`. Trust base: Arb correctness and the
cited theorems. The current self-contained CLI is `golomb_dickman_26.py`;
`40000/golomb_dickman.py` is the pipeline that produced the published certificate.

`40000/b084945.txt` is a ready OEIS b-file for A084945 with terms
$n=0\ldots39999$ (plain ASCII, LF, 40,001 lines). Run
`python 40000/verify_40000.py` for a dependency-free check of the certificate,
digits and b-file. See `40000/VERIFICATION.md` for the full checkpoint audit.

</details>
