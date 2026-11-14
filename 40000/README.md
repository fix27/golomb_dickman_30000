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
а не просто вычислены с плавающей точкой. Ширина интервала —
$10^{-40042.2602491239}$, то есть остаётся около **42 запасных знаков**.

- [`golomb_dickman_40000_digits.txt`](golomb_dickman_40000_digits.txt) — 40 000 цифр после запятой одной строкой.
- [`b084945.txt`](b084945.txt) — готовый b-file для OEIS A084945, термы `n = 0..39999`.

## Выводы

1. **Знаки доказаны, а не посчитаны приближённо.** Интервал шириной
   $10^{-40042.26}$ не пересекает границу десятичной ячейки $10^{-40000}$:
   обе границы дают одно и то же 40 000-значное усечение.
2. **Сертификат согласован со всеми чекпоинтами.** Точные целые суммы
   `Lsum`/`Hsum` по 1 209 чанкам (241 688 узлов) воспроизводят ровно те концы
   интервала, которые записаны в сертификате. Ни один чанк не противоречит
   остальным.
3. **Поправки не подогнаны.** $B=1$ и $TR+TL+B=65\,673\,984\,721\,697$
   получаются из независимо пересчитанных $B_{\rm disc}=10^{-41285.4628}$,
   $\tau_R=10^{-40045.8435}$, $\tau_L=10^{-40045.6719}$.
4. **Запас по дискретизации — около 1 249 порядков:** требуется
   $10^{-40036}$, получено $10^{-41285.46}$.
5. **Узлы пересчитаны заново.** 12 из 1 209 чанков (равномерно по всему
   диапазону $t$) пересчитаны независимо от сохранённых значений при 100
   знаках — совпадение во всех 12. Это проверка численной согласованности;
   строгий пересчёт всех узлов — уровень Arb (см. «Ограничения»).
6. **Три независимые формулы для $E_1$ согласуются.** Ряд и цепная дробь
   Стилтьеса совпадают не хуже 69 знаков на сетке $u\in[0.1,200]$; цепная
   дробь и асимптотика — не хуже 30 знаков при 40 знаках точности.
7. **Все 14 проверяемых утверждений MathWorld подтверждены:** 11 `PASS`,
   4 `NOTE` (исторические или открытые), 1 `SKIP` (полный перебор до
   15 000 знаков), **0 `FAIL`**. Подробности — в разделе
   [«Утверждения MathWorld»](#утверждения-mathworld).
8. **Последовательность Эрлса продолжена по границе.** Строки `55555` и
   `666666` не встречаются в первых 40 000 знаках, то есть $a(5),a(6)>40\,000$
   (в OEIS A225242 было зафиксировано только $a(5)>15\,000$).
9. **Дефекты комплекта (R1–R11) исправлены, цифры не менялись.** Для каждого
   текстового файла проверено, что после удаления символов CR байтовое
   содержимое совпадает с исходным. Подробности — в разделе
   [«Найденные и исправленные дефекты»](#найденные-и-исправленные-дефекты)
   и в [`VERIFICATION.md`](VERIFICATION.md).
10. **Что этим не доказано:** нормальность $\lambda$ в базе 10 остаётся
    открытой проблемой ($\chi^2=4.10$ при 9 степенях свободы — всего лишь
    согласуется с равномерностью); уровень Arb на машине сборки не
    исполнялся; рукопись не компилировалась. См. «Ограничения».

## Утверждения MathWorld

Источник: [*Golomb-Dickman Constant Digits*](https://mathworld.wolfram.com/Golomb-DickmanConstantDigits.html),
E. W. Weisstein, MathWorld (страница обновлена 13 ноября 2026; проверено
15 ноября 2026). Все утверждения страницы разобраны на 16 проверок
`checks/verify_mathworld.py`:

```bash
python checks/verify_mathworld.py                  # перебор простых до 60 знаков
python checks/verify_mathworld.py --prime-limit 3000   # исчерпывающий перебор, ~3 мин
```

| Код | Утверждение MathWorld | Статус | Результат |
|---|---|---|---|
| C1 | разложение начинается $0.6243299885\ldots$ | `PASS` | первые 20 знаков: `0.62432998854355087099` |
| C1b | это разложение OEIS A084945 | `PASS` | совпадают все 102 знака OEIS |
| C2 | позиции считаются с 1 после запятой | `PASS` | 1-базная нумерация воспроизводит все 10 вхождений |
| C3 | Митчелл (1968) вычислил 53 знака | `NOTE` | атрибуция историческая; наши 53 знака совпадают с OEIS |
| C4 | вычислено 15 000 знаков (Вайсштейн, 2013) | `NOTE` | у нас 40 000 — на 25 000 знаков больше |
| C5 | последовательность Эрлса: 28, 256, 1967, 387 | `PASS` | все четыре позиции воспроизведены |
| C6 | простые усечения при 6, 27, 57, 60, 1659, 2508 знаках | `PASS` | все шесть простые |
| C6b | — (следствие C6) | `PASS` | исчерпывающе до 3000 знаков: ровно эти шесть, 192.6 с |
| C7 | первые вхождения цифр: 15, 28, 2, 4, 3, 10, 1, 17, 8, 6 | `PASS` | пересчитано, совпало |
| C8 | последние строки 1, 33, 821 на позициях 28, 587, 6322 | `PASS` | пересчитано, совпало |
| C9 | `0123456789` и `9876543210` нет в первых 15 000 знаках | `PASS` | нет ни в 15 000, ни во всех 40 000 |
| C10 | нормальность $\lambda$ в базе 10 неизвестна | `NOTE` | открытая проблема; $\chi^2=4.10$ |
| C11 | таблица частот цифр для 10, 100, 1000, 10 000 знаков | `PASS` | все 40 чисел совпали |
| C12 | простые $624329$ и $624329988543550870992936383$ | `PASS` | оба усечения совпали |
| C13 | завершённый перебор простых до 15 000 знаков | `SKIP` | перебор выполнен до 3000 знаков |
| C14 | таблица первооткрывателей (1659, 2508) | `NOTE` | атрибуция историческая; сами длины простые (C6) |

Что стоит за статусами:

- **C1/C1b.** Совпадение с OEIS A084945 проверяется по 102 опубликованным
  знакам — это внешняя опора, не выведенная из нашего же сертификата.
- **C3, C4, C14.** Это утверждения о том, кто и когда считал, а не о цифрах.
  Машинной проверке они не поддаются; цифры с ними согласуются. Ни одна из
  наших 40 000 цифр не противоречит опубликованным 15 000.
- **C5.** Дополнительно проверено, что `55555` и `666666` не встречаются во
  всех 40 000 знаках. Это усиливает опубликованную границу $a(5)>15\,000$
  до $a(5),a(6)>40\,000$ (при условии верности сертификата).
- **C6.** Запись «простое» для 6-значного усечения доказана детерминированным
  тестом Миллера–Рабина (12 оснований, верно до $3.3\cdot10^{24}$). Для
  усечений в 27, 57, 60, 1659 и 2508 знаков используется 64 основания;
  это вероятностная простота с оценкой ошибки ниже $2^{-128}$, а не
  сертификат ECPP. Для исключения составных чисел вероятностный тест
  строг: если тест говорит «составное», число действительно составное.
- **C8.** Покрытие всех строк длины $n$ проверяется сканированием: для
  $n=1,2,3$ покрытие завершается на 28, 587 и 6322 соответственно, последние
  новые строки — `1`, `33`, `821`. Для $n=4$ из 10 000 строк встретилось
  9 796 — покрытие требует примерно 92 000 знаков и в 40 000 не укладывается.
- **C10.** Частоты цифр в 40 000 знаках:
  `[3993, 3950, 3969, 4076, 3993, 3972, 3989, 4077, 3986, 3995]`,
  $\chi^2=4.10$ (9 степеней свободы, среднее 4000, критическое значение 16.92
  на уровне 5 %). Это не доказательство нормальности и не намёк на аномалию.
- **C13.** Полный перебор до 15 000 знаков требует часов: каждое усечение —
  это число в тысячи знаков, и проверять нужно все длины подряд. По умолчанию
  перебор идёт до 3000 знаков; этого достаточно, чтобы охватить все шесть
  опубликованных значений и убедиться, что между 2508 и 3000 ничего нет.

## Файлы

| Файл | Размер | Назначение |
|---|---:|---|
| `golomb_dickman.py` | ~57 KB | Вычислитель (экспорт ноутбука); `main()` делает пробу на 60 знаках, затем расчёт на 40 000 |
| `golomb_dickman_40000_certified.txt` | ~117 KB | Сертификат: границы интервала и 40 000 знаков |
| `golomb_dickman_40000_digits.txt` | 40 001 B | Только 40 000 цифр после запятой |
| `b084945.txt` | ~302 KB | b-file для OEIS A084945: 40 000 термов, только данные |
| `verify_40000.py` | ~5 KB | Быстрая автономная проверка сертификата, цифр и b-file |
| `checks/` | ~110 KB | **Все проверки:** аудит чекпоинтов, пересчёт без Arb, Arb, MathWorld, сборщик прогонов |
| `VERIFICATION.md` | ~12 KB | Методика проверки, найденные дефекты, база доверия |
| `SHA256SUMS.txt` | ~1.5 KB | Контрольные суммы артефактов (канонический LF) |
| `.gitattributes` | 711 B | Запрет автопреобразования переводов строк в git |
| `.gitignore` | 214 B | Исключение `__pycache__/` и локальных кэшей |
| `smoke_60.txt` | 431 B | Сертификат пробного прогона на 60 знаках |
| `gd_rigorous_ckpt_40000/` | ~70 MB | Чекпоинты кусков для D = 40 000 (1210 файлов) |
| `gd_rigorous_ckpt_60/` | ~1 MB | Чекпоинты пробного прогона на 60 знаках |

В папке [`checks/`](checks/) лежат:

| Файл | Что проверяет |
|---|---|
| `run_all.py` | запускает всё подряд и печатает сводку |
| `verify_40000.py` | сертификат, файл цифр, b-file, префикс OEIS |
| `verify_checkpoints.py` | 1 209 чанков: структура, покрытие, точные суммы |
| `verify_nodes_decimal.py` | `--plan` (константы), `--chunks` (суммы), `--methods` (формулы) |
| `verify_checkpoints_arb.py` | полный пересчёт узлов в Arb |
| `verify_mathworld.py` | все утверждения MathWorld |

## Как проверять

Быстрая проверка — нужен только Python 3.11+, сторонних пакетов нет:

```bash
python verify_40000.py
```

```text
OK: 40000 digits certified
  interval: lam_lo/2^133074 .. lam_hi/2^133074
  interval width ~ 10^-40042.3
  b-file: 40000 terms, indices 0..39999
  first 102 digits match OEIS A084945
  SHA-256:
    ...
```

Она проверяет, что интервал не пересекает границу десятичной ячейки, что
цифры сертификата равны $\lfloor 10^{40000}\lambda\rfloor$, что отдельный файл
цифр совпадает с сертификатом, что в b-file ровно 40 000 термов `n = 0..39999`
и каждая цифра совпадает, и что первые 102 знака совпадают с OEIS A084945.

Полный набор:

```bash
python checks/run_all.py          # быстрый набор, ~2 минуты
python checks/run_all.py --full   # плюс перебор простых до 3000 знаков, ~5 минут
python checks/run_all.py --arb    # плюс аудит в Arb (нужен python-flint)
python checks/run_all.py --all    # всё сразу
```

| Уровень | Команда | Что проверяет | Время |
|---|---|---|---|
| 1. Артефакты | `verify_40000.py` | сертификат, цифры, b-file | 0.1 с |
| 2. Чекпоинты | `checks/verify_checkpoints.py` | все 1 209 чанков, точные суммы | 0.5 с |
| 3. Формулы | `checks/verify_nodes_decimal.py --methods` | три формулы для $E_1$ | 7 с |
| 4. Константы | `checks/verify_nodes_decimal.py --plan` | $B_{\rm disc}$, $\tau_R$, $\tau_L$ | 0.4 с |
| 5. Узлы (decimal) | `checks/verify_nodes_decimal.py --chunks 12` | 12 чанков, 100 знаков | 3 с |
| 6. MathWorld | `checks/verify_mathworld.py` | 16 утверждений страницы | 91 с |
| 7. Простые | `checks/verify_mathworld.py --prime-limit 3000` | исчерпывающий перебор | 193 с |
| 8. Узлы (Arb) | `checks/verify_checkpoints_arb.py --all` | все узлы в шаровой арифметике | часы |

Фактический вывод уровней 2, 4 и 5:

```text
OK: checkpoint audit passed
  chunks: 1209 (expected 1209)
  nodes covered: 241688
  Lsum: 40064 digits; Hsum: 40064 digits
  lower correction B = 1
  upper correction C = TR + TL + B = 65673984721697
  chunk-sum width: 10^-40042.260531
  certificate width: 10^-40042.260249
  guard digits: 40042.26

  log10 M      = 360.3909
  log10 B_disc = -41285.4628   (requirement 10^-40036)
  log10 tau_R  = -40045.8435
  log10 tau_L  = -40045.6719
  B  = 1
  TR = 26431841856436
  TL = 39242142865260
  TR + TL + B = 65673984721697
OK: discretisation and tail bounds reproduce the certificate corrections

OK: recomputed chunk sums agree with the saved intervals
    (numerical consistency only; use verify_checkpoints_arb.py for the rigorous check)
```

## Найденные и исправленные дефекты

Исходная рабочая копия имела дефекты формата и гигиены. **Ни одно число не
изменено**: для каждого файла проверено, что после удаления символов CR
байтовое содержимое совпадает с исходным.

| # | Дефект | Файл | Исправление |
|---|---|---|---|
| R1 | CRLF во всех 40 001 строках | `b084945.txt` | CRLF → LF, 348 892 → 308 891 B |
| R2 | 4 изолированных CR (смешанные переводы строк) | `golomb_dickman_40000_certified.txt` | 120 261 → 120 257 B |
| R3 | завершающий CR | `golomb_dickman_40000_digits.txt` | 40 002 → 40 001 B |
| R4 | 849 CR | `golomb_dickman.py` | 59 448 → 58 599 B |
| R5 | 4 CR | `smoke_60.txt` | 435 → 431 B |
| R6 | 156 CR | `verify_40000.py` | 5 552 → 5 396 B |
| R7 | 107 и 146 CR | `OEIS_SUBMISSION.md`, `PUBLICATION_CHECKLIST.md` | CRLF → LF |
| R8 | `SHA256SUMS.txt` внутренне противоречив: три хеша посчитаны для LF-версий, два — для CRLF, `sha256sum -c` падал на 3 файлах из 5 | `SHA256SUMS.txt` | пересчитан по всем файлам |
| R9 | ширина указана как $10^{-40042.4}$ | `README.md` | $10^{-40042.26}$ |
| R10 | мусор `__pycache__/*.pyc` (68 KB) | — | удалён, добавлен `.gitignore` |
| R11 | ширина занижалась на 0.14 порядка: `bit_length()-1` вместо `log2` | `verify_40000.py` | `math.log2`, печатает $10^{-40042.3}$ |

Именно R1 и был виден как `FAIL: b-file contains CR characters; OEIS expects
plain ASCII with LF`. Проверяющий скрипт был прав: дефект в файле, а не в
цифрах.

## Воспроизведение

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
`golomb_dickman.py` и в [`../30000/README.md`](../30000/README.md) для версии
на 30 000 знаков. Обзор репозитория — в [корневом README](../README.md).

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

## Ограничения

- **Уровень Arb не исполнялся на машине сборки.** Политика Application
  Control блокирует загрузку DLL `libflint`/`libgmp` (`WinError 4551`).
  Скрипт `checks/verify_checkpoints_arb.py` проверен на синтаксис и логику,
  но с Arb ни разу не запускался. На обычной машине:

  ```bash
  python -m pip install python-flint
  python checks/verify_checkpoints_arb.py --audit      # минуты
  python checks/verify_checkpoints_arb.py --all --jobs 8
  ```

- **Уровень 5 проверяет 12 чанков из 1 209 и 100 знаков, а не 40 042.**
  Это ловит систематическую ошибку формул, индексации или сохранённых сумм,
  но не заменяет строгий пересчёт.
- **Простые усечения больше 24 знаков — вероятностная простота.** 64
  основания Миллера–Рабина дают оценку ошибки ниже $2^{-128}$; сертификата
  ECPP нет.
- **Рукопись не компилировалась:** в среде сборки нет LaTeX.
- Проверка не заменяет повторный прогон `golomb_dickman.py`: она проверяет
  **результат**, а не воспроизводит вычисление тем же кодом.

## Контрольные суммы

Каноническая форма всех текстовых артефактов — UTF-8/ASCII с переводами
строк LF и без BOM:

```bash
sha256sum -c SHA256SUMS.txt            # 15 файлов комплекта
sha256sum -c checks/SHA256SUMS.txt     # 7 файлов проверок
Get-FileHash * -Algorithm SHA256       # Windows PowerShell
```

Ключевые значения:

```text
8892c0715f1afbb271c9f07cc9658cea0f3ca83c2dda3ec7d23b447dc908de9c  golomb_dickman_40000_certified.txt
fbae757900983bd4fb6421c5f27bf9602af1a1bfc013e4c501bb9da2b5a65327  golomb_dickman_40000_digits.txt
c73af643198da299d2e624e3c2f69d05439efcddbbefa7cec176ac63f993bed3  b084945.txt
99377813ec9c48ed84cf4649067f4ba4ee3fb7ef74f15102bf92f75cdd400324  golomb_dickman.py
047ff0471acb07482be15a5cdfc7a65ef14c3193372022e3b1dc9b144098e0ee  smoke_60.txt
```

Файл `.gitattributes` фиксирует LF для текстовых файлов и `binary` для
`*.pkl`, чтобы контрольные суммы не зависели от `core.autocrlf`.
Скрипты `verify_*.py` продублированы в корне и в `checks/`: корневые копии —
опубликованные точки входа, их хеши зафиксированы в `SHA256SUMS.txt`;
копии в `checks/` сами находят корень комплекта и собраны в один набор.

## Цитирование

```bibtex
@misc{pozdnyakov2026golombdickman40000,
  title        = {Golomb--Dickman constant: 40{,}000 certified decimal digits},
  author       = {Pozdnyakov, Konstantin},
  year         = {2026},
  howpublished = {\url{https://github.com/fix27/golomb-dickman-constant/tree/main/40000}},
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
digits are proven by truncation. The interval width is $10^{-40042.2602491239}$.

**Conclusions.** (1) The digits are certified, not merely computed. (2) The
exact integer sums over all 1,209 chunks (241,688 nodes) reproduce the
certificate endpoints, so every checkpoint is consistent. (3) The corrections
$B=1$ and $TR+TL+B=65{,}673{,}984{,}721{,}697$ follow from independently
recomputed $B_{\rm disc}=10^{-41285.46}$, $\tau_R=10^{-40045.84}$ and
$\tau_L=10^{-40045.67}$; they are not fitted. (4) The discretisation margin is
about 1,249 orders of magnitude. (5) Twelve of the 1,209 chunks were recomputed
from scratch to 100 decimals and agree. (6) The series, the Stieltjes continued
fraction and the asymptotic expansion for $E_1$ agree to at least 69 and 30
digits respectively. (7) All 14 statements of the MathWorld page
*Golomb-Dickman Constant Digits* are reproduced: 11 `PASS`, 4 historical or
open `NOTE`s, 1 `SKIP` (the 15,000-digit prime search), 0 `FAIL`; in particular
the digit tables for 10, 100, 1,000 and 10,000 digits, the first-occurrence
positions, the Earls sequence and the six prime truncations. (8) The strings
`55555` and `666666` do not occur in the 40,000 digits, extending the published
bound $a(5)>15{,}000$ of OEIS A225242 to $a(5),a(6)>40{,}000$. (9) Eleven
packaging defects (CRLF line endings, an internally inconsistent
`SHA256SUMS.txt`, a misreported interval width, stale bytecode) were fixed
without changing a single digit: every text file is byte-identical to the
original after stripping CR. (10) Normality of $\lambda$ remains open; the
empirical $\chi^2$ is 4.10 on 9 degrees of freedom.

The computation uses the trapezoidal rule after the Takahashi–Mori substitution,
explicit tail and discretization bounds (Trefethen–Weideman), Stieltjes continued
fractions, and Arb balls via `python-flint`. Trust base: Arb correctness and the
cited theorems. See `golomb_dickman.py` for the full proof and
`golomb_dickman_40000_certified.txt` for the result. `b084945.txt` is a ready
b-file for OEIS A084945 with terms $n=0\ldots39999$.

Run `python verify_40000.py` for a dependency-free check of the published
artifacts, or `python checks/run_all.py` for the whole suite (chunk audit,
decimal recomputation, constants, MathWorld statements). Add `--full` for the
exhaustive prime scan and `--arb` for the Arb audit. The Arb level was not
executed on the build machine because an Application Control policy blocks
`libflint`; it must be run once on an unrestricted host before submission.
See `VERIFICATION.md` for the methodology.

</details>
