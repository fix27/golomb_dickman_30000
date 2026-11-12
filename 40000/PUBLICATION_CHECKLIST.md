# Публикация результата: Zenodo, arXiv, SeqFan

Этот чек-лист дополняет [`OEIS_SUBMISSION.md`](OEIS_SUBMISSION.md).
OEIS — основной способ зафиксировать расширение b-file, но для цитируемости
и приоритета полезно сделать отдельный DOI и, при желании, препринт.

## 1. Zenodo DOI — самый быстрый цитируемый след

Репозиторий уже публичный:
https://github.com/fix27/golomb_dickman_30000

Порядок:

1. Убедитесь, что в репозитории лежат:
   - `40000/README.md`;
   - `40000/golomb_dickman_40000_certified.txt`;
   - `40000/golomb_dickman_40000_digits.txt`;
   - `40000/b084945.txt`;
   - `40000/verify_40000.py`;
   - `40000/golomb_dickman.py`;
   - `LICENSE`.
2. Сделайте release с тегом, например `v40000` или `v1.0.0-40000`.
3. В Zenodo: **Settings → GitHub → включить репозиторий**.
4. Создайте release в GitHub; Zenodo автоматически выпустит DOI.
5. Добавьте DOI в `40000/README.md` и в описание GitHub.
6. В метаданных Zenodo укажите:
   - Title: `Golomb–Dickman constant: 40,000 certified decimal digits`;
   - Authors: `Pozdnyakov, Konstantin`;
   - Description: 2–3 абзаца о сертификате, методе и b-file;
   - License: MIT;
   - Keywords: `Golomb-Dickman constant`, `interval arithmetic`, `Arb`, `OEIS A084945`;
   - Related identifiers: `https://oeis.org/A084945` (is documented by / is supplement to).

Проверьте, что DOI ведёт на **конкретную версию** (version DOI), а не только
на «all versions». Version DOI удобнее цитировать.

## 2. arXiv — если нужен формальный препринт

Вариант публикации: короткая заметка на 2–5 страниц.

- Раздел: **math.NA** (Numerical Analysis), cross-list **cs.MS** (Mathematical Software).
  math.NT тоже возможен, если делать акцент на самой постоянной, а не на методе.
- Для первого препринта может потребоваться **endorsement** от состоявшегося
  автора в выбранной категории. Это нормальное ограничение arXiv.
- Не называйте работу «мировым рекордом» без независимого подтверждения.
  Формулировка «extends the OEIS b-file for A084945 from 101 to 40,000 terms»
  фактична и проверяема.

Возможная структура:

1. **Introduction** — постоянная, зачем нужны цифры, что уже было в OEIS.
2. **Method** — подстановка Такахаси–Мори, правило трапеций, оценки
   Trefethen–Weideman, цепные дроби Стилтьеса, Arb.
3. **Certificate format** — `S`, `lam_lo`, `lam_hi`, проверка
   $\lfloor 10^{40000}\Lambda_-\rfloor=\lfloor 10^{40000}\Lambda_+\rfloor$.
4. **Verification** — `verify_40000.py`, совпадение первых 102 цифр с A084945,
   независимая проверка сертификата.
5. **Reproducibility** — версии Python/`python-flint`, команда запуска,
   время, чекпоинты, SHA-256.
6. **Limitations** — база доверия (Arb + цитируемые теоремы), не proof assistant.
7. **Data availability** — GitHub, Zenodo DOI, OEIS b-file.

Что не стоит делать: выдавать статистические тесты цифр за доказательство
нормальности. 40 000 цифр — это проверка выборки, а не доказательство.

## 3. SeqFan — анонс и просьба о помощи

SeqFan: https://groups.google.com/g/seqfan

Короткое сообщение (на английском, так принято в списке):

```text
Subject: 40,000 certified digits of the Golomb-Dickman constant (A084945)

I have computed and certified the first 40,000 decimal digits of the
Golomb-Dickman constant using interval arithmetic (Arb). The certificate
gives integer bounds [L-, L+] with
floor(10^40000 L-) = floor(10^40000 L+) = 624329...

Data: https://github.com/fix27/golomb_dickman_30000/tree/main/40000
OEIS: https://oeis.org/A084945

I have prepared a b-file with terms n=0..39999 (currently OEIS has 0..101).
I would be grateful if an editor could help upload it, or advise on the
preferred size (10,000 / 20,000 / 40,000 terms).

Verification: python verify_40000.py (standard library only).
Comments and corrections are welcome.
```

Не дублируйте одно и то же сообщение в нескольких тредах; одного анонса
в SeqFan достаточно.

## 4. MathOverflow / Stack Exchange — осторожно

- **MathOverflow не предназначен для проверки кода или «оцените мою
  реализацию»**. Такие вопросы закрывают.
- Для проверки кода используйте **Code Review Stack Exchange** или
  **Computational Science Stack Exchange**.
- На MathOverflow имеет смысл задавать только математический вопрос,
  например: «Известны ли нетривиальные следствия из большого числа
  сертифицированных цифр λ?» — но и он должен быть исследовательским,
  а не просьбой о рекламе результата.

## 5. Статистика цифр — только как sanity check

Если делать статистический анализ 40 000 цифр:

- частоты цифр 0–9 и критерий хи-квадрат (df = 9);
- частоты пар цифр (df = 99);
- длины серий, число повторений, автокорреляция;
- энтропия.

Но:

- 40 000 цифр мало для сильных утверждений о нормальности;
- при проверке десятков гипотез нужна поправка на множественные сравнения;
- p-value — не доказательство нормальности;
- результат стоит подавать как «эмпирические наблюдения», а не как «доказано».

## 6. Чистота репозитория перед релизом

Перед release стоит убрать из git:

- `__pycache__/` и `*.pyc` (в истории уже есть `golomb_dickman.cpython-314.pyc`);
- при желании — тяжёлые чекпоинты `gd_rigorous_ckpt_40000/` (~70 MB):
  их можно вынести в release asset или Git LFS, а не хранить в обычном git.

Добавьте `.gitignore`:

```gitignore
__pycache__/
*.py[cod]
```

Это не влияет на математику, но делает релиз аккуратнее.

## 7. Порядок действий

1. Проверить `python verify_40000.py` — должно быть `OK`.
2. Закоммитить `40000/README.md`, `verify_40000.py`, `OEIS_SUBMISSION.md`,
   `PUBLICATION_CHECKLIST.md`.
3. Сделать GitHub release `v40000`.
4. Получить Zenodo DOI и добавить его в README.
5. Отправить b-file в OEIS (или попросить редактора через SeqFan).
6. При желании — arXiv и анонс в SeqFan.
