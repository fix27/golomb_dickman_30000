# A084945: отправка b-file на 40 000 термов

## Что отправляем

| Параметр | Значение |
|---|---|
| Последовательность | OEIS [A084945](https://oeis.org/A084945) — Decimal expansion of Golomb-Dickman constant |
| Файл | [`b084945.txt`](b084945.txt) |
| Термы | `a(0)..a(39999)` — 40 000 цифр после запятой |
| Формат | `n a(n)`, только данные, ASCII, LF, завершающая пустая строка |
| Размер | 308 891 байт (лимит OEIS — около 10 MB) |
| SHA-256 | `c73af643198da299d2e624e3c2f69d05439efcddbbefa7cec176ac63f993bed3` |

Текущий b-file в A084945 содержит только `n = 0..101` и синтезирован из записи,
поэтому прежних авторов указывать не нужно: новый файл просто заменяет старый.

## Текст ссылки

Замените `YOUR NAME` на ваше имя или OEIS-логин:

```text
YOUR NAME, <a href="b084945.txt">Table of n, a(n) for n=0..39999</a>
```

## Если есть доступ к редактированию

1. Войдите на **https://oeis.org/login** (или зарегистрируйтесь через
   **https://oeis.org/account**).
2. Откройте **https://oeis.org/A084945**.
3. Нажмите **`edit`** рядом с данными последовательности.
4. В секции **Links** найдите загрузку b-file и выберите `b084945.txt`.
5. Вставьте текст ссылки выше; при необходимости добавьте комментарий:
   `Terms 0..39999 computed using certified interval arithmetic (Arb).`
6. Нажмите **Save changes**, затем **These changes are ready for review by an OEIS editor**.
7. Проверьте, что в черновике ссылка ведёт на `b084945.txt` и что видны
   `n=0..39999`.

Правки не публикуются сразу: они попадают в очередь редакторов. Это нормально.

## Если прав редактировать нет

Не нужно создавать дубликат последовательности. Попросите редактора загрузить
файл за вас:

- **SeqFan**: https://groups.google.com/g/seqfan — там читают редакторы OEIS.
- Либо найдите редактора через https://oeis.org/wiki/Editorial_Board и напишите
  через его wiki-страницу (`Email this user`).

Готовое письмо:

```text
Subject: A084945 b-file extension to 40000 terms

Dear OEIS editors,

Please upload the attached b084945.txt as a b-file for A084945
(Decimal expansion of Golomb-Dickman constant), replacing the
current synthesized b-file for n=0..101.

Link text:
YOUR NAME, <a href="b084945.txt">Table of n, a(n) for n=0..39999</a>

The 40,000 digits are backed by a certified interval-arithmetic
computation: floor(10^40000 * Lambda_-) = floor(10^40000 * Lambda_+).
Certificate, verifier, and code:
https://github.com/fix27/golomb-dickman-constant/tree/main/40000

SHA-256 of b084945.txt:
c73af643198da299d2e624e3c2f69d05439efcddbbefa7cec176ac63f993bed3

The file is plain ASCII with LF line endings, 40,000 terms,
indices 0..39999, and a trailing blank line.

Thank you,
YOUR NAME
```

Если редактор попросит сократить объём, из этого же файла можно сделать
префиксы на 10 000 или 20 000 термов:

```python
from pathlib import Path

lines = Path("b084945.txt").read_text(encoding="ascii").splitlines()
for n in (10_000, 20_000):
    Path(f"b084945_{n}.txt").write_text(
        "\n".join(lines[:n]) + "\n\n", encoding="ascii"
    )
```

Здесь первые 10 000 строк — это `n=0..9999`, первые 20 000 — `n=0..19999`.

## Почему это корректно

- 40 000 термов совпадают с `golomb_dickman_40000_digits.txt`;
- цифры совпадают с сертификатом `golomb_dickman_40000_certified.txt`;
- первые 102 цифры совпадают с официальными данными A084945;
- проверка воспроизводится командой `python verify_40000.py` без сторонних пакетов.

## Ограничения

- OEIS не является реестром «мировых рекордов»: факт фиксации — это сам b-file
  и ссылка с вашим именем.
- 40 000 термов больше обычных 10 000/20 000; это допустимо, но финальное решение
  за редакторами. Если попросят уменьшить — используйте префиксы выше.
- Не загружайте сертификат или файл цифр как b-file: b-file должен содержать
  только пары `n a(n)`.
