# Coq-формализация: по файлу на лемму

В этой папке собрана Coq-версия алгебраической части Леммы 7
(детерминантное тождество и стягивание S-дроби в J-дробь).
Каждая лемма вынесена в отдельный файл; общие определения — в `Definitions.v`.

## Структура

| Файл | Что внутри | Зависит от |
|---|---|---|
| `Definitions.v` | `aCoeff`, `bCoeff`, `An`, `Bn`, `Ddet`, `aProd`, `Ps`, `Qs`, `jb`, `ja`, `Rs`, `Ss` | — |
| `Ddet_succ.v` | `Lemma Ddet_succ` | `Definitions.v` |
| `Ddet_eq.v` | `Theorem Ddet_eq` | `Definitions.v`, `Ddet_succ.v` |
| `Ps_step3.v` | `Lemma Ps_step3` | `Definitions.v` |
| `Qs_step3.v` | `Lemma Qs_step3` | `Definitions.v` |
| `S_to_J_contraction.v` | `Theorem S_to_J_contraction` | `Definitions.v`, `Ps_step3.v`, `Qs_step3.v` |
| `All.v` | импортирует все файлы выше | все |

`Definitions.v` не содержит теорем — только определения.
`All.v` удобно использовать как единую точку входа.
Логический префикс проекта — `GolombDickmanCoq` (см. `_CoqProject`);
внутри определений используется модуль `GolombDickman`.

## Важное исправление индексации `jb`

В исходном `Формальное доказательство.coq` определение `jb` было сдвинуто на единицу:

```text
было:  jb (m+1) = alpha_{2m+1} + alpha_{2m+2}
стало: jb (m+1) = alpha_{2m+2} + alpha_{2m+3}
```

Исправленный вариант соответствует `lean/Lemma7S2J.lean`.
Без этого исправления `Ps_step3` и `S_to_J_contraction` неверны уже при `m = 2`
(коэффициент при `X` в `Ps 4` равен `alpha_2 + alpha_3`, а в `Rs 2` получался
`alpha_1 + alpha_2`).

## Исходные файлы

Три исходных файла сохранены без изменений в `original/`:

- `original/deepseek_coq_20261113_f1ac78.coq` — альтернативное доказательство
  детерминантной части;
- `original/Для S→J-стяжки.coq` — ранняя версия с `admit` и другой индексацией;
- `original/Формальное доказательство.coq` — сводный файл, из которого
  выделены новые файлы (с исправленным `jb`).

Новые файлы не перезаписывают оригиналы.

Отдельная старая версия с более общими леммами
`s_to_j_contraction_step` и `continuant_determinant_step` лежит в
`../old_version/archive/GolombDickmanCore.v`; в эту цепочку она не входит.

## Сборка

### Вариант 1: `coq_makefile` + `make`

```bash
coq_makefile -f _CoqProject -o Makefile.coq
make -f Makefile.coq
```

Или просто:

```bash
make
```

`Makefile` в этой папке делает то же самое.

### Вариант 2: Windows PowerShell

```powershell
powershell -ExecutionPolicy Bypass -File build.ps1
```

### Вариант 3: вручную

```bash
coqc -Q . GolombDickmanCoq Definitions.v
coqc -Q . GolombDickmanCoq Ddet_succ.v
coqc -Q . GolombDickmanCoq Ddet_eq.v
coqc -Q . GolombDickmanCoq Ps_step3.v
coqc -Q . GolombDickmanCoq Qs_step3.v
coqc -Q . GolombDickmanCoq S_to_J_contraction.v
coqc -Q . GolombDickmanCoq All.v
```

Требуется Coq с установленной библиотекой Mathematical Components
(`mathcomp`): используются `all_ssreflect`, `all_algebra`, `poly`.

## Проверка

Чек-лист после сборки:

- [ ] все файлы из `_CoqProject` компилируются без ошибок;
- [ ] в новых файлах нет `admit`, `Admitted`, `Axiom`, `Parameter`;
- [ ] `Print Assumptions S_to_J_contraction.` не показывает пользовательских аксиом;
- [ ] `All.v` компилируется и импортирует всю цепочку.

Пример проверки отсутствия `admit`/аксиом:

```bash
grep -RInE "admit|Admitted|Axiom|Parameter" --include="*.v" .
```

## Статус проверки

**Файлы не компилировались в среде подготовки:** Coq и mathcomp здесь не
установлены. Поэтому это структурированная версия доказательств, а не
подтверждённый сборкой результат. Перед использованием обязательно
прогоните сборку из раздела выше.

Если `Ddet_succ` не проходит из-за того, что `ring` не раскрывает `bCoeff`,
можно заменить доказательство на разбор чётности из
`original/deepseek_coq_20261113_f1ac78.coq`.
