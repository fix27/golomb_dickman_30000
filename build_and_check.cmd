@echo off
rem ---------------------------------------------------------------------------
rem Build the packaged sources and run the program's own self-check (300 digits:
rem ball radius vs. the required bound + first 99 digits against OEIS A084945).
rem
rem Needs: MSYS2 UCRT64 toolchain and FLINT (Arb inside FLINT >= 3.0).
rem Adjust MSYS below if your toolchain lives elsewhere.
rem Exit code 0 = build ok and self-check passed.
rem ---------------------------------------------------------------------------
setlocal
set MSYS=C:\msys64\ucrt64\bin
set PATH=%MSYS%;%PATH%
cd /d "%~dp0"

echo [1/3] building golomb_dickman.exe
g++ -O3 -std=c++17 -pthread src\golomb_dickman.cpp -o src\golomb_dickman.exe -lflint -lgmp -lmpfr
if errorlevel 1 goto fail

echo [2/3] building chunk_radii.exe
g++ -O3 -std=c++17 src\chunk_radii.cpp -o src\chunk_radii.exe -lflint -lgmp -lmpfr
if errorlevel 1 goto fail

echo [3/3] self-check at 300 digits
src\golomb_dickman.exe 300 --check --out verify\check_300.txt
if errorlevel 1 goto fail

echo.
echo OK: build ok, self-check passed (exit code 0).
exit /b 0

:fail
echo.
echo FAILED: see the messages above.
exit /b 1
