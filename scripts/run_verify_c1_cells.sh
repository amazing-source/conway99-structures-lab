#!/bin/bash
# Verify C1's cell-level certificates (x19k format) one at a time, each under run_capped (RAM cap 3500 MB, agreed
# with C1 on 26/09).  Arguments: base names in c1_experimental/data (without extension).
cd /c/home/researcher/Desktop/c99_structures_lab
export MSYS_NO_PATHCONV=1          # Git Bash must not rewrite the C:/... argument handed to powershell
PY="C:\Users\RESEARCHER\AppData\Local\Python\pythoncore-3.14-64\python.exe"
for b in "$@"; do
  powershell -NoProfile -ExecutionPolicy Bypass -File "C:\Users\RESEARCHER\.claude\tools\run_capped.ps1" -LimitMB 3500 \
    -TimeoutSec 900 -WorkDir "C:\Users\RESEARCHER\Desktop\c99_structures_lab" \
    -Log "C:\Users\RESEARCHER\Desktop\c99_structures_lab\logs\vc1_$b.log" -Exe "$PY" \
    -ArgList "scripts/verify_c1_cells.py C:/home/researcher/Desktop/c1_experimental/data/$b" > /dev/null
  echo "$b: $(tail -1 logs/vc1_$b.log.out) | $(grep -h 'Farkas on' logs/vc1_$b.log.out | sed 's/.*= \([0-9-]*\) .*/Farkas \1/') | $(grep -o 'max RAM=[0-9]*MB' logs/vc1_$b.log) | $(grep -o 'time=[0-9]*s' logs/vc1_$b.log)"
done
