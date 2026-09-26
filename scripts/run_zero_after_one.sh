#!/bin/sh
# Waits for the one-triangle certification to end (END line in its run log), then launches the triangle-free one,
# so that only one heavy job runs at a time (user's rule).
cd /c/home/researcher/Desktop/c99_structures_lab
until grep -q "END" logs/certify_one.runlog 2>/dev/null; do sleep 5; done
powershell -NoProfile -ExecutionPolicy Bypass -File "C:\Users\RESEARCHER\.claude\tools\run_capped.ps1" -LimitMB 4000 -CommitLimitMB 10000 -TimeoutSec 1900 -StatusSec 60 -WorkDir "C:\Users\RESEARCHER\Desktop\c99_structures_lab" -Log "C:\Users\RESEARCHER\Desktop\c99_structures_lab\logs\certify_zero.runlog" -Exe "C:\Users\RESEARCHER\AppData\Local\Python\pythoncore-3.14-64\python.exe" -ArgList "-u scripts/certify_core.py zero 1800"
