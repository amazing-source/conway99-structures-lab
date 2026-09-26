#!/bin/bash
# Orbit-level certificates for the triangle-free pure twisted classes (user's approval of 26/09: "Lancer les 13 classes").
# One class at a time, each under run_capped (RAM cap 4000 MB), CG time limit 270 s; then independent verification
# of every certificate and coverage of the 30 360 labelled triangle-free decompositions.
cd /c/home/researcher/Desktop/c99_structures_lab
PY="C:\Users\RESEARCHER\AppData\Local\Python\pythoncore-3.14-64\python.exe"
for n in 13 1 2 3 4 5 6 7 8 9 10 11 12; do
  nn=$(printf %02d $n)
  powershell -NoProfile -ExecutionPolicy Bypass -File "C:\Users\RESEARCHER\.claude\tools\run_capped.ps1" -LimitMB 4000 \
    -TimeoutSec 360 -WorkDir "C:\Users\RESEARCHER\Desktop\c99_structures_lab" \
    -Log "C:\Users\RESEARCHER\Desktop\c99_structures_lab\logs\orbit_z$nn.log" -Exe "$PY" -ArgList "scripts/orbit_p2_colgen.py $n 270" > /dev/null
  echo "zero$nn: $(grep -h 'CERTIFIED\|CONVERGED\|stopped\|master status' logs/orbit_z$nn.log.out | tail -1) | $(grep -h 'exact check' logs/orbit_z$nn.log.out | sed 's/.*->  //')"
done
python scripts/verify_orbit_certs.py coverage > logs/verify_orbit_certs.log 2>&1
tail -3 logs/verify_orbit_certs.log
