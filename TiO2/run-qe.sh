#!/bin/bash -l
export ASE_ESPRESSO_COMMAND="pw.x -in PREFIX.pwi > PREFIX.pwo" 
python3 TiO2-bs.py
