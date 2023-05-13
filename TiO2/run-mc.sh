#!/bin/bash -l
export ASE_ESPRESSO_COMMAND="pw.x -in PREFIX.pwi > PREFIX.pwo" 
python TiO2-mc.py
