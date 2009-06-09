#!/bin/bash
#ifort -g -p -traceback OpenMP-HTF.serial.f90 -o OpenMP-HTF.serial
ifort OpenMP-HTF.serial.f90 -o OpenMP-HTF.serial
