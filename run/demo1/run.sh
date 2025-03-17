#!/bin/bash

./HSIM --num 10 -M mass.mtx  -K stiffness.mtx -U P-1_csc.dat P-0_csc.dat --epsilon 1e-6 --metric I --zero 6
