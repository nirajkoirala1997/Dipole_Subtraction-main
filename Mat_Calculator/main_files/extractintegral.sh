#!/bin/bash

grep -Po 'INT.*\)' mat_amp.log | sed 's/INT(//g' | sed 's/)//g' | sed 's/,/ /g' | sort | uniq >> integral.m

grep  F.* integral.m | sort | uniq > integralsto.reduze
