#!/bin/bash

clear

make clean >/dev/null 2>&1 # Don't print the outputs

make all >/dev/null 2>&1 # Don't print the outputs

./CWSI_test 