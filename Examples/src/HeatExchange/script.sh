#!/bin/bash
echo "This script will solve the problem first evaluating the difference between two steps with L2 norm, " 
make distclean
make
./main -p parameters1.pot
./main -p parameters2.pot
