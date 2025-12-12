#!/bin/bash
srun -N 1 --ntasks-per-node=1 --cpus-per-task=1 --mem=2G --time=99:99:99 --job-name=gtdb_unpack --pty bash