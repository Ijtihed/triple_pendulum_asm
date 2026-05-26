#!/usr/bin/env bash
set -e
nasm -felf64 triple_pendulum.asm -o triple_pendulum.o
gcc  -no-pie  triple_pendulum.o   -o triple_pendulum -lm
echo "built: ./triple_pendulum"
