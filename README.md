# triple_pendulum.asm

triple pendulum in pure x86-64 nasm. lagrangian mechanics, 3x3 mass matrix
solved by cramer's rule, classical rk4, state and scratch held at x87 80-bit
extended precision. terminal animation at 30 fps.

as far as a thorough open-source search can tell, no triple pendulum in any
assembly language exists publicly. neither does double. so this exists. the
80-bit precision happens to make the asm beat scipy dop853 in energy
conservation at long times by about 10x.

## build & run

linux x86-64, or wsl ubuntu on windows.

```bash
./build.sh
./triple_pendulum
```

ctrl-c to quit. tested with nasm 2.15.05 + gcc 11.4.

## browser

`index.html` is a single-file js port of the same algorithm with a double /
triple toggle. shows the actual asm code lighting up as each phase of the
algorithm runs. `how.html` is the full derivation. no build, no deps, just
open the file.

## why

filtered `language:Assembly` + triple pendulum on github. zero. tried double
pendulum, n-pendulum, multi pendulum. zero. sourceforge, gitlab, codeberg,
rosetta code, pouët, demozoo. zero. the closest hit was an 8085 stepper-motor
demo for a single pendulum.

the math is in every classical-mechanics textbook. nobody had bothered to
write it in asm.

## what's inside

- lagrangian for three point masses on rigid rods
- 3x3 mass matrix `M_jk = α_jk · cos(θ_j - θ_k)` solved by cramer's rule
- six shared 2x2 subdeterminants, three divides per derivative call
- classical rk4 at `dt = 1e-4`
- six `fsincos` calls per derivative. one instruction, both transcendentals,
  table free
- state and cramer scratch in 80-bit `tword` storage. every memory round trip
  preserves register precision
- ansi cursor escapes for the terminal renderer

## accuracy

same physics, same initial conditions, same `dt = 1e-4` for the rk4 columns.
`dE = E(t) - E(0)`. smaller is better.

| step    | scipy dop853 | python rk4 64-bit | asm rk4 80-bit |
| ------: | -----------: | ----------------: | -------------: |
|     160 |     +3.1e-17 |          +7.2e-14 |       -7.2e-18 |
|   1 120 |     +6.4e-14 |          +2.3e-13 |       +7.4e-17 |
|  10 080 |     -4.4e-14 |          -3.1e-13 |       +1.1e-15 |
| 100 000 |     -1.7e-12 |          +6.2e-12 |       +1.4e-13 |
| 200 000 |     -5.8e-12 |          +1.2e-11 |       -5.1e-13 |

asm beats python rk4 by 4 orders of magnitude at short times and ~20x at
long times, by keeping state in 80-bit. asm also beats scipy dop853 at long
times, even though dop853 is 8th-order, because dop853 accumulates 64-bit
roundoff per internal op. dop853 is more accurate in state per unit cpu
work, which is what it's designed for. energy conservation is a different
metric.

## x87 quirk

`fld` and `fstp` accept m32fp, m64fp, m80fp. but `fadd, fsub, fmul, fdiv`
only accept m32fp and m64fp. so every 80-bit operand has to come into a
register via `fld tword` first, then combine with `faddp, fsubp, fmulp,
fdivp`. about 100 extra instructions over a 64-bit-state version. no way
around it.

## validate

```bash
# build and run the asm
./build.sh
./triple_pendulum

# python rk4 reference at the same dt (matches asm by construction)
python3 reference.py
python3 reference_double.py

# scipy dop853 ground truth (needs scipy)
python3 reference_scipy.py
python3 reference_double_scipy.py
```

## files

```
triple_pendulum.asm        the simulation
build.sh                   nasm + gcc + lm
index.html                 browser visualization, double / triple toggle
how.html                   derivation, solver, integrator, precision
reference.py               python 64-bit rk4, triple
reference_scipy.py         scipy dop853 ground truth, triple
reference_double.py        python 64-bit rk4, double
reference_double_scipy.py  scipy dop853 ground truth, double
vercel.json                static deploy config
```

## license

mit.
