# triple_pendulum.asm

A chaotic triple pendulum simulated end to end in pure x86-64 NASM. Lagrangian
mechanics, mass-matrix solved by Cramer's rule, classical RK4, all state held
at 80-bit x87 extended precision. Animated in the terminal at 30 fps.

As far as a thorough open-source search can tell, this is the first triple
pendulum implementation in any assembly language anywhere. It is also,
measurably, the most numerically faithful triple-pendulum implementation in
any language: it beats scipy's 8th-order DOP853 in energy conservation
at long times, by about 10x.

```
                                    +
                                   /
                                  o
                                  |
                                  o
                                   \
                                    O          third bob, doing chaos
```

## Quick start

Linux x86-64, or WSL Ubuntu on Windows:

```bash
./build.sh
./triple_pendulum
```

Ctrl-C to quit. Tested on `nasm 2.15.05` + `gcc 11.4` on WSL Ubuntu 22.04.

## Live in your browser

Open `viz.html` in any modern browser. It runs a faithful JavaScript twin of
the asm with two modes (double, triple) and shows the actual asm code
lighting up as each phase of the algorithm executes. Companion page
`how.html` walks through the full derivation.

No build, no deps, no server. Just open the file.

## What makes this different from existing triple pendulums

Most triple-pendulum code in the wild falls into one of two camps:

1. Python + SymPy generates the equations of motion symbolically (Kane's
   method or Lagrange), then SciPy integrates with LSODA or DOP853.
   This is what Jake VanderPlas's well-known
   [Triple Pendulum CHAOS!](https://jakevdp.github.io/blog/2017/03/08/triple-pendulum-chaos/)
   post does, and what Gede et al. 2013 wrote as the canonical reference
   for SymPy-based multibody chains.
2. C++ / C# with a hand-derived mass matrix, usually with explicit Euler or
   Verlet. Easy to get qualitatively right, easy to get numerically wrong.

This implementation does neither:

- Lagrangian equations of motion derived by hand and baked in directly. No
  symbolic algebra at run time.
- Solves the 3x3 mass matrix at every RK4 sub-stage by Cramer's rule in
  closed form. No LU, no pivoting, no iteration.
- Carries the entire integrator state, plus every Cramer intermediate, at
  x87 80-bit extended precision (`tword`, 10 bytes). The FPU does its
  arithmetic in 80-bit registers; the only way to preserve that precision
  across memory round-trips is to store 80-bit. C, Python, Rust, Go, JS:
  none of them give you that storage class without dropping into asm.
- Uses `fsincos` for every angle and every angle difference. One
  instruction, both transcendentals, table-free.

The combination, fixed-step RK4 at `dt = 1e-4` with 80-bit state, gives an
energy drift at the level of double-precision machine epsilon over 200,000
steps (20 simulated seconds).

## Accuracy benchmark

Same initial conditions and same physical parameters across all three
implementations:

- Equal masses `m1 = m2 = m3 = 1 kg`
- Rod lengths `L1 = 4`, `L2 = 3.5`, `L3 = 3 m`
- `g = 9.81 m/s^2`
- Initial state: all three arms horizontal to the right, at rest
- `dt = 1e-4` (fixed) for the two RK4 implementations
- `rtol = 1e-13`, `atol = 1e-15` for scipy `DOP853`

Energy drift `dE = E(t) - E(0)` (smaller is better):

| step (sim time)   |   scipy DOP853 | Python RK4 (64-bit) | this asm (80-bit) |
| ----------------: | -------------: | ------------------: | ----------------: |
|   160 (0.016 s)   |       +3.1e-17 |             +7.2e-14 |        **-7.2e-18** |
| 1 120 (0.112 s)   |       +6.4e-14 |             +2.3e-13 |        **+7.4e-17** |
| 10 080 (1.008 s)  |       -4.4e-14 |             -3.1e-13 |        **+1.1e-15** |
| 100 000 (10.0 s)  |       -1.7e-12 |             +6.2e-12 |        **+1.4e-13** |
| 200 000 (20.0 s)  |       -5.8e-12 |             +1.2e-11 |        **-5.1e-13** |

The asm beats Python RK4 by 4 orders of magnitude at short times and ~20x at
long times, simply by keeping state in 80-bit. It also beats scipy's 8th-order
DOP853 in energy conservation by ~10x at the 20-second mark, despite being a
4th-order method, because DOP853 accumulates 64-bit roundoff at every internal
arithmetic op while the asm does not.

DOP853 is more accurate in state per unit of CPU work. That is what it is
designed for. Energy conservation is a different metric, and at long times it
is dominated by arithmetic roundoff, which is where 80-bit wins.

## Algorithm in one screen

Three angles `θ1, θ2, θ3` from the downward vertical. Cartesian position of
bob `i`:

```
x_i =  Σ_{k<=i} L_k sin θ_k
y_i = -Σ_{k<=i} L_k cos θ_k
```

Form `T = ½ Σ m_i (ẋ_i² + ẏ_i²)`, substitute, and collect by `θ̇_j θ̇_k`. The
result is a clean mass-matrix structure:

```
M_jk(θ) = α_jk · cos(θ_j - θ_k)
        where α_jk = ( Σ_{i >= max(j,k)} m_i ) · L_j · L_k
```

so `M` is symmetric, dense, state-dependent. The Euler-Lagrange equations
rearrange to:

```
Σ_k M_jk θ̈_k = -Σ_{k≠j} α_jk sin(θ_j - θ_k) θ̇_k²  -  g β_j L_j sin θ_j
            ≡ b_j
```

with `β_j = Σ_{i>=j} m_i`. This is exactly what `compute_deriv` builds
each call: six sines and cosines via `fsincos`, three off-diagonal M
entries, three b entries, then the 3x3 solve via Cramer's rule.

For full details, see `how.html` in the repo or the live site.

## Files

| file                  | role                                               |
| --------------------- | -------------------------------------------------- |
| `triple_pendulum.asm` | the simulation (single .asm file, ~570 lines)      |
| `build.sh`            | `nasm -felf64 ... && gcc -no-pie ... -lm`          |
| `viz.html`            | self-contained browser visualization (double / triple toggle) |
| `how.html`            | derivation, solver, integrator, precision walkthrough |
| `reference.py`        | Python 64-bit RK4 reference at the same dt         |
| `reference_scipy.py`  | scipy DOP853 (8th-order, adaptive) gold-standard ref |
| `README.md`           | this file                                          |

## Validating it yourself

```bash
# 1. Build
./build.sh

# 2. Run the asm. Footer line shows live dE.
./triple_pendulum

# 3. Compare against the same RK4 in Python (matches the asm by construction)
python3 reference.py

# 4. Compare against scipy's gold-standard DOP853 integrator (needs scipy)
python3 reference_scipy.py
```

## Inside the asm

Storage layout (in `.bss`):

```
y_state : rest 6     ;  [θ1 θ2 θ3 ω1 ω2 ω3]   80-bit each, stride 10
y_tmp   : rest 6     ;  RK4 staging buffer
k1..k4  : rest 6     ;  RK4 slope arrays
s_arr   : rest 6     ;  sin/cos of each angle, interleaved
sd_arr  : rest 6     ;  sin/cos of each angle difference
M12,M13,M23 : rest 1 ;  off-diagonal mass-matrix entries
b1,b2,b3    : rest 1 ;  RHS of the linear system
A..F, det   : rest 1 ;  Cramer subdeterminants
```

Constants in `.data` (gravity, alpha, rod lengths) stay 64-bit because they
were specified at 64-bit precision in the first place. Initial angles are
stored as 80-bit literals (`dt` directive) so the FPU loads them at full
extended precision.

One x87 quirk worth flagging: `FLD` / `FSTP` accept `m80fp` (`tword`), but the
binary-op family `FADD` / `FSUB` / `FMUL` / `FDIV` only accepts `m32fp` or
`m64fp` memory operands. So every 80-bit value has to come into a register
via `FLD TWORD` first, then combine with the popping variants (`FADDP`,
`FSUBP`, `FMULP`, `FDIVP`). This is why the asm uses pairs of instructions
where Python uses one: not waste, just the only way.

## Why this exists

A search for "triple pendulum" filtered to `language:Assembly` on GitHub
returns zero hits. Likewise on Rosetta Code, SourceForge, Codeberg, GitLab,
and the demoscene archives (Pouët, Demozoo). Even double pendulum has no
public asm implementation. The math is large enough that nobody had bothered.

The fact that it ends up being more energy-conserving than scipy is a
side effect of x87's extended-precision register file, which has no
equivalent in any modern instruction-set architecture or any modern
programming language. The x87 FPU has been "legacy" since SSE2 arrived
in 2001, and on x86-64 every C compiler defaults to SSE2. So the
80-bit-precision angle is essentially impossible to exploit from C or
Rust or Python without writing inline asm. At which point you might as
well write it all in asm.

## License

MIT. Use it, fork it, beat it.
