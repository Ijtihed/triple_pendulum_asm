"""High-accuracy reference using scipy's 8th-order Dormand-Prince integrator
(DOP853) with very tight tolerances.  Treats the asm/python RK4 results as
estimates of this ground truth.

The point of having a higher-order reference is that, at the same accuracy
target, an 8th-order method makes thousands of times fewer steps than a
4th-order one, so its accumulated arithmetic roundoff is much smaller.
DOP853 at rtol=1e-13 atol=1e-15 is well-tested gold standard for chaotic
multibody systems.
"""
import math

try:
    import numpy as np
    from scipy.integrate import solve_ivp
except ImportError:
    print("This script needs scipy.  pip install scipy")
    raise SystemExit(1)

m1, m2, m3 = 1.0, 1.0, 1.0
L1, L2, L3 = 4.0, 3.5, 3.0
g = 9.81

beta = [m1 + m2 + m3, m2 + m3, m3]
L = [L1, L2, L3]
a = [[beta[max(j, k)] * L[j] * L[k] for k in range(3)] for j in range(3)]
gbL = [g * beta[j] * L[j] for j in range(3)]


def rhs(t, y):
    t1, t2, t3, w1, w2, w3 = y
    s = [math.sin(t1), math.sin(t2), math.sin(t3)]
    c12 = math.cos(t1 - t2); s12 = math.sin(t1 - t2)
    c13 = math.cos(t1 - t3); s13 = math.sin(t1 - t3)
    c23 = math.cos(t2 - t3); s23 = math.sin(t2 - t3)

    M11, M22, M33 = a[0][0], a[1][1], a[2][2]
    M12 = a[0][1] * c12
    M13 = a[0][2] * c13
    M23 = a[1][2] * c23

    b1 = -a[0][1] * s12 * w2 * w2 - a[0][2] * s13 * w3 * w3 - gbL[0] * s[0]
    b2 = +a[0][1] * s12 * w1 * w1 - a[1][2] * s23 * w3 * w3 - gbL[1] * s[1]
    b3 = +a[0][2] * s13 * w1 * w1 + a[1][2] * s23 * w2 * w2 - gbL[2] * s[2]

    M = np.array([[M11, M12, M13], [M12, M22, M23], [M13, M23, M33]])
    alpha = np.linalg.solve(M, [b1, b2, b3])
    return [w1, w2, w3, alpha[0], alpha[1], alpha[2]]


def energy(y):
    t1, t2, t3, w1, w2, w3 = y
    c12 = math.cos(t1 - t2); c13 = math.cos(t1 - t3); c23 = math.cos(t2 - t3)
    T = 0.5 * (
        a[0][0] * w1 * w1
        + a[1][1] * w2 * w2
        + a[2][2] * w3 * w3
        + 2 * a[0][1] * c12 * w1 * w2
        + 2 * a[0][2] * c13 * w1 * w3
        + 2 * a[1][2] * c23 * w2 * w3
    )
    V = -(gbL[0] * math.cos(t1) + gbL[1] * math.cos(t2) + gbL[2] * math.cos(t3))
    return T + V


if __name__ == "__main__":
    y0 = [math.pi / 2] * 3 + [0.0] * 3
    E0 = energy(y0)
    print("scipy DOP853 reference, rtol=1e-13, atol=1e-15\n")
    print(f"step       0  E = {E0:+.15e}  dE = 0.0")
    for steps in [100, 1000, 10000, 100000, 200000]:
        t_end = steps * 1e-4
        sol = solve_ivp(
            rhs, (0.0, t_end), y0,
            method="DOP853", rtol=1e-13, atol=1e-15, dense_output=False,
        )
        y_end = sol.y[:, -1].tolist()
        E = energy(y_end)
        print(
            f"step {steps:7d}  y = ({y_end[0]:.10f}, {y_end[1]:.10f}, "
            f"{y_end[2]:.10f}, {y_end[3]:.10f}, {y_end[4]:.10f}, {y_end[5]:.10f})"
            f"  E = {E:+.15e}  dE = {E - E0:+.3e}"
        )
