"""Double pendulum: scipy DOP853 (8th order, adaptive) gold-standard reference."""
import math

try:
    import numpy as np
    from scipy.integrate import solve_ivp
except ImportError:
    print("Needs scipy. pip install scipy")
    raise SystemExit(1)

m1, m2 = 1.0, 1.0
L1, L2 = 4.0, 3.0
g = 9.81

a11 = (m1 + m2) * L1 * L1
a22 = m2 * L2 * L2
a12 = m2 * L1 * L2
gb1L1 = g * (m1 + m2) * L1
gb2L2 = g * m2 * L2


def rhs(t, y):
    t1, t2, w1, w2 = y
    s1 = math.sin(t1); s2 = math.sin(t2)
    s12 = math.sin(t1 - t2); c12 = math.cos(t1 - t2)
    M11, M22 = a11, a22
    M12 = a12 * c12
    b1 = -a12 * s12 * w2 * w2 - gb1L1 * s1
    b2 = +a12 * s12 * w1 * w1 - gb2L2 * s2
    det = M11 * M22 - M12 * M12
    return [w1, w2, (b1 * M22 - M12 * b2) / det, (M11 * b2 - M12 * b1) / det]


def energy(y):
    t1, t2, w1, w2 = y
    c12 = math.cos(t1 - t2)
    T = 0.5 * (a11 * w1 * w1 + a22 * w2 * w2 + 2 * a12 * c12 * w1 * w2)
    V = -(gb1L1 * math.cos(t1) + gb2L2 * math.cos(t2))
    return T + V


if __name__ == "__main__":
    y0 = [math.pi / 2, math.pi / 2, 0.0, 0.0]
    E0 = energy(y0)
    print("double pendulum, scipy DOP853, rtol=1e-13, atol=1e-15\n")
    print(f"step       0  E = {E0:+.15e}  dE = 0.0")
    for steps in [100, 1000, 10000, 100000, 200000]:
        t_end = steps * 1e-4
        sol = solve_ivp(rhs, (0.0, t_end), y0, method="DOP853", rtol=1e-13, atol=1e-15)
        y_end = sol.y[:, -1].tolist()
        E = energy(y_end)
        print(
            f"step {steps:7d}  y = ({y_end[0]:.10f}, {y_end[1]:.10f}, "
            f"{y_end[2]:.10f}, {y_end[3]:.10f})  E = {E:+.15e}  dE = {E - E0:+.3e}"
        )
