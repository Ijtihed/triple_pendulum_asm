"""Double pendulum: Python 64-bit RK4 reference. Same equations as the JS twin
in viz.html, so a direct comparison validates the visualization."""
import math

m1, m2 = 1.0, 1.0
L1, L2 = 4.0, 3.0
g = 9.81
dt = 0.0001

a11 = (m1 + m2) * L1 * L1
a22 = m2 * L2 * L2
a12 = m2 * L1 * L2
gb1L1 = g * (m1 + m2) * L1
gb2L2 = g * m2 * L2


def deriv(y):
    t1, t2, w1, w2 = y
    s1 = math.sin(t1); s2 = math.sin(t2)
    s12 = math.sin(t1 - t2); c12 = math.cos(t1 - t2)
    M11, M22 = a11, a22
    M12 = a12 * c12
    b1 = -a12 * s12 * w2 * w2 - gb1L1 * s1
    b2 = +a12 * s12 * w1 * w1 - gb2L2 * s2
    det = M11 * M22 - M12 * M12
    return (w1, w2, (b1 * M22 - M12 * b2) / det, (M11 * b2 - M12 * b1) / det)


def rk4_step(y):
    k1 = deriv(y)
    y2 = tuple(y[i] + 0.5 * dt * k1[i] for i in range(4))
    k2 = deriv(y2)
    y3 = tuple(y[i] + 0.5 * dt * k2[i] for i in range(4))
    k3 = deriv(y3)
    y4 = tuple(y[i] + dt * k3[i] for i in range(4))
    k4 = deriv(y4)
    return tuple(y[i] + dt / 6.0 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) for i in range(4))


def energy(y):
    t1, t2, w1, w2 = y
    c12 = math.cos(t1 - t2)
    T = 0.5 * (a11 * w1 * w1 + a22 * w2 * w2 + 2 * a12 * c12 * w1 * w2)
    V = -(gb1L1 * math.cos(t1) + gb2L2 * math.cos(t2))
    return T + V


if __name__ == "__main__":
    y = (math.pi / 2, math.pi / 2, 0.0, 0.0)
    E0 = energy(y)
    print(f"step       0  E = {E0:+.15e}  dE = 0.0")
    snapshots = [100, 1000, 10000, 100000, 200000]
    for step in range(1, max(snapshots) + 1):
        y = rk4_step(y)
        if step in snapshots:
            print(
                f"step {step:7d}  y = ({y[0]:.10f}, {y[1]:.10f}, "
                f"{y[2]:.10f}, {y[3]:.10f})  "
                f"E = {energy(y):+.15e}  dE = {energy(y) - E0:+.3e}"
            )
