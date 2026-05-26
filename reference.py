"""triple pendulum, python 64-bit rk4 reference at dt = 1e-4."""
import math

m1, m2, m3 = 1.0, 1.0, 1.0
L1, L2, L3 = 4.0, 3.5, 3.0
g = 9.81
dt = 0.0001

beta1 = m1 + m2 + m3
beta2 = m2 + m3
beta3 = m3
a11 = beta1 * L1 * L1
a12 = beta2 * L1 * L2
a13 = beta3 * L1 * L3
a22 = beta2 * L2 * L2
a23 = beta3 * L2 * L3
a33 = beta3 * L3 * L3
gb1L1 = g * beta1 * L1
gb2L2 = g * beta2 * L2
gb3L3 = g * beta3 * L3


def deriv(y):
    t1, t2, t3, w1, w2, w3 = y
    s1, c1 = math.sin(t1), math.cos(t1)
    s2, c2 = math.sin(t2), math.cos(t2)
    s3, c3 = math.sin(t3), math.cos(t3)
    s12, c12 = math.sin(t1 - t2), math.cos(t1 - t2)
    s13, c13 = math.sin(t1 - t3), math.cos(t1 - t3)
    s23, c23 = math.sin(t2 - t3), math.cos(t2 - t3)

    M11, M22, M33 = a11, a22, a33
    M12 = a12 * c12
    M13 = a13 * c13
    M23 = a23 * c23

    b1 = -a12 * s12 * w2 * w2 - a13 * s13 * w3 * w3 - gb1L1 * s1
    b2 = +a12 * s12 * w1 * w1 - a23 * s23 * w3 * w3 - gb2L2 * s2
    b3 = +a13 * s13 * w1 * w1 + a23 * s23 * w2 * w2 - gb3L3 * s3

    A = M22 * M33 - M23 * M23
    B = M12 * M33 - M23 * M13
    C = M12 * M23 - M22 * M13
    D = b2 * M33 - M23 * b3
    E = b2 * M23 - M22 * b3
    F = M12 * b3 - b2 * M13

    det = M11 * A - M12 * B + M13 * C
    d1 = b1 * A - M12 * D + M13 * E
    d2 = M11 * D - b1 * B + M13 * F
    d3 = -M11 * E - M12 * F + b1 * C

    return (w1, w2, w3, d1 / det, d2 / det, d3 / det)


def rk4_step(y):
    k1 = deriv(y)
    y2 = tuple(y[i] + 0.5 * dt * k1[i] for i in range(6))
    k2 = deriv(y2)
    y3 = tuple(y[i] + 0.5 * dt * k2[i] for i in range(6))
    k3 = deriv(y3)
    y4 = tuple(y[i] + dt * k3[i] for i in range(6))
    k4 = deriv(y4)
    return tuple(
        y[i] + dt / 6.0 * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) for i in range(6)
    )


def energy(y):
    t1, t2, t3, w1, w2, w3 = y
    c12, c13, c23 = math.cos(t1 - t2), math.cos(t1 - t3), math.cos(t2 - t3)
    T = 0.5 * (
        a11 * w1 * w1
        + a22 * w2 * w2
        + a33 * w3 * w3
        + 2 * a12 * c12 * w1 * w2
        + 2 * a13 * c13 * w1 * w3
        + 2 * a23 * c23 * w2 * w3
    )
    V = -g * (beta1 * L1 * math.cos(t1) + beta2 * L2 * math.cos(t2) + beta3 * L3 * math.cos(t3))
    return T + V


if __name__ == "__main__":
    # initial: all three arms horizontal to the right, at rest
    y = (math.pi / 2, math.pi / 2, math.pi / 2, 0.0, 0.0, 0.0)
    E0 = energy(y)
    print(f"step       0  y = {y}  E = {energy(y):.6f}  dE = 0.0")
    snapshots = [100, 1000, 10000, 100000, 200000]
    for step in range(1, max(snapshots) + 1):
        y = rk4_step(y)
        if step in snapshots:
            print(
                f"step {step:7d}  y = ({y[0]:.10f}, {y[1]:.10f}, {y[2]:.10f}, "
                f"{y[3]:.10f}, {y[4]:.10f}, {y[5]:.10f})  E = {energy(y):.10f}  "
                f"dE = {energy(y) - E0:.3e}"
            )
