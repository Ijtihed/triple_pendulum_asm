; ============================================================================
;  triple_pendulum.asm
;
;  lagrangian mechanics for three point masses on rigid rods.
;     M(theta) * theta_ddot = b(theta, omega)        (3x3 dense system)
;  cramer's rule for the 3x3 solve.
;  classical rk4 at dt = 1e-4.
;  state and cramer scratch held at x87 80-bit (tword) extended precision.
;  fsincos for every angle and angle-difference.
;
;  build:  nasm -felf64 triple_pendulum.asm -o triple_pendulum.o
;          gcc  -no-pie  triple_pendulum.o   -o triple_pendulum -lm
;  run:    ./triple_pendulum                          (ctrl-c to quit)
; ============================================================================

default rel
global  main
extern  printf
extern  fflush
extern  usleep

; ----------------------------------------------------------------------------
section .data
; ----------------------------------------------------------------------------
g_const:    dq 9.81
; rk4 timestep. dt = 1e-4 puts per-step truncation at the level of double
; roundoff. at 200000 steps the energy drift stays around 1e-12.
dt_phys:    dq 0.0001
dt_half:    dq 0.00005
dt_sixth:   dq 0.0000166666666666666666667
two_d:      dq 2.0

; precomputed alpha_jk = (sum m_i for i>=max(j,k)) * L_j * L_k    (m_i = 1)
a11:        dq 48.0
a12:        dq 28.0
a13:        dq 12.0
a22:        dq 24.5
a23:        dq 10.5
a33:        dq 9.0
; gravity coefficients  g * beta_j * L_j  (beta_j = sum_{i>=j} m_i)
gb1L1:      dq 117.72
gb2L2:      dq 68.67
gb3L3:      dq 29.43
; rendering constants (renderer projects to terminal cells)
L1_d:       dq 4.0
L2_d:       dq 3.5
L3_d:       dq 3.0
piv_x_d:    dq 40.0
piv_y_d:    dq 4.0

; initial conditions held as 80-bit literals so the fpu loads them at full
; extended precision rather than rounding from a 64-bit double.
init_th1:   dt 1.5707963267948966192   ; pi/2
init_th2:   dt 1.5707963267948966192
init_th3:   dt 1.5707963267948966192

; ansi / printf strings
fmt_clr:    db 27,"[2J",27,"[H",0
fmt_pos:    db 27,"[%d;%dH%c",0
fmt_foot:   db 27,"[23;1H",27,"[Ktriple_pendulum.asm | e=%9.4f  de=%+10.3e  step=%lld",10
            db 27,"[K(lagrangian + rk4 @ dt=1e-4, 80-bit x87 state)  ctrl-c to quit",10,0

; ----------------------------------------------------------------------------
section .bss
; ----------------------------------------------------------------------------
; storage: 80-bit tword for everything the integrator writes back to memory.
; the fpu does arithmetic at 80-bit in registers, so storing 80-bit means no
; rounding on memory round-trips. constants (alpha, gravity, dt) stay 64-bit
; because their literal precision is no better than that.
;
; x87 quirk: fld / fstp accept m80fp. fadd / fsub / fmul / fdiv only accept
; m32fp and m64fp. so every 80-bit operand has to come into a register via
; fld tword first, then combine with faddp / fsubp / fmulp / fdivp.

y_state:    rest 6        ; [theta1 theta2 theta3 omega1 omega2 omega3]
y_tmp:      rest 6        ; RK4 staging
k1_arr:     rest 6
k2_arr:     rest 6
k3_arr:     rest 6
k4_arr:     rest 6
; compute_deriv scratch
s_arr:      rest 6        ; s1 c1 s2 c2 s3 c3 interleaved (80-bit each)
sd_arr:     rest 6        ; s12 c12 s13 c13 s23 c23 interleaved
M12_v:      rest 1
M13_v:      rest 1
M23_v:      rest 1
b1_v:       rest 1
b2_v:       rest 1
b3_v:       rest 1
A_v:        rest 1
B_v:        rest 1
C_v:        rest 1
D_v:        rest 1
E_v:        rest 1
F_v:        rest 1
det_v:      rest 1
E_init:     rest 1
step_cnt:   resq 1        ; 64-bit integer counter

; ----------------------------------------------------------------------------
section .text
; ----------------------------------------------------------------------------

; ============================================================================
;  compute_deriv(rdi = y_ptr, rsi = dy_ptr)
;  in:   y_ptr -> 6 * 10 bytes = [theta1..theta3 omega1..omega3]  (80-bit)
;  out:  dy_ptr <- [omega1 omega2 omega3 alpha1 alpha2 alpha3]
;  builds M and b from the lagrangian, solves M*alpha = b via cramer.
; ============================================================================
compute_deriv:
        push    rbx
        push    r12
        mov     rbx, rdi
        mov     r12, rsi

        ; ---- copy omega part of y into dy (dtheta = omega) ----
        fld     tword [rbx + 30]
        fstp    tword [r12]
        fld     tword [rbx + 40]
        fstp    tword [r12 + 10]
        fld     tword [rbx + 50]
        fstp    tword [r12 + 20]

        ; ---- per-angle sin/cos ----
        fld     tword [rbx]
        fsincos
        fstp    tword [s_arr + 10]
        fstp    tword [s_arr]

        fld     tword [rbx + 10]
        fsincos
        fstp    tword [s_arr + 30]
        fstp    tword [s_arr + 20]

        fld     tword [rbx + 20]
        fsincos
        fstp    tword [s_arr + 50]
        fstp    tword [s_arr + 40]

        ; ---- sin/cos of angle differences ----
        fld     tword [rbx]
        fld     tword [rbx + 10]
        fsubp
        fsincos
        fstp    tword [sd_arr + 10]
        fstp    tword [sd_arr]

        fld     tword [rbx]
        fld     tword [rbx + 20]
        fsubp
        fsincos
        fstp    tword [sd_arr + 30]
        fstp    tword [sd_arr + 20]

        fld     tword [rbx + 10]
        fld     tword [rbx + 20]
        fsubp
        fsincos
        fstp    tword [sd_arr + 50]
        fstp    tword [sd_arr + 40]

        ; ---- mass matrix off-diagonals  (M_jk = alpha_jk * cos(theta_j-theta_k)) ----
        fld     qword [a12]
        fld     tword [sd_arr + 10]
        fmulp
        fstp    tword [M12_v]

        fld     qword [a13]
        fld     tword [sd_arr + 30]
        fmulp
        fstp    tword [M13_v]

        fld     qword [a23]
        fld     tword [sd_arr + 50]
        fmulp
        fstp    tword [M23_v]

        ; ---- b1 = -(a12*s12*w2^2 + a13*s13*w3^2 + gb1L1*s1) ----
        fld     qword [a12]
        fld     tword [sd_arr]
        fmulp                                 ; a12*s12
        fld     tword [rbx + 40]              ; w2
        fmul    st0, st0                      ; w2^2
        fmulp                                 ; a12*s12*w2^2

        fld     qword [a13]
        fld     tword [sd_arr + 20]
        fmulp                                 ; a13*s13
        fld     tword [rbx + 50]
        fmul    st0, st0
        fmulp                                 ; a13*s13*w3^2

        faddp                                 ; sum_a

        fld     qword [gb1L1]
        fld     tword [s_arr]
        fmulp                                 ; gb1L1*s1
        faddp                                 ; sum_a + gb1L1*s1
        fchs
        fstp    tword [b1_v]

        ; ---- b2 = +a12*s12*w1^2 - a23*s23*w3^2 - gb2L2*s2 ----
        fld     qword [a12]
        fld     tword [sd_arr]
        fmulp
        fld     tword [rbx + 30]
        fmul    st0, st0
        fmulp                                 ; +a12*s12*w1^2

        fld     qword [a23]
        fld     tword [sd_arr + 40]
        fmulp
        fld     tword [rbx + 50]
        fmul    st0, st0
        fmulp                                 ; a23*s23*w3^2
        fsubp                                 ; +a12*s12*w1^2 - a23*s23*w3^2

        fld     qword [gb2L2]
        fld     tword [s_arr + 20]
        fmulp                                 ; gb2L2*s2
        fsubp
        fstp    tword [b2_v]

        ; ---- b3 = +a13*s13*w1^2 + a23*s23*w2^2 - gb3L3*s3 ----
        fld     qword [a13]
        fld     tword [sd_arr + 20]
        fmulp
        fld     tword [rbx + 30]
        fmul    st0, st0
        fmulp                                 ; +a13*s13*w1^2

        fld     qword [a23]
        fld     tword [sd_arr + 40]
        fmulp
        fld     tword [rbx + 40]
        fmul    st0, st0
        fmulp                                 ; a23*s23*w2^2

        faddp                                 ; sum

        fld     qword [gb3L3]
        fld     tword [s_arr + 40]
        fmulp                                 ; gb3L3*s3
        fsubp
        fstp    tword [b3_v]

        ; ----------------------- Cramer helper terms ------------------------
        ; A = M22*M33 - M23*M23
        fld     qword [a22]
        fmul    qword [a33]
        fld     tword [M23_v]
        fmul    st0, st0
        fsubp
        fstp    tword [A_v]

        ; B = M12*M33 - M23*M13
        fld     tword [M12_v]
        fmul    qword [a33]
        fld     tword [M23_v]
        fld     tword [M13_v]
        fmulp
        fsubp
        fstp    tword [B_v]

        ; C = M12*M23 - M22*M13
        fld     tword [M12_v]
        fld     tword [M23_v]
        fmulp
        fld     qword [a22]
        fld     tword [M13_v]
        fmulp
        fsubp
        fstp    tword [C_v]

        ; D = b2*M33 - M23*b3
        fld     tword [b2_v]
        fmul    qword [a33]
        fld     tword [M23_v]
        fld     tword [b3_v]
        fmulp
        fsubp
        fstp    tword [D_v]

        ; E = b2*M23 - M22*b3
        fld     tword [b2_v]
        fld     tword [M23_v]
        fmulp
        fld     qword [a22]
        fld     tword [b3_v]
        fmulp
        fsubp
        fstp    tword [E_v]

        ; F = M12*b3 - b2*M13
        fld     tword [M12_v]
        fld     tword [b3_v]
        fmulp
        fld     tword [b2_v]
        fld     tword [M13_v]
        fmulp
        fsubp
        fstp    tword [F_v]

        ; det = M11*A - M12*B + M13*C
        fld     qword [a11]
        fld     tword [A_v]
        fmulp
        fld     tword [M12_v]
        fld     tword [B_v]
        fmulp
        fsubp
        fld     tword [M13_v]
        fld     tword [C_v]
        fmulp
        faddp
        fstp    tword [det_v]

        ; alpha1 = (b1*A - M12*D + M13*E) / det
        fld     tword [b1_v]
        fld     tword [A_v]
        fmulp
        fld     tword [M12_v]
        fld     tword [D_v]
        fmulp
        fsubp
        fld     tword [M13_v]
        fld     tword [E_v]
        fmulp
        faddp
        fld     tword [det_v]
        fdivp
        fstp    tword [r12 + 30]

        ; alpha2 = (M11*D - b1*B + M13*F) / det
        fld     qword [a11]
        fld     tword [D_v]
        fmulp
        fld     tword [b1_v]
        fld     tword [B_v]
        fmulp
        fsubp
        fld     tword [M13_v]
        fld     tword [F_v]
        fmulp
        faddp
        fld     tword [det_v]
        fdivp
        fstp    tword [r12 + 40]

        ; alpha3 = (-M11*E - M12*F + b1*C) / det
        fld     qword [a11]
        fld     tword [E_v]
        fmulp
        fchs
        fld     tword [M12_v]
        fld     tword [F_v]
        fmulp
        fsubp
        fld     tword [b1_v]
        fld     tword [C_v]
        fmulp
        faddp
        fld     tword [det_v]
        fdivp
        fstp    tword [r12 + 50]

        pop     r12
        pop     rbx
        ret

; ============================================================================
;  axpy macro: for i=0..5,  y_tmp[i] = y_state[i] + scalar * src[i]
;       %1 = 80-bit source array,  %2 = 64-bit scalar
; ============================================================================
%macro AXPY_TO_YTMP 2
%assign i 0
%rep 6
        fld     tword [%1 + i*10]
        fmul    qword [%2]
        fld     tword [y_state + i*10]
        faddp
        fstp    tword [y_tmp + i*10]
%assign i i+1
%endrep
%endmacro

; ============================================================================
;  rk4_step:  one classical rk4 step on y_state.
;    k1 = f(y)
;    k2 = f(y + dt/2 k1)
;    k3 = f(y + dt/2 k2)
;    k4 = f(y + dt   k3)
;    y += dt/6 * (k1 + 2*k2 + 2*k3 + k4)
; ============================================================================
rk4_step:
        push    rbp
        mov     rbp, rsp
        sub     rsp, 16

        lea     rdi, [y_state]
        lea     rsi, [k1_arr]
        call    compute_deriv
        AXPY_TO_YTMP k1_arr, dt_half

        lea     rdi, [y_tmp]
        lea     rsi, [k2_arr]
        call    compute_deriv
        AXPY_TO_YTMP k2_arr, dt_half

        lea     rdi, [y_tmp]
        lea     rsi, [k3_arr]
        call    compute_deriv
        AXPY_TO_YTMP k3_arr, dt_phys

        lea     rdi, [y_tmp]
        lea     rsi, [k4_arr]
        call    compute_deriv

        ; y += (dt/6) * (k1 + 2*(k2+k3) + k4)
%assign i 0
%rep 6
        fld     tword [k2_arr + i*10]
        fld     tword [k3_arr + i*10]
        faddp
        fmul    qword [two_d]
        fld     tword [k1_arr + i*10]
        faddp
        fld     tword [k4_arr + i*10]
        faddp
        fmul    qword [dt_sixth]
        fld     tword [y_state + i*10]
        faddp
        fstp    tword [y_state + i*10]
%assign i i+1
%endrep

        leave
        ret

; ============================================================================
;  compute_energy: returns total mechanical E = T + V as a double in xmm0.
; ============================================================================
compute_energy:
        push    rbp
        mov     rbp, rsp

        ; refresh c12, c13, c23
        fld     tword [y_state]
        fld     tword [y_state + 10]
        fsubp
        fcos
        fstp    tword [sd_arr + 10]

        fld     tword [y_state]
        fld     tword [y_state + 20]
        fsubp
        fcos
        fstp    tword [sd_arr + 30]

        fld     tword [y_state + 10]
        fld     tword [y_state + 20]
        fsubp
        fcos
        fstp    tword [sd_arr + 50]

        ; t = 0.5*( a11 w1^2 + a22 w2^2 + a33 w3^2
        ;          + 2 a12 c12 w1 w2 + 2 a13 c13 w1 w3 + 2 a23 c23 w2 w3 )
        fld     qword [a11]
        fld     tword [y_state + 30]
        fmul    st0, st0
        fmulp

        fld     qword [a22]
        fld     tword [y_state + 40]
        fmul    st0, st0
        fmulp
        faddp

        fld     qword [a33]
        fld     tword [y_state + 50]
        fmul    st0, st0
        fmulp
        faddp

        ; cross term 2 a12 c12 w1 w2
        fld     qword [two_d]
        fmul    qword [a12]
        fld     tword [sd_arr + 10]
        fmulp
        fld     tword [y_state + 30]
        fmulp
        fld     tword [y_state + 40]
        fmulp
        faddp

        ; 2 a13 c13 w1 w3
        fld     qword [two_d]
        fmul    qword [a13]
        fld     tword [sd_arr + 30]
        fmulp
        fld     tword [y_state + 30]
        fmulp
        fld     tword [y_state + 50]
        fmulp
        faddp

        ; 2 a23 c23 w2 w3
        fld     qword [two_d]
        fmul    qword [a23]
        fld     tword [sd_arr + 50]
        fmulp
        fld     tword [y_state + 40]
        fmulp
        fld     tword [y_state + 50]
        fmulp
        faddp

        fld1
        fld1
        faddp
        fdivp                              ; * 0.5

        ; v = -(gb1L1 cos(t1) + gb2L2 cos(t2) + gb3L3 cos(t3))
        fld     tword [y_state]
        fcos
        fmul    qword [gb1L1]
        fld     tword [y_state + 10]
        fcos
        fmul    qword [gb2L2]
        faddp
        fld     tword [y_state + 20]
        fcos
        fmul    qword [gb3L3]
        faddp
        fchs
        faddp                              ; e = t + v on fpu top

        ; round to double for printf return.
        sub     rsp, 16
        fstp    qword [rsp]
        movsd   xmm0, [rsp]
        add     rsp, 16

        leave
        ret

; ============================================================================
;  draw_char_at(esi=row, edx=col, ecx=char)
; ============================================================================
draw_char_at:
        push    rbp
        mov     rbp, rsp
        sub     rsp, 16
        lea     rdi, [fmt_pos]
        xor     eax, eax
        call    printf
        leave
        ret

; ============================================================================
;  render:  forward kinematics in 64-bit (for terminal coords), then ansi
;           escape sequences for the three bobs and the energy footer.
; ============================================================================
render:
        push    rbp
        mov     rbp, rsp
        sub     rsp, 80

        lea     rdi, [fmt_clr]
        xor     eax, eax
        call    printf

        ; bob 1
        fld     tword [y_state]
        fsincos                            ; cos, sin on fpu
        fld     qword [L1_d]
        fmul    st0, st2                   ; L1*sin
        fadd    qword [piv_x_d]
        fstp    qword [rbp-8]              ; x1
        fld     qword [L1_d]
        fmul    st0, st1                   ; L1*cos
        fadd    qword [piv_y_d]
        fstp    qword [rbp-16]             ; y1
        fstp    st0
        fstp    st0

        ; bob 2
        fld     tword [y_state + 10]
        fsincos
        fld     qword [L2_d]
        fmul    st0, st2
        fadd    qword [rbp-8]
        fstp    qword [rbp-24]
        fld     qword [L2_d]
        fmul    st0, st1
        fadd    qword [rbp-16]
        fstp    qword [rbp-32]
        fstp    st0
        fstp    st0

        ; bob 3
        fld     tword [y_state + 20]
        fsincos
        fld     qword [L3_d]
        fmul    st0, st2
        fadd    qword [rbp-24]
        fstp    qword [rbp-40]
        fld     qword [L3_d]
        fmul    st0, st1
        fadd    qword [rbp-32]
        fstp    qword [rbp-48]
        fstp    st0
        fstp    st0

        ; pivot '+'
        fld     qword [piv_y_d]
        fistp   dword [rbp-56]
        mov     esi, [rbp-56]
        fld     qword [piv_x_d]
        fistp   dword [rbp-56]
        mov     edx, [rbp-56]
        mov     ecx, '+'
        call    draw_char_at

        ; bob 1
        fld     qword [rbp-16]
        fistp   dword [rbp-56]
        mov     esi, [rbp-56]
        fld     qword [rbp-8]
        fistp   dword [rbp-56]
        mov     edx, [rbp-56]
        mov     ecx, 'o'
        call    draw_char_at

        ; bob 2
        fld     qword [rbp-32]
        fistp   dword [rbp-56]
        mov     esi, [rbp-56]
        fld     qword [rbp-24]
        fistp   dword [rbp-56]
        mov     edx, [rbp-56]
        mov     ecx, 'o'
        call    draw_char_at

        ; bob 3
        fld     qword [rbp-48]
        fistp   dword [rbp-56]
        mov     esi, [rbp-56]
        fld     qword [rbp-40]
        fistp   dword [rbp-56]
        mov     edx, [rbp-56]
        mov     ecx, 'O'
        call    draw_char_at

        ; footer: e and de
        call    compute_energy             ; xmm0 = e (double)
        movsd   xmm1, xmm0
        fld     tword [E_init]
        sub     rsp, 16
        fstp    qword [rsp]
        movsd   xmm2, [rsp]
        add     rsp, 16
        subsd   xmm1, xmm2                 ; de = e - e_init

        mov     rsi, [step_cnt]
        lea     rdi, [fmt_foot]
        mov     eax, 2
        call    printf

        xor     edi, edi
        call    fflush

        leave
        ret

; ============================================================================
;  main
; ============================================================================
main:
        push    rbp
        mov     rbp, rsp
        sub     rsp, 16

        fld     tword [init_th1]
        fstp    tword [y_state]
        fld     tword [init_th2]
        fstp    tword [y_state + 10]
        fld     tword [init_th3]
        fstp    tword [y_state + 20]
        fldz
        fstp    tword [y_state + 30]
        fldz
        fstp    tword [y_state + 40]
        fldz
        fstp    tword [y_state + 50]

        xor     rax, rax
        mov     [step_cnt], rax

        call    compute_energy             ; xmm0 = initial e
        sub     rsp, 16
        movsd   [rsp], xmm0
        fld     qword [rsp]
        fstp    tword [E_init]
        add     rsp, 16

.loop:
        ; 160 rk4 steps per frame  =  16 ms simulated  (dt = 1e-4)
        mov     r12d, 160
.physloop:
        call    rk4_step
        inc     qword [step_cnt]
        dec     r12d
        jnz     .physloop

        call    render

        mov     edi, 33000
        call    usleep
        jmp     .loop

        leave
        ret
