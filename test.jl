include("MACHINE_SIMULATOR.jl")
using JuMP
using OSQP
using LinearAlgebra
using Plots
using Printf
using Ipopt
# using PyPlot
# pygui(true)


# Initialize the AC machine
# I.M. has KE = 0
# P.M. has Rreq = 0
npp,   IN,    R,    Ld,    Lq,  KE, Rreq,   Js = (
  2, 10.0, 10.0, 0.015, 0.025, 0.5,  0.0, 0.05)
global ACM = The_AC_Machine(npp, IN, R, Ld, Lq, KE, Rreq, Js)


# # # Build Controller
# # 1 - Nested Loops Controller
# CLBW_iQ = 2 * π * 100s
# K = npp * KE / Js
# delta = 20
# # Q-axis 
# KP_i = Lq * CLBW_iQ
# KI_i = R / Lq
# KP_v = CLBW_iQ / (K * delta)
# KI_v = CLBW_iQ / (delta^2)
# # D-axis
# CLBW_iD = 2 * π * 100
# KP_id = Ld * CLBW_iD
# KI_id = R / Ld
# Nested_CRTL = NestedLoopsController(KP_i, KI_i, KP_v, KI_v, KP_id, KI_id)

# # 2 - MPC Controller
Ω_star, T_L = (0.0, 0.0)
dt = 0.001 # Sampling time
N = 400 # Prediction horizon
Q1, Q2, Q3 = (50.0, 0.0, 100.0) # Q1 for iD; Q2 for iQ; Q3 for Ω
R1, R2 = (1.0, 1.0) # R1 for uD; R2 for uQ
F1, F2 = (400.0, 50.0) # F1 for final Ω； F2 for iD
uD_max, uQ_max = (100.0, 100.0)
MPC_CRTL = MPControl(dt, Ld, Lq, R, npp, KE, KE, Js, Ω_star, T_L,
    N, Q1, Q2, Q3, R1, R2, F1, F2,
    uD_max, uQ_max)


println("ACM")

MACHINE_TS = 1e-4
t0 = 0 # Initial time
TIME = 4 # Simulation time duration
machine_times, watch_data = ACMSimPyIncremental(MACHINE_TS, t0, TIME, ACM, MPC_CRTL)
myplot(machine_times, watch_data)



# OMEGA_r_mech_star = watch_data[12, :]
# plot(machine_times, OMEGA_r_mech_star, xlabel="Time (s)", ylabel="OMEGA_r_mech_star", title="OMEGA_r_mech_star vs Time")






