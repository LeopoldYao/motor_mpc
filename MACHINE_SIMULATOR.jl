include("mpc_pmsm_nonlinear.jl") # choose the controller to use
using LinearAlgebra
using Plots
using ProgressMeter
using PyPlot

mutable struct ACMachine
    # Name plate data
    npp::Int32
    npp_inv::Float64
    IN::Float64
    # Electrical parameters
    R::Float64
    Ld::Float64
    Lq::Float64
    KE::Float64
    Rreq::Float64
    # Mechanical parameters
    Js::Float64
    Js_inv::Float64
    # States
    NS::Int32
    x::Vector{Float64}
    # Inputs
    uab::Vector{Float64}
    udq::Vector{Float64}
    TLoad::Float64
    # Output
    omega_slip::Float64
    omega_r_elec::Float64
    omega_r_mech::Float64
    omega_syn::Float64
    theta_d::Float64
    theta_d_mech::Float64
    KA::Float64
    iD::Float64
    iQ::Float64
    iAlfa::Float64
    iBeta::Float64
    ia::Float64
    ib::Float64
    ic::Float64
    Tem::Float64
    cosT::Float64
    sinT::Float64
    # Simulation settings
    MACHINE_SIMULATIONS_PER_SAMPLING_PERIOD::Int32
    bool_apply_load_model::Int32
end
# ACMachine constructor
function The_AC_Machine(npp, IN, R, Ld, Lq, KE, Rreq, Js)
    npp_inv = 1.0 / npp
    NS = 5 # Number of states
    x = zeros(Float64, NS)
    KA = KE
    x[3] = KE # Initial condition for KA is KE
    uab = zeros(Float64, 2)
    udq = zeros(Float64, 2)
    return ACMachine(npp, npp_inv, IN, R, Ld, Lq, KE, Rreq, Js, 1.0 / Js, NS, x, uab, udq, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, KA, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1, 0)
end



# Machine dynamics
# states:
# x[1] = theta_d_mech
# x[2] = omega_r_mech
# x[3] = KA
# x[4] = iD
# x[5] = iQ
# inputs:
# u[1] = ud
# u[2] = uq
function DYNAMICS_MACHINE(t, x, ACM, CLARKE_TRANS_TORQUE_GAIN=1.5)
    fx = zeros(Float64, ACM.NS)

    theta_d_mech = x[1]
    omega_r_mech = x[2]
    ACM.theta_d = x[1] * ACM.npp
    ACM.omega_r_elec = x[2] * ACM.npp
    KA = x[3]
    iD = x[4]
    iQ = x[5]

    if KA == 0.0
        ACM.omega_slip = 0.0
    else
        ACM.omega_slip = ACM.Rreq * iQ / KA
    end

    ACM.omega_syn = x[2] * ACM.npp + ACM.omega_slip 

    # Electrical Subsystem
    if ACM.Rreq > 0
        fx[3] = ACM.Rreq * iD - ACM.Rreq / (ACM.Ld - ACM.Lq) * KA
        fx[4] = (ACM.udq[1] - ACM.R * iD + ACM.omega_syn * ACM.Lq * iQ - fx[3]) / ACM.Lq
    elseif ACM.Rreq < 0
        error("ACM.Rreq is used to calculate slip so it must be zero for PMSM.")
    else
        fx[4] = (ACM.udq[1] - ACM.R * iD + ACM.omega_syn * ACM.Lq * iQ) / ACM.Ld
        fx[3] = (ACM.Ld - ACM.Lq) * fx[4] + 0.0 # The equation has taken derivative, 0.0 is the derivative of KA
    end
    fx[5] = (ACM.udq[2] - ACM.R * iQ - ACM.omega_syn * ACM.Lq * iD - ACM.omega_syn * ACM.KA) / ACM.Lq

    # Mechanical Subsystem
    ACM.Tem = CLARKE_TRANS_TORQUE_GAIN * ACM.npp * KA * iQ
    fx[1] = x[2] + ACM.omega_slip / ACM.npp
    fx[2] = (ACM.Tem - ACM.TLoad) / ACM.Js

    return fx
end


function RK4_MACHINE(t, ACM, hs)
    NS = ACM.NS
    k1, k2, k3, k4 = zeros(Float64, NS), zeros(Float64, NS), zeros(Float64, NS), zeros(Float64, NS)
    xk, fx = zeros(Float64, NS), zeros(Float64, NS)

    fx = DYNAMICS_MACHINE(t, ACM.x, ACM)
    for i in 1:NS
        k1[i] = fx[i] * hs
        xk[i] = ACM.x[i] + k1[i] * 0.5
    end

    fx = DYNAMICS_MACHINE(t, xk, ACM)
    for i in 1:NS
        k2[i] = fx[i] * hs
        xk[i] = ACM.x[i] + k2[i] * 0.5
    end

    fx = DYNAMICS_MACHINE(t, xk, ACM)
    for i in 1:NS
        k3[i] = fx[i] * hs
        xk[i] = ACM.x[i] + k3[i]
    end

    fx = DYNAMICS_MACHINE(t, xk, ACM)
    for i in 1:NS
        k4[i] = fx[i] * hs
        ACM.x[i] = ACM.x[i] + (k1[i] + 2 * (k2[i] + k3[i]) + k4[i]) / 6.0
    end
end


# Machine simulation
function ACMSimPyIncremental(MACHINE_TS, t0, TIME, ACM, CTRL)

    machine_times = t0:MACHINE_TS:(t0 + TIME - MACHINE_TS)
    watch_data = zeros(Float64, 20, length(machine_times))
    watch_index = 1
    control_counter = 0

    progress = Progress(length(machine_times), 1)  # 1% update
    for t in machine_times
        next!(progress)

        # Numerical Integration (ode4) with 5 states
        RK4_MACHINE(t, ACM, MACHINE_TS)

        # Machine Simulation Output @ MACHINE_TS
        ACM.theta_d_mech = ACM.x[1]
        ACM.omega_r_mech = ACM.x[2]
        ACM.KA = ACM.x[3]
        ACM.iD = ACM.x[4]
        ACM.iQ = ACM.x[5]
        ACM.theta_d = ACM.theta_d_mech * ACM.npp
        ACM.omega_r_elec = ACM.omega_r_mech * ACM.npp
        ACM.omega_syn = ACM.omega_r_elec + ACM.omega_slip

        # Inverse Park transformation
        ACM.cosT = cos(ACM.theta_d)
        ACM.sinT = sin(ACM.theta_d)
        ACM.iAlfa = ACM.iD * ACM.cosT + ACM.iQ * -ACM.sinT
        ACM.iBeta = ACM.iD * ACM.sinT + ACM.iQ * ACM.cosT

        # TEST SETUP - TARGET SPEED
        # 1. Fixed Ω_star
        # OMEGA_r_mech_star = 100.0
        # 2. Change TLoad and Ω_star 
        # ACM.TLoad = t < 2 ? 0.0 : (t > 2 && t < 6 ? 3.0 : (t > 6 ? 0.0 : ACM.TLoad))
        OMEGA_r_mech_star = t < 1.5 ? 100.0 : -100.0
        # 3. Sinusoidal Ω_star
        # OMEGA_r_mech_star = 100.0 * sin(2 * π * 10 * t)

        # Calculate control action using nested PI controller
        dt = MACHINE_TS
        if control_counter % 100 == 0
            uD, uQ = control_action(CTRL, OMEGA_r_mech_star, ACM.omega_r_mech, ACM.iQ, ACM.iD, t)
            ACM.udq[1] = uD
            ACM.udq[2] = uQ
        end
        control_counter += 1
        # uQ, uD = control_action(CTRL, OMEGA_r_mech_star, ACM.omega_r_mech, ACM.iQ, ACM.iD, dt)
        # ACM.udq[1] = uD
        # ACM.udq[2] = uQ
        # if isnothing(CTRL)
        #     ACM.uab[1] = 10 * cos(2 * π * 0.5 * t)
        #     ACM.uab[2] = 10 * sin(2 * π * 0.5 * t)
        # end
        # # Park transformation
        # ACM.udq[1] = ACM.uab[1] * ACM.cosT + ACM.uab[2] * ACM.sinT
        # ACM.udq[2] = ACM.uab[1] * -ACM.sinT + ACM.uab[2] * ACM.cosT

        # # Initial condition uD
        # if t == t0
        #     uD = 0.0001
        # end

        # Watch @ MACHINE_TS
        watch_data[1, watch_index] = mod(ACM.theta_d, 2 * π)
        watch_data[2, watch_index] = ACM.omega_r_mech
        watch_data[3, watch_index] = ACM.KA
        watch_data[4, watch_index] = ACM.iD
        watch_data[5, watch_index] = ACM.iQ
        watch_data[6, watch_index] = ACM.Tem

        watch_data[7, watch_index] = ACM.udq[1] # uD
        watch_data[8, watch_index] = ACM.udq[2] # uQ
        watch_data[9, watch_index] = ACM.TLoad
        watch_data[10, watch_index] = OMEGA_r_mech_star # speed command
        
        watch_data[11, watch_index] = t

        watch_index += 1
    end

    return machine_times, watch_data
end


function myplot(machine_times, watch_data)
    ylabels = ["theta",
               "speed",
               "KA",
               "iD",
               "iQ / q-axis current",
               "Tem",
               "uD",
               "uQ / q-axis voltage",
               "TLoad",
               "OMEGA_r_mech_star",
                "time"]
        plots = []

    # Iterate over the labels and the corresponding watch_data index
    # watch_data = watch_data[:, 50000:50500]
    # machine_times = machine_times[50000:50500]


    for (index, ylabel) in enumerate(ylabels)
        trace = watch_data[index, :]
        p = plot(machine_times, trace, xlabel="Time (s)", ylabel=ylabel, title="$ylabel vs Time", legend=false, linewidth=1)
        push!(plots, p)
    end
    
    # Create a layout of 3 rows and 4 columns
    plot_com = plot(plots..., layout=(4, 4), size=(1200, 800))
    savefig( plot_com, "plots/motor_plots.svg")
    
end
