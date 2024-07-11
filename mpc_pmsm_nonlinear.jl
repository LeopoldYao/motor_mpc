using JuMP
using OSQP
using LinearAlgebra
using Plots
using Printf
using Ipopt
using KNITRO





mutable struct MPControl
    # MPC parameters (fixed)
    dt::Float64 # Sampling time
    N::Int # Prediction horizon
    nx::Int # Number of states
    nu::Int # Number of control inputs

    # Cost function (fixed)
    Q::Matrix{Float64} # State error weight
    R::Matrix{Float64} # Control input 
    F::Matrix{Float64} # Terminal cost
    # Constraints (fixed)
    u_max::Vector{Float64}
    

    # System matrix (A[1,2], A[2,1] change during MPCSim)
    A::Matrix{Float64} 
    B::Matrix{Float64}
    Ka12::Float64 # A[1,2] = Ka12 * Ω_star
    Ka21::Float64 # A[2,1] = Ka21 * Ω_star
    Ka32_id::Float64 # A[3,2] = A[3,2]' * (1 + Ka32_id * iD)  # ΨA = ΨPM * (1 + Ka32_id * iD)
    K_QL::Float64 # T_L * K_QL = iQL
    # C::Matrix{Float64}

    # Intermediate variables(Ω_star, T_L change during MPCSim)
    Ω_star::Float64
    T_L::Float64


    # Initial state & Reference
    # x[1] = iD
    # x[2] = iQ
    # x[3] = Ω
    # x[4] = uD[k-1]
    # x[5] = uQ[k-1]
    # x[6] = T_L
    # u[1] = ΔuD
    # u[2] = ΔuQ
    x_0::Vector{Float64} # current motor state, initial state of a MPC roll out, [x[1], x[2], x[3], x[4], x[5], x[6]]
    u_pre::Vector{Float64} # previous motor input, [uD[k-1], uQ[k-1]]
    # x_ref::Vector{Float64}



    function MPControl(dt::Float64, Ld::Float64, Lq::Float64, R::Float64, npp::Int, ΨPM::Float64, ΨA::Float64, Js::Float64, Ω_star::Float64, T_L::Float64,
        N::Int, Q1::Float64, Q2::Float64, Q3::Float64, R1::Float64, R2::Float64,  F1::Float64, F2::Float64, 
        uD_max::Float64=10.0, uQ_max::Float64=10.0) 
        """
        constructor for MPControl
        """
        A = [
            (Ld - R * dt) / Ld    dt * (Lq / Ld) * npp * Ω_star  0.0                         dt / Ld  0.0       0.0;
            -dt * npp * Ω_star    (Lq - R * dt) / Lq             -dt * npp * ΨPM / Lq        0.0      dt / Lq   0.0;
            0.0                   dt * (3/2) * npp * ΨA / Js      1.0                         0.0      0.0       -dt / Js;
            0.0                   0.0                            0.0                         1.0      0.0       0.0;
            0.0                   0.0                            0.0                         0.0      1.0       0.0;
            0.0                   0.0                            0.0                         0.0      0.0       1.0
        ] # set ΨA_0 = ΨPM
        Ka12 = dt * (Lq / Ld) * npp
        Ka21 = -dt * npp
        Ka32_id = (Ld - Lq) / ΨPM 
        K_QL = 1 / (3/2 * npp * ΨPM) # For equavalent TLoad inverse current
        B = [
            dt / Ld  0.0;
            0.0      dt / Lq;
            0.0      0.0;
            1.0      0.0;
            0.0      1.0;
            0.0      0.0
        ]

        Q = Diagonal([Q1, Q2, Q3]) # Q[1,1] for iD; Q[2,2] for ( iQ - iQL ); Q[3,3] for (Ω - Ω_star). where iQL = T_L / ( 3/2 * npp * ΨPM) = T_L * K_QL
        R = Diagonal([R1, R2]) # R[1,1] for uD ; R[2,2] for uQ
        F = Diagonal([F1, F2]) # F[1,1] for iD; F[2,2] for Ω
        u_max = [uD_max, uQ_max]

        Ω_star = 0.0
        T_L = 0.0
        x_0 = [0.0, 0.0, 0.0, 0.0, 0.0, 0,0]
        u_pre = [0.0, 0.0]

        nx = size(A, 1)
        nu = size(B, 2)


        new(dt, N, nx, nu, 
            Q, R, F, u_max, 
            A, B, Ka12, Ka21, Ka32_id, K_QL,
            Ω_star, T_L,
            x_0, u_pre)
    end

end


# function control_action(Ad, Bd, Cd, Q, R, N, x_ref, x0, u_max)
function control_action(CTRL::MPControl, Ω_star::Float64, Ω::Float64, iQ::Float64, iD::Float64, t::Float64)
    """
    To calculate the control action uD, uQ using MPC
    Target speed Ω_star is given
    """
    # Update MPControl's content (with new Ω_star)
    CTRL.Ω_star = Ω_star
    CTRL.A[1,2], CTRL.A[2,1] = CTRL.Ka12 * Ω_star, CTRL.Ka21 * Ω_star
    CTRL.x_0 = [iD, iQ, Ω, CTRL.u_pre[1], CTRL.u_pre[2], CTRL.T_L]
    

    # Create a non-linear optimization model
    model = Model(Ipopt.Optimizer)
    set_silent(model)

    # Define variables x, u
    @variable(model, x[1:CTRL.nx, 1:CTRL.N+1])
    @variable(model, u[1:CTRL.nu, 1:CTRL.N])

    # Initial state constraint
    @constraint(model, x[:, 1] .== CTRL.x_0) 
    # A. Simplified dynamics constraints
    # for k in 1:CTRL.N
    #     @constraint(model, x[:, k+1] .== CTRL.A * x[:, k] + CTRL.B * u[:, k])
    # end
    # B. Non-linear dynamics constraints
    for k in 1:CTRL.N
        @NLconstraint(model, x[1, k+1] == CTRL.A[1,1]         * x[1, k] + CTRL.Ka12 * x[3, k] * x[2, k]  + CTRL.A[1,4] * x[4, k]                            + CTRL.B[1,1] * u[1, k] )
        @NLconstraint(model, x[2, k+1] == CTRL.Ka21 * x[3, k] * x[1, k] + CTRL.A[2,2]         * x[2, k]  + CTRL.A[2,3] * x[3, k] + CTRL.A[2,5] * x[5, k]    + CTRL.B[2,2] * u[2, k] )
        @NLconstraint(model, x[3, k+1] == x[3, k] + CTRL.A[3, 2] * (1 + CTRL.Ka32_id * x[1, k])  * x[2, k] + CTRL.A[3,6] * x[6, k])
        for i in [4, 5, 6]
            @constraint(model, x[i, k+1] == dot(CTRL.A[i, :], x[:, k]) + dot(CTRL.B[i, :], u[:, k]))
        end
    end

    # Input and state process constraints
    # @constraint(model, u .<= 100.0)
    # @constraint(model, u .>= -100.0) # Smooth
    # # uD, uQ constraints
    @constraint(model, x[4:5, 1:CTRL.N] .<= reshape(CTRL.u_max, 2, 1))
    @constraint(model, x[4:5, 1:CTRL.N] .>= -reshape(CTRL.u_max, 2, 1)) # Assuming symmetry in constraints
    # # T_L constraint
    # @constraint(model, x[6, 1:CTRL.N+1] .<= 0.1)
    # @constraint(model, x[6, 1:CTRL.N+1] .>= -0.1)

    # Output and state target constraints
    # iD constraint in last 50 of N steps
    @constraint(model, x[1, CTRL.N-50 :CTRL.N] .<= 5.0) 
    @constraint(model, x[1, CTRL.N-50 :CTRL.N] .>= -5.0) 
    # Ω constraint in last 50 of N steps
    @constraint(model, x[3, CTRL.N-50 :CTRL.N] .<= Ω_star + 5.0)
    @constraint(model, x[3, CTRL.N-50 :CTRL.N] .>= Ω_star - 5.0)



    # Objective: Minimize cost over trajectory
    @NLobjective(model, Min,    sum(CTRL.Q[1,1] * (x[1, k])^2 + 
                                    CTRL.Q[3,3] * (x[3, k] - Ω_star)^2 for k in 1:CTRL.N+1) +
                                sum(CTRL.R[1,1] * (u[1, k])^2 + CTRL.R[2,2] * (u[2, k])^2 for k in 1:CTRL.N) +
                                CTRL.F[1,1] * (x[3, CTRL.N+1] - Ω_star)^2 + CTRL.F[2,2] * (x[1, CTRL.N+1])^2)
    
    # Solve the optimization problem
    set_optimizer_attribute(model, "max_iter", 100000)
    optimize!(model)

    # Extract control action
    u1, u2 = value.(u[:, 1]) # u1 is ΔuD; u2 is ΔuQ
    uD, uQ = CTRL.u_pre[1] + u1, CTRL.u_pre[2] + u2

    CTRL.u_pre = [uD, uQ]

    # Plot 
    if t == 1.0
        ylabels = ["iD",
               "iQ",
               "speed",
               "uD_pre",
               "uQ_pre",
               "TLoad",
                "ΔuD",
                "ΔuQ"]
        plots = []
        for (index, ylabel) in enumerate(ylabels[1:6])
            trace = [value(x[index, k]) for k in 1:CTRL.N+1]
            # trace = [value(x[index, k]) for k in 1:10] # Prediction from t = 1.0 to t = 1.0 + 10*0.001 = 1.01
            p = plot( trace, xlabel="Time (s)", ylabel=ylabel, title="$ylabel vs Time", legend=false, linewidth=2)
            push!(plots, p)
        end
        for (index, ylabel) in enumerate(ylabels[7:8])
            trace = [value(u[index, k]) for k in 1:CTRL.N]
            # trace = [value(u[index, k]) for k in 1:9]
            p = plot( trace, xlabel="Time (s)", ylabel=ylabel, title="$ylabel vs Time", legend=false, linewidth=2)
            push!(plots, p)
        end
        plot(plots..., layout=(3, 3), size=(1600, 1600))
        # savefig("mpc_plots_short.png")
        savefig("plots/mpc_plots.png")
    end

    # # print the control action
    # print("Control action delta1: ", value.(u[1, 1]), " ", value.(u[2, 1]), "\n")
    # print("Control action delta2: ", value.(u[1, 2]), " ", value.(u[2, 2]), "\n")
    # print("Control action: ", uD, " ", uQ, "\n")
    return uD, uQ
end





















