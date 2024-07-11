include("mpc_pmsm_nonlinear.jl") # choose the controller to use
using JuMP
using OSQP
using LinearAlgebra
using Plots
using Printf
using Ipopt

struct ACMachine
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


npp,   IN,    R,    Ld,    Lq,  KE, Rreq,   Js = (
  2, 10.0, 10.0, 0.015, 0.025, 0.5,  0.0, 0.05)
global ACM = The_AC_Machine(npp, IN, R, Ld, Lq, KE, Rreq, Js)


ACM.theta_d = 2.0