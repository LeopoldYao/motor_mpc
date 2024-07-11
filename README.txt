# Run include("test.jl") in Julia REPL

# Motor and model parameters are setup in "test.jl"

# "MACHINE_SIMULATOR.jl" is motor simulator

# Currently we use "mpc_pmsm_nonlinear.jl" as controller, which is chosen from first line of MACHINE_SIMULATOR.jl



# Controller List
# 1. nestedloop_pmsm.jl is nestedloop PI controller
# 2. mpc_pmsm_linear is simplified model mpc, with linear problem. It setups a linear model dynammics
# 3. mpc_pmsm_nonlinear is full model mpc, with non-linear problem. It linearizes the non-linear model dynammics around each state.

# Plots
# "mpc_plots.png" from script mpc_pmsm.jl shows -- one model predictive roll out traj at time t=1s.
# "motor_plots.png" from function myplot and output "Julia Plots" shows -- real simulator's states over simulation TIME.