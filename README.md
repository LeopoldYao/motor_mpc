# Motor Simulator and Controller Setup

Motor and model parameters setup in `test.jl`. 
The motor simulator is in `MACHINE_SIMULATOR.jl`. 
The current controller is selected from the first line of this file.

## Controllers
1. **nestedloop_pmsm.jl** - Nested loop PI controller.
2. **mpc_pmsm_linear** - Simplified model MPC with linear dynamics.
3. **mpc_pmsm_nonlinear** - Full model MPC with non-linear dynamics, linearized around each state.

* Controller:* `mpc_pmsm_nonlinear.jl`

## Running the Simulation
Run the simulation by including the `test.jl` script in Julia REPL:
```julia
include("test.jl")
```

Alternatively, you can run the script directly from the terminal:
```shell
julia test.jl
```
