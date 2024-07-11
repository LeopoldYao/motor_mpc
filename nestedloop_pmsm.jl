using LinearAlgebra

mutable struct NestedLoopsController
    Kp_i::Float64
    Ki_i::Float64
    Kp_v::Float64
    Ki_v::Float64
    int_error_i::Float64
    int_error_v::Float64
    Kp_id::Float64
    Ki_id::Float64
    int_error_id::Float64
    function NestedLoopsController(Kp_i, Ki_i, Kp_v, Ki_v, Kp_id, Ki_id)
        new(Kp_i, Ki_i, Kp_v, Ki_v, 0.0, 0.0, Kp_id, Ki_id, 0.0001)
    end
end

function control_action(controller::NestedLoopsController, omega_r_mech_star::Float64, omega_r_mech::Float64, iQ::Float64, iD::Float64, dt::Float64)
    # Speed control error
    error_v = omega_r_mech_star - omega_r_mech
    controller.int_error_v += error_v * dt
    iQ_star = controller.Kp_v * (error_v + controller.Ki_v * controller.int_error_v)  # Desired iQ from speed control
    error_i = iQ_star - iQ
    controller.int_error_i += error_i * dt
    uQ = controller.Kp_i * (error_i + controller.Ki_i * controller.int_error_i)  # Compute q-axis voltage command
    
    # D-axis Current control error
    iD_star = 0
    # Current control error
    error_iD = iD_star - iD
    controller.int_error_id += error_iD * dt
    uD = controller.Kp_id * (error_iD + controller.Ki_id * controller.int_error_id)
    # controller.last_omega_r_mech = omega_r_mech
    return uQ, uD
end



