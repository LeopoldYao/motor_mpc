using LinearAlgebra

# System dynamics model
function dynamics(x, u, dt)
    A = [1.0 dt; 0 1.0]
    B = [0.5*dt^2; dt]
    return A * x + B * u[1]
end

# Stage cost function
function stage_cost(x, u, Q, R)
    return 0.5 * (x' * Q * x + u' * R * u)
end

# Terminal cost function
function terminal_cost(x, Qf)
    return 0.5 * x' * Qf * x
end

# iLQR algorithm implementation
function ilqr(x0, U, dt, N, Q, R, Qf, max_iter)
    n = length(x0)  # State dimension
    m = length(U[1])  # Control input dimension

    X = [x0]  # Initialize state trajectory
    for k in 1:N
        push!(X, dynamics(X[end], U[k], dt))
    end

    for iter in 1:max_iter
        # Linearization and quadratic approximation
        A = Array{Matrix{Float64}}(undef, N)
        B = Array{Matrix{Float64}}(undef, N)
        L = Array{Matrix{Float64}}(undef, N)
        Vx = Qf * X[end] # Terminal state cost Qf is [n, n], X[end] is [n, 1] -> Vx is [n, 1]
        Vxx = Qf # Terminal state cost Qf is [n, n] -> Vxx is [n, n]

        for k in N:-1:1
            A[k] = [1.0 dt; 0 1.0] # A[k] is [n, n]
            B[k] = [0.5*dt^2 dt]' # B[k] is [n, m]
            Qx = Q * X[k] + A[k]' * Vx # Q is [n, n], X[k] is [n, 1] -> Qx is [n, 1]
            Qu = R * U[k] + B[k]' * Vx # R is a scalar [m, m], U[k] is [m, 1] -> Qu is [m, 1]
            Qxx = Q + A[k]' * Vxx * A[k] # Q is [n, n], A[k] is [n, n], Vxx is [n, n] -> Qxx is [n, n]
            Quu = R + B[k]' * Vxx * B[k] # R is a scalar [m, m], B[k] is [n, m], Vxx is [n, n] -> Quu is [m, m]
            Qux = B[k]' * Vxx * A[k] # B[k] is [n, m], Vxx is [n, n], A[k] is [n, n] -> Qux is [m, n]

            K = -inv(Quu) * Qux # Quu is [m, m], Qux is [m, n] -> K is [m, n]
            k_vec = -inv(Quu) * Qu # Quu is [m, m], Qu is [m, 1] -> k_vec is [m, 1]

            L[k] = K
            U[k] += k_vec

            Vx = Qx + K' * Quu * k_vec + K' * Qu + Qux' * k_vec
            Vxx = Qxx + K' * Quu * K + K' * Qux + Qux' * K
        end

        # Update state trajectory
        X = [x0]
        for k in 1:N
            x_next = dynamics(X[end], U[k], dt)
            print("x_next: ", x_next)
            push!(X, x_next)
        end
    end

    return X, U
end

# Optimize control inputs using iLQR and simulate the system
function main()
    x0 = [0.0; 0.0]  # Initial state
    N = 50  # Number of time steps
    dt = 0.1  # Time step size
    max_iter = 100  # Maximum number of iterations

    U = [[0.1 * randn()] for _ in 1:N]  # Initial control inputs

    Q = Diagonal([1.0, 1.0])
    R = Diagonal([0.1])
    Qf = Diagonal([10.0, 10.0])

    X, U_opt = ilqr(x0, U, dt, N, Q, R, Qf, max_iter)

    # Plot results
    x1_values = [x[1] for x in X]
    x2_values = [x[2] for x in X]
    # Plot x1
    p1 = plot(x1_values, label="x1", xlabel="Time", ylabel="x1", title="State x1 vs. Time")
    # Plot x2
    p2 = plot(x2_values, label="x2", xlabel="Time", ylabel="x2", title="State x2 vs. Time")
    # Plot U
    p3 = plot(U_opt, label="u", xlabel="Time", ylabel="u", title="Control input u vs. Time")

    plot(p1, p2, p3, layout=(3, 1))

    # println("Optimized state trajectory:")
    # println(X)

    # println("Optimized control inputs:")
    # println(U_opt)
end

main()
