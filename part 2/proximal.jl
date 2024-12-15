include("model.jl")
using JuMP

function proximal_gradient_descent(
                                z_lower::Float64,
                                data::OptVaxData,
                                max_iter::Int64 = 1000,)

    # Unpack data

    M = data.M
    m = data.m
    n = data.n
    valid_pairs = [(i, j) for i in N, j in J if i != j]

    u_subtours = zeros(m, n, M)
    u_capacities = zeros(M)
    u = [u_subtours; u_capacities]

    g_subtours = zeros(m, n, M)
    g_capacities = zeros(M)
    
    for _ in 1:max_iter

        # Solve the dual problem
        model = OptVax_Dual(u_subtours, u_capacities, data)
        z_dual = objective_value(model)

        z = value.(model[:z])
        beta = value.(model[:beta])

        # Compute the gradient
        
        for (i, j) in valid_pairs, k in 1:M
            g_subtours[i, j, l] = (n - 2) - (beta[i] - beta[j] + (n - 1) * z[i, j, k])
        end

        for k in 1:M
            g_capacities[k] = Q - sum(q[j] * z[i, j, k] for j in J, i in N)
        end
        
        g = [vec(g_subtours); g_capacities]

        # Compute the step size

        alpha = (z_dual - z_lower) / (norm(g)^2)

        # Update the variables

        u = max(0, u + alpha * g)
        u_subtours = reshape(u[1:(m * n * M)], m, n, M)
        u_capacities = u[(m * n * M + 1):end]

    end

    return z_dual
end