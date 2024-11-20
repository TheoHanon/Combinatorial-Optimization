using JuMP
using HiGHS
using LinearAlgebra

function solve_OptVax1(
    n::Int64, m::Int64, D::Array{Float64, 2}, A::Array{Int64, 2}, 
    Q::Int64, C::Vector{Int64}, q::Vector{Int64}, f::Vector{Int64}, 
    p::Int64, B::Float64, R::Vector{Float64}, J_prime::Vector{Int64}, M::Int64
)
    # Indices
    I = 1:m  # VC locations
    J = 1:n  # Localities
    N = 1:(n + m)  # Combined set

    # Initialize model
    model = Model(HiGHS.Optimizer)
    set_optimizer_attribute(model, "time_limit", 1800.0)

    # Variables
    @variable(model, y[I], Bin)                           # VC locations
    @variable(model, u[J, 1:M], Bin)                     # MMT-VC assignments
    @variable(model, z[N, N, 1:M], Bin)                  # Routing
    @variable(model, v[J], Bin)                          # Coverage indicator
    @variable(model, delta[1:M], Bin)                    # MMT usage
    @variable(model, beta[N], Int)                       # Subtour elimination

    # Objective: Maximize coverage and vaccination
    @expression(model, sum_uq, sum(sum(u[:, k] for k in 1:M) .* q))
    @objective(model, Max, sum_uq + sum(q .* v))

    # Constraints
    # VC coverage
    @constraint(model, [j in J], v[j] == sum(A[i, j] * y[i] for i in I))

    # Routing and assignments
    @constraint(model, [j in J, k in 1:M], sum(z[i, j, k] for i in N) == u[j, k])
    @constraint(model, [j in J], sum(u[j, k] for k in 1:M) <= 1)

    # Subtour elimination constraints
    @constraint(model, beta_low[J], beta[J] .>= 2)
    @constraint(model, [j in J], beta[j] <= n)
    
    @constraint(model, [i in I], beta[i + n] == y[i])
    valid_pairs = [(i, j) for i in N, j in J if i != j]
    @constraint(model, [k = 1:M, (i, j) in valid_pairs], 
        beta[i] - beta[j] + (n - 1) * z[i, j, k] <= n - 2)

    # Total VC locations must equal 1
    @constraint(model, sum(y) == 1)

    # Coverage and locality vaccination
    @constraint(model, sum(u[:, k] for k in 1:M) .+ v .<= 1)

    # Flow conservation constraints
    @constraint(model, [k = 1:M, h in N], 
        sum(z[i, h, k] for i in N) - sum(z[h, j, k] for j in N) == 0)
    @constraint(model, [i in N, k in 1:M], z[i, i, k] == 0)

    # Coverage and MMT-VC relationships
    @expression(model, sum_z[l in I, k in 1:M], sum(z[l + n, i, k] for i in N))
    @constraint(model, [l in I, k in 1:M], sum_z[l, k] <= y[l])
    @constraint(model, [l in I, k in 1:M], sum_z[l, k] <= delta[k])
    @constraint(model, [l in I, k in 1:M], sum_z[l, k] >= y[l] + delta[k] - 1)

    # Ensure z is only active if MMT is used
    @constraint(model, [i in N, j in N], z[i, j, :] .<= delta)
    @constraint(model, [k in 1:M, j in J], u[j, k] <= delta[k])

    # Cost constraints
    @constraint(model, sum_uq <= sum(C .* y))

    # Enforced coverage for specific localities
    @constraint(model, [j in J_prime], sum(u[j, k] for k in 1:M) + v[j] == 1)

    # Budget constraint
    @constraint(model,
        sum(D .* sum(z[:, :, k] for k in 1:M)) + 
        sum(delta) * p + 
        sum(y .* f) <= B
    )

    # Order constraints for MMT usage
    @constraint(model, [k in 1:M-1], delta[k] >= delta[k + 1])

    # Capacity constraint
    @constraint(model, [k in 1:M], 
        sum(q[j] * z[i, j, k] for j in J, i in N) <= Q)

    return model
end

function solve_OptVax2(n::Int64, m::Int64, D::Array{Float64, 2}, A::Array{Int64, 2}, Q::Int64, C::Vector{Int64}, q::Vector{Int64}, f::Vector{Int64}, p::Int64, B::Float64, R::Vector{Float64}, J_prime::Vector{Int64}, M::Int64) 
    # Indices
    I = 1:m  # VC locations
    J = 1:n  # Localities
    N = 1:(n + m)  # Combined set

    # Initialize model
    model = Model(HiGHS.Optimizer)
    set_optimizer_attribute(model, "time_limit", 1800.0)

    # Variables
    @variable(model, y[I], Bin)                           # VC locations
    @variable(model, u[J, 1:M], Bin)                     # MMT-VC assignments
    @variable(model, z[N, N, 1:M], Bin)                  # Routing
    @variable(model, v[J], Bin)                          # Coverage indicator
    @variable(model, delta[1:M], Bin)                    # MMT usage
    @variable(model, beta[N], Int)                       # Subtour elimination

    # Objective: Maximize coverage and vaccination
    @expression(model, sum_uq, sum(sum(u[:, k] for k in 1:M) .* q))
    @objective(model, Max, sum_uq + sum(q .* v))

    # Constraints
    # VC coverage
    @constraint(model, [j in J], v[j] == sum(A[i, j] * y[i] for i in I))

    # Routing and assignments
    @constraint(model, [j in J, k in 1:M], sum(z[i, j, k] for i in N) == u[j, k])
    @constraint(model, [j in J], sum(u[j, k] for k in 1:M) <= 1)

    # Subtour elimination constraints
    @constraint(model, beta_low[J], beta[J] .>= 2)
    @constraint(model, beta_high[J], beta[J] .<= n)
    @constraint(model, [i in I], beta[i + n] == y[i])
    valid_pairs = [(i, j) for i in N, j in J if i != j]
    @constraint(model, [k = 1:M, (i, j) in valid_pairs], 
        beta[i] - beta[j] + (n - 1) * z[i, j, k] <= n - 2)

    # Total VC locations must equal 1
    @constraint(model, sum(y) == 1)

    # Coverage and locality vaccination
    @constraint(model, sum(u[:, k] for k in 1:M) .+ v .<= 1)

    # Flow conservation constraints
    @constraint(model, [k = 1:M, h in N], 
        sum(z[i, h, k] for i in N) - sum(z[h, j, k] for j in N) == 0)
    @constraint(model, [i in N, k in 1:M], z[i, i, k] == 0)

    # Coverage and MMT-VC relationships
    @expression(model, sum_z[l in I, k in 1:M], sum(z[l + n, i, k] for i in N))
    @constraint(model, [l in I, k in 1:M], sum_z[l, k] <= y[l])
    @constraint(model, [l in I, k in 1:M], sum_z[l, k] <= delta[k])
    @constraint(model, [l in I, k in 1:M], sum_z[l, k] >= y[l] + delta[k] - 1)

    # Ensure z is only active if MMT is used
    @constraint(model, [i in N, j in N], z[i, j, :] .<= delta)
    @constraint(model, [k in 1:M, j in J], u[j, k] <= delta[k])

    # Cost constraints
    @constraint(model, sum_uq <= sum(C .* y))

    # Enforced coverage for specific localities
    @constraint(model, [j in J_prime], sum(u[j, k] for k in 1:M) + v[j] == 1)

    # Budget constraint
    @constraint(model,
        sum(D .* sum(z[:, :, k] for k in 1:M)) + 
        sum(delta) * p + 
        sum(y .* f) <= B
    )

    # Order constraints for MMT usage
    @constraint(model, [k in 1:M-1], delta[k] >= delta[k + 1])

    # Capacity constraint
    @constraint(model, [k in 1:M], 
        sum(q[j] * z[i, j, k] for j in J, i in N) <= Q)

    #Additional constraint
    @expression(model, sum_delta, sum(delta))
    @constraint(model, sum_delta >= sum_uq/Q)

    return model
end


function solve_OptVax2LP(n::Int64, m::Int64, D::Array{Float64, 2}, A::Array{Int64, 2}, Q::Int64, C::Vector{Int64}, q::Vector{Int64}, f::Vector{Int64}, p::Int64, B::Float64, R::Vector{Float64}, J_prime::Vector{Int64}, M::Int64) 
    I = 1:m
    J = 1:n
    N = 1:(n + m)
    
    model = Model(HiGHS.Optimizer)
    set_optimizer_attribute(model, "time_limit", 1800.0)


    # Variables
    @variable(model, 0 .<= y[I] .<= 1)
    @variable(model, 0 .<= u[J, 1:M] .<= 1)
    @variable(model, 0 .<= v[J] .<= 1)  
    @variable(model, 0 .<= z[N, N, 1:M] .<= 1)
    @variable(model, 0 .<= delta[1:M] .<= 1)
    @variable(model, beta[N])
    
    # Objective: Maximize coverage and vaccination
    @expression(model, sum_uq, sum(sum(u[:, k] for k in 1:M) .* q))
    @objective(model, Max, sum_uq + sum(q .* v))

    # Constraints
    # VC coverage
    @constraint(model, [j in J], v[j] == sum(A[i, j] * y[i] for i in I))

    # Routing and assignments
    @constraint(model, [j in J, k in 1:M], sum(z[i, j, k] for i in N) == u[j, k])
    @constraint(model, [j in J], sum(u[j, k] for k in 1:M) <= 1)

    # Subtour elimination constraints
    @constraint(model, beta_low[J], beta[J] .>= 2)
    @constraint(model, beta_high[J], beta[J] .<= n)
    @constraint(model, [i in I], beta[i + n] == y[i])
    valid_pairs = [(i, j) for i in N, j in J if i != j]
    @constraint(model, [k = 1:M, (i, j) in valid_pairs], 
        beta[i] - beta[j] + (n - 1) * z[i, j, k] <= n - 2)

    # Total VC locations must equal 1
    @constraint(model, sum(y) == 1)

    # Coverage and locality vaccination
    @constraint(model, sum(u[:, k] for k in 1:M) .+ v .<= 1)

    # Flow conservation constraints
    @constraint(model, [k = 1:M, h in N], 
        sum(z[i, h, k] for i in N) - sum(z[h, j, k] for j in N) == 0)
    @constraint(model, [i in N, k in 1:M], z[i, i, k] == 0)

    # Coverage and MMT-VC relationships
    @expression(model, sum_z[l in I, k in 1:M], sum(z[l + n, i, k] for i in N))
    @constraint(model, [l in I, k in 1:M], sum_z[l, k] <= y[l])
    @constraint(model, [l in I, k in 1:M], sum_z[l, k] <= delta[k])
    @constraint(model, [l in I, k in 1:M], sum_z[l, k] >= y[l] + delta[k] - 1)

    # Ensure z is only active if MMT is used
    @constraint(model, [i in N, j in N], z[i, j, :] .<= delta)
    @constraint(model, [k in 1:M, j in J], u[j, k] <= delta[k])

    # Cost constraints
    @constraint(model, sum_uq <= sum(C .* y))

    # Enforced coverage for specific localities
    @constraint(model, [j in J_prime], sum(u[j, k] for k in 1:M) + v[j] == 1)

    # Budget constraint
    @constraint(model,
        sum(D .* sum(z[:, :, k] for k in 1:M)) + 
        sum(delta) * p + 
        sum(y .* f) <= B
    )

    # Order constraints for MMT usage
    @constraint(model, [k in 1:M-1], delta[k] >= delta[k + 1])

    # Capacity constraint
    @constraint(model, [k in 1:M], 
        sum(q[j] * z[i, j, k] for j in J, i in N) <= Q)

    #Additional constraint
    @expression(model, sum_delta, sum(delta))
    @constraint(model, sum_delta >= sum_uq/Q)

    return model
end


