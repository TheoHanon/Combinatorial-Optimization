using JuMP
using HiGHS
using LinearAlgebra

function solve_OptVax1(n::Int64, m::Int64, D::Array{Float64, 2}, A::Array{Int64, 2}, Q::Int64, C::Vector{Int64}, q::Vector{Int64}, f::Vector{Int64}, p::Int64, B::Float64, R::Vector{Float64}, J_prime::Vector{Int64}, M::Int64)
    I = 1:m
    J = 1:n
    N = 1:(n + m)
    
    model = Model(HiGHS.Optimizer)

    # Variables
    @variable(model, y[I], Bin)
    @variable(model, u[J], Bin)
    #@variable(model, v[J], Bin)
    @variable(model, z[N, N, 1:M], Bin)
    @variable(model, delta[1:M], Bin)
    @variable(model, beta[N], Int)
    @expression(model, v[j in J], dot(A[I, j], y[I]))
    #Expression to avoid multiple computations
    @expression(model, sum_uq, sum(u .* q))

    # Objective function
    @objective(model, Max, sum_uq + sum(q .* v))

    # Constraints

    # Subtour elimination constraints
    @constraint(model, beta[J] .>= q)
    @constraint(model, beta[J] .<= Q)
    valid_pairs = [(i, j) for i in N, j in J if i != j]
    @constraint(model, [k = 1:M, (i, j) in valid_pairs], beta[i] - beta[j] >= q[j] - (1 - z[i, j, k]) * Q)

    # Total VC locations must equal 1
    @constraint(model, sum(y) == 1)

    # Coverage and locality vaccination relationships
    #@constraint(model, [j in J], v[j] == sum(A[i, j] * y[i] for i in I))
    @constraint(model, u .+ v .<= 1.0)

    # Routing constraints
    @constraint(model, [j in J], sum(z[i, j, k] for k in 1:M, i in N) == u[j])

    # Flow conservation constraints
    @constraint(model, [k=1:M, h in N], sum(z[i, h, k] for i in N) - sum(z[h, j, k] for j in N) == 0.0)

    # Avoid double arrows
    @constraint(model, [i in N, k in 1:M], z[i, i, k] == 0)

    # Coverage and MMT-VC relationships
    @expression(model, sum_z[l in I, k in 1:M], sum(z[l + n, i, k] for i in N))
    @constraint(model, [l in I, k in 1:M], sum_z[l, k] <= y[l])
    @constraint(model, [l in I, k in 1:M], sum_z[l, k] <= delta[k])
    @constraint(model, [l in I, k in 1:M], sum_z[l, k] >= y[l] + delta[k] - 1)

    # Ensure z is only active if MMT is used
    @constraint(model, [i in N, j in N], z[i, j, :] .<= delta)

    # Cost constraint
    @constraint(model, sum_uq <= sum(C .* y))

    # Enforced coverage for specific localities
    @constraint(model, [j in J_prime], u[j] + v[j] == 1)

    # Budget constraint
    @constraint(model,
        sum(D .* sum(z[:,:,k] for k in 1:M)) +
        sum(delta) * p +
        sum(y .* f) <= B
    )

    optimize!(model)
    return model
end



function solve_OptVax2(n::Int64, m::Int64, D::Array{Float64, 2}, A::Array{Int64, 2}, Q::Int64, C::Vector{Int64}, q::Vector{Int64}, f::Vector{Int64}, p::Int64, B::Float64, R::Vector{Float64}, J_prime::Vector{Int64}, M::Int64) 
    I = 1:m
    J = 1:n
    N = 1:(n + m)
    
    model = Model(HiGHS.Optimizer)

    # Variables
    @variable(model, y[I], Bin)
    @variable(model, u[J], Bin)
    #@variable(model, v[J], Bin)
    @variable(model, z[N, N, 1:M], Bin)
    @variable(model, delta[1:M], Bin)
    @variable(model, beta[N], Int)
    @expression(model, v[j in J], dot(A[I, j], y[I]))
    #Expression to avoid multiple computations
    @expression(model, sum_uq, sum(u .* q))

    # Objective function
    @objective(model, Max, sum_uq + sum(q .* v))

    # Constraints

    # Subtour elimination constraints
    @constraint(model, beta[J] .>= q)
    @constraint(model, beta[J] .<= Q)
    valid_pairs = [(i, j) for i in N, j in J if i != j]
    @constraint(model, [k = 1:M, (i, j) in valid_pairs], beta[i] - beta[j] >= q[j] - (1 - z[i, j, k]) * Q)

    # Total VC locations must equal 1
    @constraint(model, sum(y) == 1)

    # Coverage and locality vaccination relationships
    #@constraint(model, [j in J], v[j] == sum(A[i, j] * y[i] for i in I))
    @constraint(model, u .+ v .<= 1.0)

    # Routing constraints
    @constraint(model, [j in J], sum(z[i, j, k] for k in 1:M, i in N) == u[j])

    # Flow conservation constraints
    @constraint(model, [k=1:M, h in N], sum(z[i, h, k] for i in N) - sum(z[h, j, k] for j in N) == 0.0)

    # Avoid double arrows
    @constraint(model, [i in N, k in 1:M], z[i, i, k] == 0)

    # Coverage and MMT-VC relationships
    @expression(model, sum_z[l in I, k in 1:M], sum(z[l + n, i, k] for i in N))
    @constraint(model, [l in I, k in 1:M], sum_z[l, k] <= y[l])
    @constraint(model, [l in I, k in 1:M], sum_z[l, k] <= delta[k])
    @constraint(model, [l in I, k in 1:M], sum_z[l, k] >= y[l] + delta[k] - 1)

    # Ensure z is only active if MMT is used
    @constraint(model, [i in N, j in N, k in 1:M], z[i, j, :] .<= delta)

    # Cost constraint
    @constraint(model, sum_uq <= sum(C .* y))

    # Enforced coverage for specific localities
    @constraint(model, [j in J_prime], u[j] + v[j] == 1)

    @expression(model, sum_delta, sum(delta))

    # Budget constraint
    @constraint(model,
        sum(D .* sum(z[:,:,k] for k in 1:M)) +
        sum_delta * p +
        sum(y .* f) <= B
    )

    #Additional constraint
    @constraint(model, sum_delta >= sum_uq/Q)


    optimize!(model)
    return model
end


