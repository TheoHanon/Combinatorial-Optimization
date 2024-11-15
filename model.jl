using JuMP
using HiGHS

function solve_OptVax1(n::Int64, m::Int64, D::Array{Float64, 2}, A::Array{Int64, 2}, Q::Int64, C::Vector{Int64}, q::Vector{Int64}, f::Vector{Int64}, p::Int64, B::Float64, R::Vector{Float64}, J_prime::Vector{Int64}, M::Int64) 
    
    I = 1:m
    J = 1:n
    N = 1:(n+m)
 
    model = Model(HiGHS.Optimizer)

    # Creation of the variables
    
    @variable(model, y[i in I], Bin)
    @variable(model, u[j in J], Bin)
    @variable(model, v[j in J], Bin)
    @variable(model, z[i in N, j in N, k in 1:M], Bin)
    @variable(model, delta[k in 1:M], Bin)

    

    # Objective function
    @objective(model, Max, sum(q[j] * (u[j] + v[j]) for j in J))

    # Constraints


    ## SUBTOUR ELIMINATION CONSTRAINTS---------------
    @variable(model, beta[i in N], Int)

    for j in J
        @constraint(model, beta[j] >= q[j])
        @constraint(model, beta[j] <= Q)
    end

    for k in 1:M, i in N, j in J
        if i != j
            @constraint(model, beta[i] - beta[j] >=  q[j] - (1-z[i,j,k]) * Q)
        end
    end

    #------------------------------------------------

    # 1. Total VC locations must equal 1
    @constraint(model, sum(y[i] for i in I) == 1)
    
    # 2. Coverage and locality vaccination relationships
    for j in J
        @constraint(model, v[j] == sum(A[i, j] * y[i] for i in I))
        @constraint(model, u[j] + v[j] <=1.0)
    end

    #  3. Routing constraints
    for j in J 
        @constraint(model, sum(z[i, j, k] for k in 1:M, i in N) == u[j])
    end

    # 4.  Flow conservation constraints
    for k in 1:M
        for h in N
            @constraint(model, sum(z[i, h, k] for i in N) - sum(z[h, j, k] for j in N) == 0.0)
        end
    end

    # 5. Capacity constraints for MMTs
    # for k in 1:M
    #     @constraint(model, sum(q[j] * z[i,j,k] for j in J, i in N) <= Q)
    # end

    # 5.5. Avoid double arrows
    for i in N, k in 1:M            
        @constraint(model, z[i,i,k] == 0)
    end

    # 6. Coverage and MMT-VC relationships
    for l in I
        for k in 1:M
            @constraint(model, sum(z[l+n, i, k] for i in N) <= y[l])
            @constraint(model, sum(z[l+n, i, k] for i in N) <= delta[k])
            @constraint(model, sum(z[l+n, i, k] for i in N) >= y[l] + delta[k] - 1)
        end
    end

    # 7. Ensure z is only active if MMT is used
    for i in N, j in N, k in 1:M
        @constraint(model, z[i, j, k] <= delta[k])
    end

    # 8. Cost constraint
    @constraint(model, sum(q[j] * u[j] for j in J) <= sum(C[i] * y[i] for i in I))

    # 9. Enforced coverage for specific localities
    for j in J_prime
        @constraint(model, u[j] + v[j] == 1)
    end

    # 10. Budget constraint
    @constraint(model, (sum(z[i,j,k]*D[i,j] for k in 1:M, i in N, j in N)) + (sum(delta[k] for k in 1:M) * p) + sum(y[i] * f[i] for i in I) <= B)

    optimize!(model)

    return model
end


function solve_OptVax2(n::Int64, m::Int64, D::Array{Float64, 2}, A::Array{Int64, 2}, Q::Int64, C::Vector{Int64}, q::Vector{Int64}, f::Vector{Int64}, p::Int64, B::Float64, R::Vector{Float64}, J_prime::Vector{Int64}, M::Int64) 
    
    I = 1:m
    J = 1:n
    N = 1:(n+m)
 
    model = Model(HiGHS.Optimizer)

    # Creation of the variables
    
    @variable(model, y[i in I], Bin)
    @variable(model, x[j in J], Bin)
    @variable(model, u[j in J], Bin)
    @variable(model, v[j in J], Bin)
    @variable(model, z[i in N, j in N, k in 1:M], Bin)
    @variable(model, delta[k in 1:M], Bin)

    # Objective function

    @objective(model, Max, sum(q[j] * x[j] for j in J))

    # Constraints

    ## SUBTOUR ELIMINATION CONSTRAINTS---------------
    @variable(model, beta[i in N], Int)

    for i in I
        @constraint(model, beta[i+n] == y[i])
    end

    for j in J
        @constraint(model, beta[j] >=2)
        @constraint(model, beta[j] <= n)
    end

    for k in 1:M, i in N, j in J
        if i != j
            @constraint(model, beta[i] - beta[j] + (n-1)*z[i,j,k] <= n-2)
        end
    end

    #------------------------------------------------

    # 1. Total VC locations must equal 1
    @constraint(model, sum(y[i] for i in I) == 1)
    
    # 2. Coverage and locality vaccination relationships
    for j in J
        @constraint(model, v[j] == sum(A[i, j] * y[i] for i in I))
        @constraint(model, x[j] <= u[j] + v[j])
        @constraint(model, u[j] + v[j] <=1.0)
    end

    #  3. Routing constraints
    for j in J 
        @constraint(model, sum(z[i, j, k] for k in 1:M, i in N) == u[j])
    end

    # 4.  Flow conservation constraints
    for k in 1:M
        for h in N
            @constraint(model, sum(z[i, h, k] for i in N) - sum(z[h, j, k] for j in N) == 0.0)
        end
    end

    # 5. Capacity constraints for MMTs
    for k in 1:M
        @constraint(model, sum(q[j] * z[i,j,k] for j in J, i in N) <= Q)
    end

    # 5.5. Avoid double arrows
    for i in N, k in 1:M            
        @constraint(model, z[i,i,k] == 0)
    end

    # 6. Coverage and MMT-VC relationships
    for l in I
        for k in 1:M
            @constraint(model, sum(z[l+n, i, k] for i in N) <= y[l])
            @constraint(model, sum(z[l+n, i, k] for i in N) <= delta[k])
            @constraint(model, sum(z[l+n, i, k] for i in N) >= y[l] + delta[k] - 1)
        end
    end

    # 7. Ensure z is only active if MMT is used
    for i in N, j in N, k in 1:M
        @constraint(model, z[i, j, k] <= delta[k])
    end

    # 8. Cost constraint
    @constraint(model, sum(q[j] * u[j] for j in J) <= sum(C[i] * y[i] for i in I))

    # 9. Enforced coverage for specific localities
    for j in J_prime
        @constraint(model, x[j] == 1)
    end

    # 10. Budget constraint
    @constraint(model, (sum(z[i,j,k]*D[i,j] for k in 1:M, i in N, j in N)) + (sum(delta[k] for k in 1:M) * p) + sum(y[i] * f[i] for i in I) <= B)


    optimize!(model)

    return model
end


