using JuMP
using HiGHS

function solve_OptVax1(n::Int64, m::Int64, D::Array{Float64, 2}, A::Array{Int64, 2}, Q::Int64, C::Vector{Int64}, q::Vector{Int64}, f::Vector{Int64}, p::Int64, B::Float64, R::Vector{Float64}, J_prime::Vector{Int64}, M::Int64) 
    
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
    for i in N
        for j in N
            for k in 1:M
                @constraint(model, z[i,j,k]+z[j,i,k] <= 1)
            end
        end
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

    # 11. Total MMT limit
    @constraint(model, sum(delta[k] for k in 1:M) <= M)

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

    # 1. Total VC locations must equal 1
    @constraint(model, sum(y[i] for i in I) == 1)
    
    # 2. Coverage and locality vaccination relationships
    for j in J
        @constraint(model, v[j] == sum(A[i, j] * y[i] for i in I))
        @constraint(model, x[j] <= u[j] + v[j])
        @constraint(model, u[j] + v[j] <=1)

    end

    #  3. Routing constraints
    for j in J
        @constraint(model, sum(z[i, j, k] for k in 1:M, i in N) == u[j])
    end

    # 4.  Flow conservation constraints
    for k in 1:M
        for h in N
            @constraint(model, sum(z[i, h, k] for i in N) - sum(z[h, j, k] for j in N) == 0)
        end
    end

    # 5. Capacity constraints for MMTs
    for k in 1:M
        @constraint(model, sum(q[j] * z[i,j,k] for j in J, i in N) <= Q)
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

    # 11. Total MMT limit
    @constraint(model, sum(delta[k] for k in 1:M) <= M)

    # 12. Additionnal constraints on the number of MMTs
    
    @constraint(model, sum(delta[k] for k in 1:M) >= sum(u[j]*q[j]/Q for j in J))

    optimize!(model)

    return model
end


