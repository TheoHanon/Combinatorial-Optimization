struct OptVaxData
    n::Int64
    m::Int64
    D::Array{Float64, 2}
    A::Array{Int64, 2}
    Q::Int64
    C::Vector{Int64}
    q::Vector{Int64}
    f::Vector{Int64}
    p::Int64
    B::Float64
    J_prime::Vector{Int64}
    M::Int64
end


function initOptVaxData(
    n::Int64, 
    m::Int64, 
    D::Array{Float64, 2}, 
    A::Array{Int64, 2}, 
    Q::Int64, 
    C::Vector{Int64}, 
    q::Vector{Int64}, 
    f::Vector{Int64}, 
    p::Int64, 
    B::Float64, 
    J_prime::Vector{Int64}, 
    M::Int64
)
    return OptVaxData(n, m, D, A, Q, C, q, f, p, B, J_prime, M)
end



function OptVax_Dual(
                    u_subtours::Vector{Float64}, 
                    u_capacities::Vector{Float64},
                    data::OptVaxData
                    )

        # Unpack data

        n = data.n
        m = data.m
        D = data.D
        A = data.A
        Q = data.Q
        C = data.C
        q = data.q
        f = data.f
        p = data.p
        B = data.B
        J_prime = data.J_prime
        M = data.M

        # Indices
        I = 1:m  # VC locations
        J = 1:n  # Localities
        N = 1:(n + m)  # Combined set
    
        # Initialize model
        model = Model(HiGHS.Optimizer)
        
        # Variables
        @variable(model, y[I], Bin)                           # VC locations
        @variable(model, u[J, 1:M], Bin)                     # MMT-VC assignments
        @variable(model, z[N, N, 1:M], Bin)                  # Routing
        @variable(model, v[J], Bin)                          # Coverage indicator
        @variable(model, delta[1:M], Bin)                    # MMT usage
        @variable(model, beta[N], Int)                       # Subtour elimination
    
        ############################
        #   Lagrangean relaxation  #
        ############################
        
        valid_pairs = [(i, j) for i in N, j in J if i != j]
        @expression(model, subtours, sum(u_subtours[i, j, l] * ((n - 2) - ( beta[i] - beta[j] + (n - 1) * z[i, j, k])) for k = 1:M, (i, j) in valid_pairs))
        @expression(model, capacities, sum( u_capacities[k] * (Q - sum(q[j] * z[i, j, k] for j in J, i in N)) for k in 1:M))
        @expression(model, sum_uq, sum(sum(u[:, k] for k in 1:M) .* q))

        @objective(model, Max, sum_uq + sum(q .* v) + subtours + capacities)

        ############################
        #   Constraints            #
        ############################


        @constraint(model, [j in J], v[j] == sum(A[i, j] * y[i] for i in I))
        @constraint(model, [j in J, k in 1:M], sum(z[i, j, k] for i in N) == u[j, k])
        @constraint(model, [j in J], sum(u[j, k] for k in 1:M) <= 1)
    
        # Subtour elimination constraints
        @constraint(model, beta_low[J], beta[J] .>= 2)
        @constraint(model, beta_high[J], beta[J] .<= n)
        @constraint(model, [i in I], beta[i + n] == y[i])
    
        # Total VC locations must equal 1
        @constraint(model, sum(y) == 1)
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
        
        # Solve the model
        optimize!(model)

        return model
end

