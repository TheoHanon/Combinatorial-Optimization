using JuMP


function local_search_2Opt(MMTs, D)

    new_MMTs = deepcopy(MMTs)

    for (idx, MMT_route) in enumerate(new_MMTs)

        best_route = [i for (i, j) in MMT_route[2:end]] 
        push!(best_route, MMT_route[end][2])

        best_cost = sum(D[best_route[i], best_route[i+1]] for i in 1:length(best_route)-1)
        improved = true

        while improved
            improved = false
            for i in 2:length(best_route)-2
                for j in i+1:length(best_route)-1
                    new_route = deepcopy(best_route)
                    new_route[i:j] = reverse(new_route[i:j])
                    new_cost = sum(D[new_route[i], new_route[i+1]] for i in 1:length(new_route)-1)
                    if new_cost < best_cost 
                        best_route = deepcopy(new_route)
                        best_cost = new_cost
                        improved = true
                    end
                end
            end
        end

        new_MMTs[idx] = [(best_route[i], best_route[i+1]) for i in 1:length(best_route)-1]
        pushfirst!(new_MMTs[idx], (best_route[1], best_route[1]))
    end

    return new_MMTs
end

function greedy_MMTs(VC::Int64, n::Int64, m::Int64, D::Array{Float64, 2}, A::Array{Int64, 2}, Q::Int64, C::Vector{Int64}, q::Vector{Int64}, f::Vector{Int64}, p::Int64, B::Float64, J_prime::Vector{Int64}, M::Int64)

    MMTs = [[(VC+n, VC+n)]]
    Q_MMTs = [0]

    Budget = p + f[VC]
    Q_tot = 0
    
    localities_to_visit = [j for j in 1:n if A[VC, j] == 0]
    
    while true
        
        assigned = false
        utilities = [(q[j] / D[MMT[end][2], j], idx, j) for (idx, MMT) in enumerate(MMTs) for j in localities_to_visit if !(j in J_prime)]
        utilities_prio = [(q[j]/D[MMT[end][2], j], idx, j) for (idx, MMT) in enumerate(MMTs) for j in localities_to_visit if j in J_prime]

        sorted_utilities = vcat(sort(utilities_prio, by = x -> -x[1]), sort(utilities, by = x -> -x[1])) # Put high priority localities first

        for (utility, idx, j) in sorted_utilities
        
            if (Q_MMTs[idx] + q[j] <= Q) && (Budget + D[MMTs[idx][end][2], j] + D[j, VC+n] <= B) && (Q_tot + q[j] <= C[VC])
                push!(MMTs[idx], (MMTs[idx][end][2], j))
                Q_MMTs[idx] += q[j]
                Budget += D[MMTs[idx][end][2], j]
                Q_tot += q[j]
                deleteat!(localities_to_visit, findfirst(==(j), localities_to_visit))
                assigned = true
                break

            elseif (q[j] <= Q) && (Budget + D[VC+n, j] + D[j, VC+n] + p <= B) && (Q_tot + q[j] <= C[VC]) && (length(MMTs) < M)
                push!(MMTs, [(VC+n, VC+n), (VC+n, j)])
                push!(Q_MMTs, q[j])
                Budget += D[VC+n, j] + p
                Q_tot += q[j]
                deleteat!(localities_to_visit, findfirst(==(j), localities_to_visit))
                assigned = true
                break
            end
    
        end

        if !assigned
            break
        end
    end

    for MMT in MMTs
        push!(MMT, (MMT[end][2], VC+n)) # Return to VC
        Budget += D[MMT[end][2], VC+n] # Add it to the budget; no need to check if it exceeds B already checked in the loop
    end
    
    return MMTs, Budget, Q_MMTs, Q_tot
end


function greedy_OptVax(n::Int64, m::Int64, D::Array{Float64, 2}, A::Array{Int64, 2}, Q::Int64, C::Vector{Int64}, q::Vector{Int64}, f::Vector{Int64}, p::Int64, B::Float64, R::Vector{Float64}, J_prime::Vector{Int64}, M::Int64) 

    best_MMTs = []
    best_Budget = 0.0
    best_Q_MMTs = []
    best_Q_tot = 0
    best_VC = nothing
    
    for VC in 1:m
        MMTs, Budget, Q_MMTs, Q_tot = greedy_MMTs(VC, n, m, D, A, Q, C, q, f, p, B, J_prime, M)

        if Q_tot > best_Q_tot
            best_MMTs = MMTs
            best_Budget = Budget
            best_Q_MMTs = Q_MMTs
            best_Q_tot = Q_tot
            best_VC = VC
        end
    end

    best_Q_tot += sum(A[best_VC, j]*q[j] for j in 1:n)

    return best_VC, best_MMTs, best_Budget, best_Q_MMTs, best_Q_tot
end



function greedy_init(
    model::Model, MMTs, VC::Int64, M::Int64, n::Int64, m::Int64, A
)
    # Indices
    I = 1:m  # VC locations
    J = 1:n  # Localities

    # Extract variables from the model
    z = model[:z]
    delta = model[:delta]
    y = model[:y]
    u = model[:u]
    beta = model[:beta]
    v = model[:v]

    # Initialize VC locations
    for i in I
        set_start_value(y[i], 0)
    end
    set_start_value(y[VC], 1)

    # Initialize MMTs
    num_MMTs = length(MMTs)
    for idx in 1:M
        set_start_value(delta[idx], idx <= num_MMTs ? 1 : 0)
    end

    # Initialize routes
    for i in 1:(n + m), j in 1:(n + m), k in 1:M
        set_start_value(z[i, j, k], 0)
    end

    # Initialize MMT assignments
    for j in J, k in 1:M
        set_start_value(u[j, k], 0)
    end

    # Initialize beta for VC locations
    for i in I
        set_start_value(beta[i + n], i == VC ? 1 : 0)
    end

    # Initialize beta for localities
    for j in J
        set_start_value(beta[j], 2)
    end

    # Initialize coverage
    for j in J
        set_start_value(v[j], A[VC, j])
    end

    # Set up MMT routes
    for (idx, MMT) in enumerate(MMTs)
        Q_tot = 1
        for (i, j) in MMT[2:end-1]
            Q_tot += 1
            set_start_value(z[i, j, idx], 1)
            set_start_value(u[j, idx], 1)
            set_start_value(beta[j], Q_tot)
        end
        (iend, jend) = MMT[end]
        set_start_value(z[iend, jend, idx], 1)
    end
end



