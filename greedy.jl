function greedy_MMTs(VC::Int64, n::Int64, m::Int64, D::Array{Float64, 2}, A::Array{Int64, 2}, Q::Int64, C::Vector{Int64}, q::Vector{Int64}, f::Vector{Int64}, p::Int64, B::Float64, J_prime::Vector{Int64}, M::Int64)

    MMTs = [[(VC+n, VC+n)]]
    Q_MMTs = [0]

    Budget = p + f[VC]
    Q_tot = 0
    
    localities_to_visit = [j for j in 1:n if A[VC, j] == 0]
    
    while true

        # println("Localities to visit: ", localities_to_visit)
        
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

            elseif (q[j] <= Q) && (Budget + D[VC+n, j] + D[j, VC+n] + p <= B) && (Q_tot + q[j] <= C[VC]) && (length(MMTs) <= M)
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
    
    for VC in 1:m
        MMTs, Budget, Q_MMTs, Q_tot = greedy_MMTs(VC, n, m, D, A, Q, C, q, f, p, B, J_prime, M)

        if Q_tot > best_Q_tot
            best_MMTs = MMTs
            best_Budget = Budget
            best_Q_MMTs = Q_MMTs
            best_Q_tot = Q_tot
        end
    end


    return best_MMTs, best_Budget, best_Q_MMTs, best_Q_tot
end