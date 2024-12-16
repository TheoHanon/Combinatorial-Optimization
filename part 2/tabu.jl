using JuMP
include("model.jl")
include("tabu_moves.jl")

# --------------------------------------------------------------------------
# Evaluate Q: The total vaccination quantity delivered by the MMTs plus that 
# obtained from the VC assignment.
function evaluate_Q(MMTs, data)
    VC = MMTs[1][1][1] - data.n
    Q_tot = 0.0

    # Sum quantities delivered on all routes
    for route in MMTs
        node_seq = extract_node_seq(route)
        for i in 2:(length(node_seq)-1)
            Q_tot += data.q[node_seq[i]]
        end
    end

    # Add quantity from assigned localities (A matrix)
    Q_tot += sum(data.A[VC, j] * data.q[j] for j in 1:data.n)
    return Q_tot
end

# --------------------------------------------------------------------------
# Evaluate solution: Returns total_cost - Q_tot for given MMTs configuration.
function evaluate_solution(MMTs, data)
    VC = MMTs[1][1][1] - data.n
    total_cost = 0.0
    Q_tot = 0.0

    # Sum travel costs and quantities
    for route in MMTs
        for (i, j) in route
            total_cost += data.D[i, j]
            if j <= data.n
                Q_tot += data.q[j]
            end
        end
    end

    Q_tot += sum(data.A[VC, j] * data.q[j] for j in 1:data.n)
    return total_cost - Q_tot
end

# --------------------------------------------------------------------------
# Add a move to the tabu list and maintain tabu tenure size.
function add_to_tabu_list(tabu_list, move, tabu_tenure)
    push!(tabu_list, (move, tabu_tenure))
    if length(tabu_list) > tabu_tenure
        popfirst!(tabu_list)
    end
end

# --------------------------------------------------------------------------
# Check if a move is in the tabu list.
function is_tabu(tabu_list, move)
    for (tabu_move, tenure) in tabu_list
        if tabu_move == move
            return true
        end
    end
    return false
end

# --------------------------------------------------------------------------
# Extract node sequence from a route defined as pairs (i,j).
function extract_node_seq(route)
    node_seq = [route[1][1]]
    for (i, j) in route
        push!(node_seq, j)
    end
    return node_seq
end

# --------------------------------------------------------------------------
# Generate a route from a node sequence by pairing consecutive nodes.
function generate_route(node_seq)
    route = []
    for i in 1:(length(node_seq)-1)
        push!(route, (node_seq[i], node_seq[i+1]))
    end
    return route
end

# --------------------------------------------------------------------------
# Check feasibility of MMTs solution given data constraints.
function is_feasible(MMTs, data)
    VC = MMTs[1][1][1] - data.n

    B = data.p * length(MMTs) + data.f[VC]
    Q_tot = 0.0
    Q_MMTs = zeros(length(MMTs))

    # Accumulate budget and quantities
    for (idx, route) in enumerate(MMTs)
        for (i, j) in route
            B += data.D[i, j]
            if j < data.n
                Q_tot += data.q[j]
                Q_MMTs[idx] += data.q[j]
            end
        end
    end

    # Check overall constraints
    return B <= data.B && Q_tot <= data.C[VC] && all(Q_MMTs .<= data.Q)
end

# --------------------------------------------------------------------------
# Tabu Search: Attempts to improve solution using tabu moves.
function tabu_search(data::OptVaxData;
                     alpha::Float64 = 0.1,
                     max_iter::Int64 = 1000,
                     no_improve_limit::Int64 = 100)

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

    # Initialize tabu structure
    tabu_list = Vector{Tuple{Any,Int}}()
    tabu_tenure = round(Int, alpha * n)

    # Get initial solution from greedy approach
    best_VC, best_MMTs = greedy_OptVax(data) 
    best_cost = evaluate_solution(best_MMTs, data)
    best_Q = evaluate_Q(best_MMTs, data)

    current_MMTs = deepcopy(best_MMTs)
    current_cost = best_cost
    current_VC = best_VC
    current_Q = best_Q

    best_solution = (best_VC, best_MMTs, best_cost, best_Q)
    no_improve_count = 0

    # Main search loop
    for iteration in 1:max_iter
        two_opt_moves = generate_2opt_moves(current_MMTs)
        cross_moves = generate_cross_moves(current_MMTs)

        best_move = nothing
        best_move_cost = Inf
        best_move_Q = Inf
        best_move_solution = nothing

        # Update tabu tenure countdown
        for i in reverse(eachindex(tabu_list))
            (tabu_move, tenure) = tabu_list[i]
            if tenure - 1 <= 0
                deleteat!(tabu_list, i)
            else
                tabu_list[i] = (tabu_move, tenure - 1)
            end
        end

        # Explore moves (two-opt and cross moves)
        for move in [two_opt_moves; cross_moves]
            # Check if move is tabu
            if is_tabu(tabu_list, move)
                # Aspiration criteria: allow if it improves best known cost
                if length(move) == 3
                    new_MMTs = apply_2opt_move(current_MMTs, move)
                else
                    new_MMTs = apply_cross_move(current_MMTs, move)
                    if !is_feasible(new_MMTs, data)
                        continue
                    end
                end

                new_cost = evaluate_solution(new_MMTs, data)
                new_Q = evaluate_Q(new_MMTs, data)
                if new_cost < best_cost
                    if new_cost < best_move_cost
                        best_move = move
                        best_move_cost = new_cost
                        best_move_Q = new_Q
                        best_move_solution = new_MMTs
                    end
                end

            else
                # Move is not tabu
                if length(move) == 3
                    new_MMTs = apply_2opt_move(current_MMTs, move)
                else
                    new_MMTs = apply_cross_move(current_MMTs, move)
                    if !is_feasible(new_MMTs, data)
                        continue
                    end
                end

                new_cost = evaluate_solution(new_MMTs, data)
                new_Q = evaluate_Q(new_MMTs, data)
                if new_cost < best_move_cost
                    best_move = move
                    best_move_cost = new_cost
                    best_move_Q = new_Q
                    best_move_solution = new_MMTs
                end
            end
        end

        if best_move === nothing
            # No moves found or none are better
            no_improve_count += 1
            if no_improve_count > no_improve_limit
                break
            end
        else
            # Apply best move
            current_MMTs = best_move_solution
            current_cost = best_move_cost
            current_Q = best_move_Q

            # Add move to tabu list
            add_to_tabu_list(tabu_list, best_move, tabu_tenure)

            # Update global best if improved
            if current_cost < best_cost
                best_cost = current_cost
                best_Q = current_Q
                best_MMTs = current_MMTs
                best_solution = (best_VC, best_MMTs, best_cost, best_Q)
                no_improve_count = 0
            else
                no_improve_count += 1
                if no_improve_count > no_improve_limit
                    break
                end
            end
        end

        # Apply Greedy starting from current best solution
        new_MMTs, _, _, _ = greedy_MMTs(best_VC, data, best_MMTs)
        new_Q = evaluate_Q(new_MMTs, data)
        new_cost = evaluate_solution(new_MMTs, data)

        if new_cost < best_cost
            best_cost = new_cost
            best_Q = new_Q
            best_MMTs = new_MMTs
            best_solution = (best_VC, best_MMTs, best_cost, best_Q)
            no_improve_count = 0
        end
    end

    print("Is feasible: ", is_feasible(best_MMTs, data))
    return best_solution
end

# --------------------------------------------------------------------------
# Greedy MMTs construction given a VC and possibly an existing MMTs structure.
function greedy_MMTs(VC::Int64, data::OptVaxData, MMTs = nothing)
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

    if isnothing(MMTs)
        # Start with a single empty MMT containing just VC start/end
        MMTs = [[(VC+n, VC+n)]]
        Q_MMTs = [0.0]
        Budget = p + f[VC]
        Q_tot = 0
        localities_to_visit = [j for j in 1:n if A[VC, j] == 0]
    else
        # Remove empty MMTs of form [(x,x)]
        del_idx = []
        for (idx, MMT) in enumerate(MMTs)
            if MMT[1][1] == MMT[1][2]
                push!(del_idx, idx)
            end
        end
        for idx in reverse(del_idx)
            deleteat!(MMTs, idx)
        end

        Q_MMTs = [sum(data.q[j] for (i,j) in route[1:end-1]) for route in MMTs]
        Q_tot = sum(Q_MMTs) + sum(data.A[VC, j]*data.q[j] for j in 1:data.n)

        MMTs = deepcopy(MMTs)
        for MMT in MMTs
            deleteat!(MMT, lastindex(MMT)) # Remove return to VC for the construction phase
        end

        Budget = length(MMTs)*p + f[VC] + sum(D[i,j] for route in MMTs for (i,j) in route)
        localities_to_visit = [j for j in 1:n if A[VC,j] == 0]

        # Remove already visited localities
        for route in MMTs
            node_seq = extract_node_seq(route)
            for i in 2:(length(node_seq)-1)
                deleteat!(localities_to_visit, findfirst(==(node_seq[i]), localities_to_visit))
            end
        end
    end

    # Assign localities to routes based on a greedy heuristic
    while true
        assigned = false

        utilities = [(q[j]/D[MMT[end][2], j], idx, j) 
                     for (idx, MMT) in enumerate(MMTs)
                     for j in localities_to_visit if !(j in J_prime)]

        utilities_prio = [(q[j]/D[MMT[end][2], j], idx, j) 
                          for (idx, MMT) in enumerate(MMTs) 
                          for j in localities_to_visit if j in J_prime]

        # Put high priority localities first
        sorted_utilities = vcat(sort(utilities_prio, by = x -> -x[1]),
                                sort(utilities, by = x -> -x[1]))

        for (utility, idx, j) in sorted_utilities
            route_end = MMTs[idx][end][2]

            # Try to add locality j to existing route if feasible
            if (Q_MMTs[idx] + q[j] <= Q) &&
               (Budget + D[route_end, j] + D[j, VC+n] <= B) &&
               (Q_tot + q[j] <= C[VC]) &&
               (route_end != j)

                push!(MMTs[idx], (route_end, j))
                Q_MMTs[idx] += q[j]
                Budget += D[route_end, j]
                Q_tot += q[j]
                deleteat!(localities_to_visit, findfirst(==(j), localities_to_visit))
                assigned = true
                break

            # If not possible in existing routes, try creating a new MMT if feasible
            elseif (q[j] <= Q) &&
                  (Budget + D[VC+n, j] + D[j, VC+n] + p <= B) &&
                  (Q_tot + q[j] <= C[VC]) &&
                  (length(MMTs) < M)

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

    # Close routes by returning to VC
    for MMT in MMTs
        push!(MMT, (MMT[end][2], VC+n))
        if MMT[1] == (VC+n, VC+n)
            deleteat!(MMT, 1)
        elseif MMT[end] == (VC+n, VC+n)
            deleteat!(MMT, lastindex(MMT))
        end
        Budget += D[MMT[end][2], VC+n]
    end

    return MMTs, Budget, Q_MMTs, Q_tot
end

# --------------------------------------------------------------------------
# Greedy OptVax solution: finds best VC by applying greedy_MMTs and choosing 
# the VC with max Q_tot.
function greedy_OptVax(data::OptVaxData)
    best_MMTs = []
    best_Budget = 0.0
    best_Q_MMTs = []
    best_Q_tot = 0
    best_VC = nothing

    for VC in 1:data.m
        MMTs, Budget, Q_MMTs, Q_tot = greedy_MMTs(VC, data)
        if Q_tot > best_Q_tot
            best_MMTs = MMTs
            best_Budget = Budget
            best_Q_MMTs = Q_MMTs
            best_Q_tot = Q_tot
            best_VC = VC
        end
    end

    best_Q_tot += sum(data.A[best_VC, j] * data.q[j] for j in 1:data.n)
    return best_VC, best_MMTs
end