###############################################################################
# Generate "cross moves" for Tabu Search:
# A cross move involves taking a node (locality) from one route (from_route) 
# and inserting it into another route (to_route) at a specific position.
###############################################################################
function generate_cross_moves(MMTs)
    moves = []

    for from_idx in eachindex(MMTs)
        from_route = MMTs[from_idx]
        from_node_seq = extract_node_seq(from_route)

        # from_node_seq = [VC+n, j1, j2, ..., jK, VC+n]
        # Localities are at positions from 2 to length(node_seq)-1.
        for i in 2:(length(from_node_seq)-1)
            from_loc = from_node_seq[i]

            # Attempt to insert from_loc into all other routes
            for to_idx in eachindex(MMTs)
                # Skip if it's the same route
                if to_idx == from_idx
                    continue
                end

                to_route = MMTs[to_idx]
                to_node_seq = extract_node_seq(to_route)

                # Possible insertion positions: 
                # between 2 and length(to_node_seq)-1 (locality positions)
                for to_insert in 2:(length(to_node_seq)-1)
                    push!(moves, (from_idx, from_loc, to_idx, to_insert))
                end
            end
        end
    end

    return moves
end


###############################################################################
# Apply a cross move to the solution:
# Removes the specified locality from one route and inserts it into another.
###############################################################################
function apply_cross_move(MMTs, move)
    # Move format: (from_idx, from_loc, to_idx, to_insert)
    (from_idx, from_loc, to_idx, to_insert) = move
    new_MMTs = deepcopy(MMTs)

    #####################
    # Update from_route #
    #####################
    from_route = new_MMTs[from_idx]
    from_node_seq = extract_node_seq(from_route)

    # Remove locality from 'from' route
    deleteat!(from_node_seq, findfirst(==(from_loc), from_node_seq))
    new_MMTs[from_idx] = generate_route(from_node_seq)

    ###################
    # Update to_route #
    ###################
    to_route = new_MMTs[to_idx]
    to_node_seq = extract_node_seq(to_route)

    # Insert locality into 'to' route
    insert!(to_node_seq, to_insert, from_loc)
    new_MMTs[to_idx] = generate_route(to_node_seq)

    return new_MMTs
end


###############################################################################
# Generate 2-opt moves:
# A 2-opt move selects a route and reverses a contiguous segment of it to 
# potentially reduce travel distance.
###############################################################################
function generate_2opt_moves(MMTs)
    moves = []

    # Each move is of the form (route_idx, start_pos, end_pos) 
    # indicating which segment of the route to reverse.
    for route_idx in eachindex(MMTs)
        node_seq = extract_node_seq(MMTs[route_idx])

        # Possible segments to reverse are between node_seq[2] and node_seq[end-1]
        for i in 2:(length(node_seq)-2)
            for j in (i+1):(length(node_seq)-1)
                push!(moves, (route_idx, i, j))
            end
        end
    end

    return moves
end


###############################################################################
# Apply a 2-opt move:
# Reverses the segment of the route defined by (i:j).
###############################################################################
function apply_2opt_move(MMTs, move)
    # Move format: (route_idx, i, j)
    (route_idx, i, j) = move
    new_MMTs = deepcopy(MMTs)

    route = new_MMTs[route_idx]
    node_seq = extract_node_seq(route)

    # Reverse the [i:j] segment
    new_node_seq = deepcopy(node_seq)
    new_node_seq[i:j] = reverse(node_seq[i:j])

    # Rebuild the route from the modified node sequence
    new_route = []
    for idx in 1:(length(new_node_seq)-1)
        push!(new_route, (new_node_seq[idx], new_node_seq[idx+1]))
    end

    new_MMTs[route_idx] = new_route
    return new_MMTs
end
