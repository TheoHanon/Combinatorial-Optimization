function generate_cross_moves(MMTs)

    moves = []

    for from_idx in eachindex(MMTs)

        from_route = MMTs[from_idx]
        from_node_seq = extract_node_seq(from_route)
        

        # Route looks like: (VC+n, VC+n), (VC+n, j1), (j1, j2), ..., (jK, VC+n)
        # node_seq = [VC+n, j1, j2, ..., jK, VC+n]
        # Localities are j1...jK
        # i from 2 to length(node_seq)-1 to skip start/end

        for i in 2:length(from_node_seq)-1
            from_loc = from_node_seq[i]

            for to_idx in eachindex(MMTs)
                if to_idx == from_idx
                    continue
                end

                to_route = MMTs[to_idx]
                to_node_seq = extract_node_seq(to_route)

                for to_insert in 2:length(to_node_seq)-1
                    # We will consider a move: remove `locality` from from_route, insert into to_route at to_insert
                    push!(moves, (from_idx, from_loc, to_idx, to_insert))
                end
            end
        end
    end

    return moves

end


function apply_cross_move(MMTs, move)

    new_MMTs = deepcopy(MMTs)

    # Move format: (from_idx, from_loc, to_idx, to_insert)
    (from_idx, from_loc, to_idx, to_insert) = move

    #####################
    # Update from_route #
    #####################

    from_route = new_MMTs[from_idx]
    from_node_seq = extract_node_seq(from_route)
    deleteat!(from_node_seq, findfirst(==(from_loc), from_node_seq))
    new_MMTs[from_idx] = generate_route(from_node_seq)

    ###################
    # Update to_route #
    ###################

    to_route = new_MMTs[to_idx]
    to_node_seq = extract_node_seq(to_route)
    insert!(to_node_seq, to_insert, from_loc)
    new_MMTs[to_idx] = generate_route(to_node_seq)

    return new_MMTs

end


function generate_2opt_moves(MMTs) 

    moves = []

    # Each move: (route_idx, start_pos, end_pos)
    # where route_idx is the index of the MMT route in MMTs,
    # and start_pos, end_pos define the segment of the route to reverse.

    for route_idx in eachindex(MMTs)

        node_seq = extract_node_seq(MMTs[route_idx])
    
        for i in 2:length(node_seq)-2
            for j in i+1:length(node_seq)-1
                push!(moves, (route_idx, i, j))
            end
        end
    end

    return moves
end



function apply_2opt_move(MMTs, move)

    # move is (route_idx, i, j)

    (route_idx, i, j) = move
    new_MMTs = deepcopy(MMTs)
    route = new_MMTs[route_idx]
    node_seq = extract_node_seq(route)

    # Perform the 2-opt reversal on [i:j] segment of node_seq
    new_node_seq = deepcopy(node_seq)
    new_node_seq[i:j] = reverse(node_seq[i:j])

    # Rebuild route from new_node_seq
    # Each consecutive pair forms an edge
    new_route = []
    for idx in 1:length(new_node_seq)-1
        push!(new_route, (new_node_seq[idx], new_node_seq[idx+1]))
    end

    new_MMTs[route_idx] = new_route
    return new_MMTs
end

