function preprocess(n, m, x_VC, y_VC, x_loc, y_loc, R)

    distance_matrix = zeros(Float64, m, n)

    for i in 1:m 
        for j in 1:n
            distance_matrix[i, j] = sqrt((x_VC[i] - x_loc[j])^2 + (y_VC[i] - y_loc[j])^2)
        end
    end

    A = zeros(Int, m, n)


    for i in 1:m
        for j in 1:n
            if distance_matrix[i, j] <= R[i]
                A[i, j] = 1
            end
        end
    end


    D = zeros(Float64, n+m, n+m)
    #The n first rows/cols of D relate to the localities, the m last to the VC
    #In other words, calling D[i,j] for i in I, j in J corresponds to call D[i+n,j]
    for i in 1:m+n
        for j in 1:m+n
            if i <= n && j <= n
                # Both indices in `loc`
                D[i, j] = sqrt((x_loc[i] - x_loc[j])^2 + (y_loc[i] - y_loc[j])^2)
            elseif i > n && j > n
                # Both indices in `VC`
                D[i, j] = sqrt((x_VC[i - n] - x_VC[j - n])^2 + (y_VC[i - n] - y_VC[j - n])^2)
            elseif i > n && j <= n
                # i in `VC`, j in `loc`
                D[i, j] = sqrt((x_VC[i - n] - x_loc[j])^2 + (y_VC[i - n] - y_loc[j])^2)
            else
                # i in `loc`, j in `VC`
                D[i, j] = sqrt((x_loc[i] - x_VC[j - n])^2 + (y_loc[i] - y_VC[j - n])^2)
            end
        end
    end

    return A, D
end

# Input: Path of the instance 
function parse_instance(path::String)
    # Open the file in reading mode
    file = open(path, "r") 

    # Reading the first-variables to get the dimension of the problem
    # Number of localities
    n = parse(Int, readline(file))

    # Number of VC
    m = parse(Int, readline(file))

    # Trick: Empty line --> Going to the next line
    readline(file)
    
    # Declaring and pre-allocating data variables
    # coordinates for the VC
    x_VC, y_VC = Vector{Float64}(undef,m), Vector{Float64}(undef,m)  
    # coordinates for the localities
    x_loc, y_loc = Vector{Float64}(undef,n), Vector{Float64}(undef,n) 
    # vehicles capacity
    Q = nothing
    # depot capacities
    C = Vector{Int}(undef, m)
    # Population for each locality 
    q = Vector{Int}(undef, n)
    # Opening costs for the VC 
    f = Vector{Float64}(undef, m) 
    # Opening cost of a route / cost of a vechicle
    p = nothing
    # Type of cost: 0 or 1 (0 means that the costs are integer - 1 that costs are real)
    tc = nothing 
    # Budget
    B = nothing
    # Radius for each VC
    R = Vector{Float64}(undef, m)
    # Localities where patients should be immunized in priority
    localities_with_high_priorities = Vector{Int}()
    # Maximum number of vehicles
    M = nothing

    # coordinates for the VC (x and y)
    for i in 1:m 
        str_mixed_xy = split(readline(file))
        x_VC[i], y_VC[i] = parse(Float64,str_mixed_xy[1]), parse(Float64,str_mixed_xy[2])
    end
    readline(file)

    # coordinates for the localities
    for i in 1:n
        str_mixed_xy = split(readline(file))
        x_loc[i], y_loc[i] = parse(Float64,str_mixed_xy[1]), parse(Float64,str_mixed_xy[2])
    end
    readline(file)

    # vehicle capacity
    Q = parse(Int,readline(file))
    readline(file)

    # depot capacities
    for i in 1:m
        C[i] = parse(Int, readline(file))
    end
    readline(file)

    # number of people in localities
    for i in 1:n
        q[i] = parse(Int, readline(file)) 
    end
    readline(file)

    # opening costs for the VC
    for i in 1:m
        f[i] = parse(Float64, readline(file))
    end
    readline(file)

    # opening costs of a route (cost of a vehicle)
    p = parse(Float64, readline(file))
    readline(file)

    # 0 or 1 (0 means that the costs are integer - 1 that costs are real)
    tc = parse(Int, readline(file))
    readline(file); readline(file) # two line breaks

    # Budget
    B = parse(Float64, readline(file))
    readline(file)

    # Radius for each VC
    for k=1:m
        R[k] = parse(Float64,readline(file))    
    end 
    readline(file)

    # Localities where patients should be immunized in priority
    there_is_data = true
    while there_is_data 
        line = readline(file)
        there_is_data = ! isempty(line)
        if there_is_data
            push!(localities_with_high_priorities,parse(Int,line))
        end
    end

    # Maximum number of vehicles
    M = parse(Int,readline(file))

    close(file)

    if tc == 1
        return n, m, x_VC, y_VC, x_loc, y_loc, Q, C, q, f, p, tc, B, R, localities_with_high_priorities, M
    else # tc == 0
        return  n, m, x_VC, y_VC, x_loc, y_loc, Q, C, q, convert(Vector{Int},f), convert(Int,p), tc, B, R, localities_with_high_priorities, M # Converting data types of costs from float to integers
    end
end



function parse_to_struct(path::String)

    n, m, x_VC, y_VC, x_loc, y_loc, Q, C, q, f, p, tc, B, R, localities_with_high_priorities, M = parse_instance(path)
    A, D = preprocess(n, m, x_VC, y_VC, x_loc, y_loc, R)


    return OptVaxData(n, m, D, A, Q, C, q, f, p, B, localities_with_high_priorities, M)

end


#=
The structure of the instance files is as follows:

number of customers
number of available VC

coordinates for the VC (x and y)

coordinates for the customers

vehicle capacity

depot capacities (for Tuzun instances, each one is equal to the total demand as there is no capacity on the VC)

customers demands

opening costs for the VC

opening cost of a route (cost of a vehicle)

0 or 1 (0 means that the costs are integer - 1 that costs are real)

Budget

Radius for each VC

Localities where patients should be immunized in priority

Maximum number of vehicles

\\end of file
- - - - - - - - - - 
To calculate the matrix distance (or the cost to link any 2 points A and B in the graph), we use the mathematical formula:

sqrt( (xA-xB)² + (yA-yB)² )
=#