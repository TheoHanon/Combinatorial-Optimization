include("model.jl")
include("Parser.jl")
include("greedy.jl")
include("utils.jl")


using JSON
using CPUTime
using JuMP
using HiGHS


function baseline_solve(n, m, D, A, Q, C, q, f, p, B, R, localities_with_high_priorities, M)

    model = solve_OptVax2(n, m, D, A, Q, C, q, f, p, B, R, localities_with_high_priorities, M)

    optimize!(model)
    termination_status = JuMP.termination_status(model)

    if !(termination_status == MOI.FEASIBLE_POINT)
        obj_value = 0.0
     else
         obj_value = objective_value(model)
     end
    
    return obj_value

end

function solve(n, m, D, A, Q, C, q, f, p, B, R, localities_with_high_priorities, M, model_type = "1", warm_start = false)

    if model_type == "1"
        model = solve_OptVax1(n, m, D, A, Q, C, q, f, p, B, R, localities_with_high_priorities, M)
    elseif model_type == "2"
        model = solve_OptVax2(n, m, D, A, Q, C, q, f, p, B, R, localities_with_high_priorities, M)
    else
        model = solve_OptVax2LP(n, m, D, A, Q, C, q, f, p, B, R, localities_with_high_priorities, M)
    end

    if warm_start 
        VC, MMTs, Budget, Q_MMTs, Q_tot = greedy_OptVax(n, m, D, A, Q, C, q, f, p, B, R, localities_with_high_priorities, M)
        new_MMTs = local_search_2Opt(MMTs, D)
        greedy_init(model, new_MMTs, VC, M, n, m, A)
    end

    optimize!(model)

    termination_status = JuMP.termination_status(model)

    if !(termination_status == MOI.FEASIBLE_POINT)
       obj_value = 0.0
    else
        obj_value = objective_value(model)
    end

    return obj_value
end

function solve_greedy(n, m, D, A, Q, C, q, f, p, B, R, localities_with_high_priorities, M)

    VC, MMTs, Budget, Q_MMTs, Q_tot = greedy_OptVax(n, m, D, A, Q, C, q, f, p, B, R, localities_with_high_priorities, M)
    new_MMTs = local_search_2Opt(MMTs, D)

    return Q_tot

end

source= "Combinatorial-Optimization/Instances"
files = readdir(source)

for i in reverse(1:length(files))
    file = files[i]
    if occursin("100", file) || occursin("200", file)
            deleteat!(files, i)
    end
end


z_star = zeros(length(files))
z_tilde = zeros(6, length(files))
time = zeros(6, length(files))


for (k,file) in enumerate(files) 

    println("Solving instance $file")
    n, m, x_VC, y_VC, x_loc, y_loc, Q, C, q, f, p, tc, B, R, localities_with_high_priorities, M = parse_instance(joinpath(source, file))
    A, D = preprocess(n, m, x_VC, y_VC, x_loc, y_loc, R)

    # Get results for model 1 with/without warm start

    t0 = CPUtime_us()
    z_tilde[1, k] = solve(n, m, D, A, Q, C, q, f, p, B, R, localities_with_high_priorities, M, "1", false)
    t1 = CPUtime_us()
    time[1, k] = (t1 - t0) / 1e6

    t0 = CPUtime_us()
    z_tilde[5, k] = solve(n, m, D, A, Q, C, q, f, p, B, R, localities_with_high_priorities, M, "1", true)
    t1 = CPUtime_us()
    time[5, k] = (t1 - t0) / 1e6

    # Get results for model 2 with/without warm start
    t0 = CPUtime_us()
    z_star[k] = baseline_solve(n, m, D, A, Q, C, q, f, p, B, R, localities_with_high_priorities, M)
    z_tilde[2, k] = z_star[k]
    t1 = CPUtime_us()
    time[2, k] = (t1 - t0) / 1e6

    t0 = CPUtime_us()
    z_tilde[6, k] = solve(n, m, D, A, Q, C, q, f, p, B, R, localities_with_high_priorities, M, "2", true)
    t1 = CPUtime_us()
    time[6, k] = (t1 - t0) / 1e6


    # Get results for model 2LP

    t0 = CPUtime_us()
    z_tilde[3, k] = solve(n, m, D, A, Q, C, q, f, p, B, R, localities_with_high_priorities, M, "LP", false)
    t1 = CPUtime_us()
    time[3, k] = (t1 - t0) / 1e6

    # Get results of Greedy + warm start

    t0 = CPUtime_us()
    z_tilde[4, k] = solve_greedy(n, m, D, A, Q, C, q, f, p, B, R, localities_with_high_priorities, M)
    t1 = CPUtime_us()
    time[4, k] = (t1 - t0) / 1e6

end

results = Dict(
    "z_star" => z_star,
    "z_tilde" => z_tilde,
    "time" => time
)

open("Combinatorial-Optimization/results.json", "w") do io
    JSON.print(io, results)
end
