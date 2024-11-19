using Plots

function plot_solution(model, x_VC, y_VC, x_loc, y_loc, R, n, m, M, localities_with_high_priorities)
    # Retrieve the optimized values of y, z, and delta from the model
    y_values = value.(model[:y])
    z_values = value.(model[:z])  # z_values is a 3D array with dimensions (i, j, k)
    delta_values = value.(model[:delta])

    # Create a new plot
    color = :blue
    plt = plot()
        # Plot loc points without adding them to the legend (but with single labels for groups)
    first_black = true
    first_green = true
    for j in 1:n
        if j in localities_with_high_priorities
            # Plot green square points and add label only for the first one
            if first_green
                scatter!(plt, [x_loc[j]], [y_loc[j]], color=:green, marker=:square, aspect_ratio=:equal, label="High priority")
                first_green = false
            else
                scatter!(plt, [x_loc[j]], [y_loc[j]], color=:green, marker=:square, aspect_ratio=:equal, label="")
            end
        else
            # Plot black circle points and add label only for the first one
            if first_black
                scatter!(plt, [x_loc[j]], [y_loc[j]], color=:black, marker=:circle, aspect_ratio=:equal, label="Low priority")
                first_black = false
            else
                scatter!(plt, [x_loc[j]], [y_loc[j]], color=:black, marker=:circle, aspect_ratio=:equal, label="")
            end
        end
    end

    # For each VC point, check if y[i] is 1; if so, plot in red and add a shaded circle
    for i in 1:m 
        if y_values[i] == 1
            # Plot the VC point in red
            scatter!(plt, [x_VC[i]], [y_VC[i]], color=:red, marker=:star5, label="VC")
            
            # Draw a shaded circle of radius R[i] around the VC location
            θ = range(0, 2π, length=100)
            circle_x = x_VC[i] .+ R[i] .* cos.(θ)
            circle_y = y_VC[i] .+ R[i] .* sin.(θ)
            plot!(plt, circle_x, circle_y, lw=1.5, linecolor=:red, label="", fillalpha=0.1, fillcolor=:red)
        end
    end

    # Create a color gradient (you can use other gradient types too)
    color_map = [:blue,:teal,:green,:darkred,:orange]

    # For all z[i,j,k], draw arrows only if z[i,j,k] == 1 for any k
    for i in 1:(m+n), j in 1:(m+n), k in 1:M
        if abs(z_values[i, j, k] - 1.0) < 1e-1
            # Scale the color based on the value of k, you can normalize it if needed
            color = color_map[k]  # This gives a color based on the value of k

            # Create the arrows with the color mapped to k
            if i <= n && j <= n
                # Arrow from loc[i] to loc[j]
                dx = x_loc[j] - x_loc[i]
                dy = y_loc[j] - y_loc[i]
                quiver!(plt, [x_loc[i]], [y_loc[i]], quiver=([dx], [dy]), color=color, linewidth=1, arrowhead=2)
            elseif i > n && j > n
                # Arrow from VC[i] to VC[j]
                dx = x_VC[j - n] - x_VC[i - n]
                dy = y_VC[j - n] - y_VC[i - n]
                quiver!(plt, [x_VC[i - n]], [y_VC[i - n]], quiver=([dx],[dy]), color=color, linewidth=1, arrowhead=2)
            elseif i <= n && j > n
                # Arrow from loc[i] to VC[j]
                dx = x_VC[j - n] - x_loc[i]
                dy = y_VC[j - n] - y_loc[i]
                quiver!(plt, [x_loc[i]], [y_loc[i]], quiver=([dx], [dy]), color=color, linewidth=1, arrowhead=2)
            elseif i > n && j <= n
                # Arrow from VC[i] to loc[j]
                dx = x_loc[j] - x_VC[i - n]
                dy = y_loc[j] - y_VC[i - n]
                quiver!(plt, [x_VC[i - n]], [y_VC[i - n]], quiver=([dx], [dy]), color=color, linewidth=1, arrowhead=2)
            end
        end
    end

    # Display the final plot
    savefig(plt, "test.pdf")
    display(plt)
end


function plot_solution_greedy(MMTs, x_VC, y_VC, x_loc, y_loc, R, n, m, M, localities_with_high_priorities)

    # Create a new plot
    color = :blue
    plt = plot()
        # Plot loc points without adding them to the legend (but with single labels for groups)
    first_black = true
    first_green = true

    for j in 1:n
        if j in localities_with_high_priorities
            # Plot green square points and add label only for the first one
            if first_green
                scatter!(plt, [x_loc[j]], [y_loc[j]], color=:green, marker=:square, aspect_ratio=:equal, label="High priority")
                first_green = false
            else
                scatter!(plt, [x_loc[j]], [y_loc[j]], color=:green, marker=:square, aspect_ratio=:equal, label="")
            end
        else
            # Plot black circle points and add label only for the first one
            if first_black
                scatter!(plt, [x_loc[j]], [y_loc[j]], color=:black, marker=:circle, aspect_ratio=:equal, label="Low priority")
                first_black = false
            else
                scatter!(plt, [x_loc[j]], [y_loc[j]], color=:black, marker=:circle, aspect_ratio=:equal, label="")
            end
        end
    end

    VC = MMTs[1][1][1]
    
    
    # Plot the VC point in red
    scatter!(plt, [x_VC[VC-n]], [y_VC[VC-n]], color=:red, marker=:star5, label="VC")
    
    # Draw a shaded circle of radius R[i] around the VC location
    θ = range(0, 2π, length=100)
    circle_x = x_VC[VC-n] .+ R[VC-n] .* cos.(θ)
    circle_y = y_VC[VC-n] .+ R[VC-n] .* sin.(θ)
    plot!(plt, circle_x, circle_y, lw=1.5, linecolor=:red, label="", fillalpha=0.1, fillcolor=:red, axis = nothing)
    
    
    # Create a color gradient (you can use other gradient types too)
    color_map = [:blue,:teal,:green,:darkred,:orange, :purple, :pink, :yellow, :brown, :cyan]
    
    # For all z[i,j,k], draw arrows only if z[i,j,k] == 1 for any k
        
    
    for (k, MMT) in enumerate(MMTs)
        for (i, j) in MMT
            color = color_map[k]  
            # Create the arrows with the color mapped to k
            if i <= n && j <= n
                # Arrow from loc[i] to loc[j]
                dx = x_loc[j] - x_loc[i]
                dy = y_loc[j] - y_loc[i]
                quiver!(plt, [x_loc[i]], [y_loc[i]], quiver=([dx], [dy]), color=color, linewidth=1, arrowhead=2)
            elseif i > n && j > n
                # Arrow from VC[i] to VC[j]
                dx = x_VC[j - n] - x_VC[i - n]
                dy = y_VC[j - n] - y_VC[i - n]
                quiver!(plt, [x_VC[i - n]], [y_VC[i - n]], quiver=([dx],[dy]), color=color, linewidth=1, arrowhead=2)
            elseif i <= n && j > n
                # Arrow from loc[i] to VC[j]
                dx = x_VC[j - n] - x_loc[i]
                dy = y_VC[j - n] - y_loc[i]
                quiver!(plt, [x_loc[i]], [y_loc[i]], quiver=([dx], [dy]), color=color, linewidth=1, arrowhead=2)
            elseif i > n && j <= n
                # Arrow from VC[i] to loc[j]
                dx = x_loc[j] - x_VC[i - n]
                dy = y_loc[j] - y_VC[i - n]
                quiver!(plt, [x_VC[i - n]], [y_VC[i - n]], quiver=([dx], [dy]), color=color, linewidth=1, arrowhead=2)
            end
        end
    end

    # Display the final plot
    # savefig(plt, "test.pdf")
    display(plt)
end