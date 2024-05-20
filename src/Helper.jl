module Helper
    import IterTools as itr
    using ..AcceleratedDubins

    function read_tsp_file(filename::String)
        coordinates::Array{Tuple{Float64, Float64}, 1} = []

        open(filename) do file
            in_coordinates = false
            
            for line in eachline(file)
                if in_coordinates
                    if line == "EOF"
                        break
                    end

                    tokens = split(line)

                    x_coord = parse(Float64, tokens[2])
                    y_coord = parse(Float64, tokens[3])
                    
                    # Store the coordinates as a tuple
                    push!(coordinates, (x_coord, y_coord))
                else
                    # Check if we've reached the NODE_COORD_SECTION
                    if line == "NODE_COORD_SECTION"
                        in_coordinates = true
                    end
                end
            end
        end

        return coordinates
    end

    function find_lowest_edge(graph::Array{Float64, 6}, starting::Int64, ending::Int64, num_headings::Int64, num_speeds::Int64)
        best = Inf

        best_starting::Tuple{Int64, Int64, Int64} = (starting, -1, -1)
        best_ending::Tuple{Int64, Int64, Int64} = (ending, -1, -1)

        for (v_i, v_f) in itr.product(1:num_speeds, 1:num_speeds)
            for (h_i, h_f) in itr.product(1:num_headings, 1:num_headings)
                total_time = graph[starting, ending, v_i, v_f, h_i, h_f] + graph[ending, starting, v_f, v_i, h_f, h_i]
                if total_time < best
                    best = total_time
                    best_starting = (starting, v_i, h_i)
                    best_ending = (ending, v_f, h_f)
                end
            end
        end

        return best_starting, best_ending, best
    end

    function find_lowest_edge_between(graph::Array{Float64}, starting::Tuple{Int64, Int64, Int64}, ending::Tuple{Int64, Int64, Int64}, middle::Int64, num_headings::Int64, num_speeds::Int64)
        best_time = Inf
        best_middle::Tuple{Int64, Int64, Int64} = (middle, -1, -1)
        start_to_end = graph[starting[1], ending[1], starting[2], ending[2], starting[3], ending[3]]

        #Starting and ending configurations are already defined, need to define middle configurations
        for speed in 1:num_speeds
            for heading in 1:num_headings
                #dist(start,middle) + dist(middle, ending) - dist(start,ending)
                dist = graph[starting[1], middle, starting[2], speed, starting[3], heading] + graph[middle, ending[1], speed, ending[2], heading, ending[3]] - start_to_end
                if dist < best_time
                    best_time = dist
                    best_middle = (middle, speed, heading)
                end
            end
        end

        best_middle, best_time
    end

    function retrieve_path(locations::Vector{Tuple{Float64, Float64}}, configurations::Vector{Tuple{Int64, Int64, Int64}}, parameters, speeds::Vector{Float64}, headings::Vector{Float64}, radii::Vector{Float64})
        #(dubinspath, starting_speed, ending_speed)
        full_path::Vector{Tuple{AcceleratedDubins.DubinsPathR2, Float64, Float64}} = []
        total_time = 0
        params = [parameters.v_min, parameters.v_max, parameters.a_max, -parameters.a_min]

        for i in eachindex(configurations)
            next = i == length(configurations) ? 1 : i + 1

            n_i, v_i, h_i = configurations[i]
            n_f, v_f, h_f = configurations[next]

            starting::Vector{Float64} = [locations[n_i][1], locations[n_i][2], headings[h_i]]
            ending::Vector{Float64} = [locations[n_f][1], locations[n_f][2], headings[h_f]]

            path, time, _ = AcceleratedDubins.fastest_path(starting, ending, radii, params, [speeds[v_i], speeds[v_f]])

            push!(full_path, (path, speeds[configurations[i][2]], speeds[configurations[next][2]]))
            total_time += time
        end

        return full_path, total_time
    end

    function configuration_time(graph::Array{Float64,6}, configurations::Vector{Tuple{Int64, Int64, Int64}})
        total_time = 0
        for i in eachindex(configurations)
            next = i == length(configurations) ? 1 : i + 1
            n_i, v_i, h_i = configurations[i]
            n_f, v_f, h_f = configurations[next]
            total_time += graph[n_i, n_f, v_i, v_f, h_i, h_f]
        end
        return total_time
    end

end