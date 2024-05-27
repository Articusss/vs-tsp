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

    function find_lowest_2_tour(graph::Array{Float64, 6}, starting::Int64, ending::Int64, num_headings::Int64, num_speeds::Int64)
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

    function find_lowest_edge(graph::Array{Float64, 6}, starting::Int64, ending::Int64, num_headings::Int64, num_speeds::Int64)
        best = Inf

        best_starting::Tuple{Int64, Int64, Int64} = (starting, -1, -1)
        best_ending::Tuple{Int64, Int64, Int64} = (ending, -1, -1)

        for (v_i, v_f) in itr.product(1:num_speeds, 1:num_speeds)
            for (h_i, h_f) in itr.product(1:num_headings, 1:num_headings)
                total_time = graph[starting, ending, v_i, v_f, h_i, h_f]
                if total_time < best
                    best = total_time
                    best_starting = (starting, v_i, h_i)
                    best_ending = (ending, v_f, h_f)
                end
            end
        end

        return best_starting, best_ending, best
    end

    #Lowest edge with starting config
    function find_lowest_edge(graph::Array{Float64, 6}, starting::Tuple{Int64, Int64, Int64}, ending::Int64, num_headings::Int64, num_speeds::Int64)
        best_time = Inf
        best_ending::Tuple{Int64, Int64, Int64} = (ending, -1, -1)

        for v_f in 1:num_speeds
            for h_f in 1:num_headings
                total_time = graph[starting[1], ending, starting[2], v_f, starting[3], h_f]
                if total_time < best_time
                    best_time = total_time
                    best_ending = (ending, v_f, h_f)
                end
            end
        end

        return best_ending, best_time
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

    function shortest_time_by_sequence(graph::Array{Float64, 6}, sequence::Vector{Int64})
        num_headings = size(graph,5)
        num_speeds = size(graph, 3)

        best = Inf
    
        #All possible starting speed/heading angles
        for (speed_start, heading_start) in itr.product(1:num_speeds, 1:num_headings)
            
            #Holds best path of arbitrary position considering (speed, headingAngle)
            prev::Array{Tuple{Float64, Bool}, 2} = [(graph[sequence[1], sequence[2], speed_start, i, heading_start, j], false) for i in 1:num_speeds, j in 1:num_headings]
            curr::Array{Tuple{Float64, Bool}, 2} = fill((0., false), (num_speeds, num_headings))

            #Validates value, to remove necessity of creating new matrix everytime this runs
            valid = true
            local_best_time = Inf

            #Update best path for pos + 1, starting from second position
            for pos in 2:(length(sequence)-1)
                curr_node = sequence[pos + 1]
                prev_node = sequence[pos]
                for (prev_speed, curr_speed) in itr.product(1:num_speeds, 1:num_speeds)
                    for (prev_heading, curr_heading) in itr.product(1:num_headings, 1:num_headings)
                        #Not valid, assign first value for future comparisons
                        if curr[curr_speed, curr_heading][2] != valid
                            curr[curr_speed, curr_heading] = (prev[prev_speed,prev_heading][1] + graph[prev_node, curr_node, prev_speed, curr_speed, prev_heading, curr_heading], valid)
                        else
                            #Valid, compare with previously calculated value
                            val = prev[prev_speed,prev_heading][1] + graph[prev_node, curr_node, prev_speed, curr_speed, prev_heading, curr_heading]
                            if val < curr[curr_speed, curr_heading][1]
                                curr[curr_speed, curr_heading] = (val, valid)
                            end
                        end
                    end
                end
                
                #Swap, swap validity if necessary
                curr, prev = prev, curr
                valid = pos % 2 == 1 ? !valid : valid
            end

            #Now connect to start again
            for prev_speed in 1:num_speeds
                for prev_heading in 1:num_headings
                    val = prev[prev_speed, prev_heading][1] + graph[sequence[length(sequence)], sequence[1], prev_speed, speed_start, prev_heading, heading_start]
                    if val < local_best_time
                        local_best_time = val
                    end
                end
            end

            #Check with global best
            if local_best_time < best
                best = local_best_time
            end
        end

        return best
    end

    #Same as shortest_time_by_sequence, but returns configuration too
    function shortest_configuration_by_sequence(graph::Array{Float64, 6}, sequence::Vector{Int64})
        num_headings = size(graph,5)
        num_speeds = size(graph, 3)

        best = Inf
        best_config = []
    
        #All possible starting speed/heading angles
        for (speed_start, heading_start) in itr.product(1:num_speeds, 1:num_headings)
            
            #Holds best path of arbitrary position considering (speed, headingAngle)
            prev::Array{Tuple{Float64, Bool}, 2} = [(graph[sequence[1], sequence[2], speed_start, i, heading_start, j], false) for i in 1:num_speeds, j in 1:num_headings]
            curr::Array{Tuple{Float64, Bool}, 2} = fill((0., false), (num_speeds, num_headings))

            prev_configs::Matrix{Vector{Tuple{Int64, Int64, Int64}}} = [[(sequence[1], speed_start, heading_start), (sequence[2], i, j)] for i in 1:num_speeds, j in 1:num_headings]
            curr_configs = Matrix{Vector{Tuple{Int64, Int64, Int64}}}(undef, num_speeds, num_headings)

            #Validates value, to remove necessity of creating new matrix everytime this runs
            valid = true
            local_best_time = Inf
            local_best_config = []

            #Update best path for pos + 1, starting from second position
            for pos in 2:(length(sequence)-1)
                curr_node = sequence[pos + 1]
                prev_node = sequence[pos]
                for (prev_speed, curr_speed) in itr.product(1:num_speeds, 1:num_speeds)
                    for (prev_heading, curr_heading) in itr.product(1:num_headings, 1:num_headings)
                        #Not valid, assign first value for future comparisons
                        if curr[curr_speed, curr_heading][2] != valid
                            curr[curr_speed, curr_heading] = (prev[prev_speed,prev_heading][1] + graph[prev_node, curr_node, prev_speed, curr_speed, prev_heading, curr_heading], valid)
                        
                            curr_configs[curr_speed, curr_heading] = vcat(prev_configs[prev_speed, prev_heading], [(curr_node, curr_speed, curr_heading)])
                        else
                            #Valid, compare with previously calculated value
                            val = prev[prev_speed,prev_heading][1] + graph[prev_node, curr_node, prev_speed, curr_speed, prev_heading, curr_heading]
                            if val < curr[curr_speed, curr_heading][1]
                                curr[curr_speed, curr_heading] = (val, valid)
                                curr_configs[curr_speed, curr_heading] = vcat(prev_configs[prev_speed, prev_heading], [(curr_node, curr_speed, curr_heading)])
                            end
                        end
                    end
                end
                
                #Swap, swap validity if necessary
                curr, prev = prev, curr
                curr_configs, prev_configs = prev_configs, curr_configs
                valid = pos % 2 == 1 ? !valid : valid
            end

            #Now connect to start again
            for prev_speed in 1:num_speeds
                for prev_heading in 1:num_headings
                    val = prev[prev_speed, prev_heading][1] + graph[sequence[length(sequence)], sequence[1], prev_speed, speed_start, prev_heading, heading_start]
                    if val < local_best_time
                        local_best_time = val
                        local_best_config = prev_configs[prev_speed, prev_heading]
                    end
                end
            end

            #Check with global best
            if local_best_time < best
                best = local_best_time
                best_config = local_best_config
            end
        end

        return best, best_config
    end
end