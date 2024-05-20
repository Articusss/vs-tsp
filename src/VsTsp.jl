module VsTsp
    import IterTools as itr
    import Plots as plt
    using Random
    include("AcceleratedDubins.jl")
    include("Visual.jl")
    include("Helper.jl")

    struct VehicleParameters
        v_min::Float64
        v_max::Float64
        a_min::Float64
        a_max::Float64
        r_min::Float64
        r_max::Float64
    end

    function test(graph, speed, heading, radii)
        points = Helper.read_tsp_file("./instances/berlin52.tsp")
        params = VehicleParameters(30., 67., -3., 2., 65.7, 264.2)

        config, sequence, base_time = cheapest_insertion(graph)
        path, path_time = Helper.retrieve_path(points, config, params, speed, heading, radii)
        calc_time = Helper.configuration_time(graph, config)
        println((base_time, path_time, calc_time))
    end

    function test2(graph, speeds, headings, radii)
        locations = Helper.read_tsp_file("./instances/berlin52.tsp")
        parameters = VehicleParameters(30., 67., -3., 2., 65.7, 264.2)
        params = [parameters.v_min, parameters.v_max, parameters.a_max, -parameters.a_min]

        node_i, node_f = 1,2
        h_i, h_f = 1,2
        v_i, v_f = 1,2
        start::Vector{Float64} = [locations[node_i][1], locations[node_i][2], headings[h_i]]
        stop::Vector{Float64} = [locations[node_f][1], locations[node_f][2], headings[h_f]]
        path, time, _ = AcceleratedDubins.fastest_path(start, stop, radii, params, [speeds[v_i], speeds[v_f]])
        println((time, graph[node_i, node_f, v_i, v_f, h_i, h_f]))
        println(path.lengths)

        path, time, _ = AcceleratedDubins.fastest_path(stop, start, radii, params, [speeds[v_f], speeds[v_i]])
        println((time, graph[node_f, node_i, v_f, v_i, h_f, h_i]))
    
    end

    function getGraph()
        points = Helper.read_tsp_file("./instances/berlin52.tsp")
        params = VehicleParameters(30., 67., -3., 2., 65.7, 264.2)
        num_speeds = 4

        return compute_trajectories(points, params, num_speeds)
    end

    function get_params()
        return VehicleParameters(30., 67., -3., 2., 65.7, 264.2)
    end

    #TODO not calculate matrix diagonals (distance from the point to itself) can speed up
    #Also maybe take advantage of paths being borderline simmetric (only issue is radii and acceleration)
    function compute_trajectories(locations::Array{Tuple{Float64, Float64}, 1}, parameters::VehicleParameters, num_speeds::Int64, num_headings::Int64 = 8, radii_samples::Int64 = 8)
        #Build a 6-dimensional graph (starting node, ending node, starting speed, ending speed, starting heading angle, ending heading angle)
        num_locations = length(locations)
        graph::Array{Float64, 6} = fill(-1, (num_locations, num_locations, num_speeds, num_speeds, num_headings, num_headings))

        #Split heading angles and speeds according to their numbers
        #NOTE - USING MAX SPEED OF R_MAX, otherwise all maximum speeds will be invalid
        speeds = collect(range(parameters.v_min, AcceleratedDubins.speed_by_radius(parameters.r_max), num_speeds))
        headings = collect(range(0, 2 * pi, num_headings))

        radii = AcceleratedDubins.radii_samples_exp(parameters.r_min, parameters.r_max, radii_samples)

        # v_min v_max a_max -a_min
        params = [parameters.v_min, parameters.v_max, parameters.a_max, -parameters.a_min]
        
        for (node_i, node_f) in itr.product(1:num_locations, 1:num_locations)
            for (v_i, v_f) in itr.product(1:num_speeds, 1:num_speeds)
                for (h_i, h_f) in itr.product(1:num_headings, 1:num_headings)
                    start::Vector{Float64} = [locations[node_i][1], locations[node_i][2], headings[h_i]]
                    stop::Vector{Float64} = [locations[node_f][1], locations[node_f][2], headings[h_f]]
                    path, time, _ = AcceleratedDubins.fastest_path(start, stop, radii, params, [speeds[v_i], speeds[v_f]])
                    #Set to edge, note that we can't take advantage of simmetry because acceleration max/min is not necessarily the same
                    graph[node_i,node_f,v_i,v_f,h_i,h_f] = path === nothing ? Inf : time
                end
            end
        end

        return graph, speeds, headings, radii
    end

    function shortest_time_by_sequence(graph::Array{Float64, 6}, sequence::Vector{Int64})
        num_headings = size(graph,5)
        num_speeds = size(graph, 3)

        best = Inf
    
        #All possible starting speed/heading angles
        for (speed_start, heading_start) in itr.product(1:num_speeds, 1:num_headings)
            
            #Holds best path of arbitrary position considering (speed, headingAngle)
            prev::Array{Tuple{Float64, Bool}, 2} = [(graph[sequence[1], sequence[2], speed_start, i, heading_start, j], false) for i in 1:num_speeds, j in 1:num_headings]
            curr::Array{Tuple{Float64, Bool}, 2} = fill((0, false), (num_speeds, num_headings))

            #Validates value, to remove necessity of creating new matrix everytime this runs
            valid = true
            local_best_time = Inf

            for (starting_speed, starting_heading) in itr.product(1:num_speeds, 1:num_headings)
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
                    valid = pos % 2 == 0 ? !valid : valid
                end

                #Now connect to start again
                for prev_speed in 1:num_speeds
                    for prev_heading in 1:num_headings
                        val = prev[prev_speed, prev_heading][1] + graph[sequence[length(sequence)], sequence[1], prev_speed, starting_speed, prev_heading, starting_heading]
                        if val < local_best_time
                            local_best_time = val
                        end
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

    function cheapest_insertion(graph::Array{Float64,6})
        #Cheapest insertion algorithm
        num_headings = size(graph,5)
        num_speeds = size(graph, 3)
        
        #TODO follow standard cheapest insertion and start with vertex with lowest time?
        #Connect first two points in the best configuration possible
        best_starting, best_ending, total_time = Helper.find_lowest_edge(graph, 1, 2, num_headings, num_speeds)
        to_add = Set{Int64}(3:size(graph,1))

        #Vector corresponding to (point, speed, heading angle) of each point in the sequence
        full_path::Vector{Tuple{Int64, Int64, Int64}} = [best_starting, best_ending]

        while !isempty(to_add)
            best_configuration::Tuple{Int64, Int64, Int64} = (-1, -1, -1)
            best_time = Inf
            best_position::Int64 = -1

            #Find cheapest insertion -> Try to insert each node between every point in solution, select the best position and insert
            for candidate in to_add
                for pos in eachindex(full_path)
                    configuration, time = Helper.find_lowest_edge_between(graph, full_path[pos], full_path[pos == length(full_path) ? 1 : pos + 1], candidate, num_headings, num_speeds)
                    if time < best_time
                        best_time = time
                        best_position = pos
                        best_configuration = configuration
                    end
                end
            end

            #Insert best insertion
            delete!(to_add, best_configuration[1])
            insert!(full_path, best_position + 1, best_configuration)
            total_time += best_time
        end

        #Get sequence too to facilitate future algorithms
        sequence::Vector{Int64} = [x[1] for x in full_path]
        return full_path, sequence, total_time
    end

end # module VsTsp