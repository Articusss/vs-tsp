module VsTsp
    import PyPlot as plt
    import IterTools as itr
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

    function test()
        params = [5., 15., 2.6, 4] ## v_min v_max a_max -a_min
        speeds = [7.5, 12.5]
        r_min = 20. # minimum turning radius
        r_max = 90. # maximum turning radius
        start = [0, 0, 2*pi*rand()]
        stop = [100, 100, 2*pi*rand()]
        path, errcode = AcceleratedDubins.fastest_path(start, stop, AcceleratedDubins.radii_samples_exp(r_min, r_max, 3), params, speeds)
        #Visual.plot_dubins_curve(path)
        Visual.plot_speeds(path, params, speeds)
    end

    function test2()
        points = Helper.read_tsp_file("./instances/berlin52.tsp")
        params = VehicleParameters(30., 67., -3., 2., 65.7, 264.2)
        num_speeds = 4

        return compute_trajectories(points[1:3], params, num_speeds)
    end

    function compute_trajectories(locations::Array{Tuple{Float64, Float64}, 1}, parameters::VehicleParameters, num_speeds::Int64, num_headings::Int64 = 8, radii_samples::Int64 = 8)
        #Build a 6-dimensional graph (starting node, ending node, starting speed, ending speed, starting heading angle, ending heading angle)
        num_locations = length(locations)
        graph::Array{Float64, 6} = fill(-1, (num_locations, num_locations, num_speeds, num_speeds, num_headings, num_headings))

        #Split heading angles and speeds according to their numbers
        speeds = collect(range(parameters.v_min, parameters.v_max, num_speeds))
        headings = collect(range(0, 2 * pi, num_headings))

        radii = AcceleratedDubins.radii_samples_exp(parameters.r_min, parameters.r_max, radii_samples)

        # v_min v_max a_max -a_min
        params = [parameters.v_min, parameters.v_max, parameters.a_max, -parameters.a_min]
        
        for (node_i, node_f) in itr.product(1:num_locations, 1:num_locations)
            for (v_i, v_f) in itr.product(1:num_speeds, 1:num_speeds)
                for (h_i, h_f) in itr.product(1:num_headings, 1:num_headings)
                    #Check if already calculated
                    if(graph[node_i,node_f,v_i,v_f,h_i,h_f] == -1)
                        start::Vector{Float64} = [locations[node_i][1], locations[node_i][2], headings[h_i]]
                        stop::Vector{Float64} = [locations[node_f][1], locations[node_f][2], headings[h_f]]
                        path, _ = AcceleratedDubins.fastest_path(start, stop, radii, params, [speeds[v_i], speeds[v_f]])
                        #Set to edge and take advantage of simmetry
                        if path === nothing
                            graph[node_i,node_f,v_i,v_f,h_i,h_f] = Inf
                            graph[node_f,node_i,v_f,v_i,h_f,h_i] = Inf
                        else
                            time = AcceleratedDubins.path_time(path, params, [speeds[v_i], speeds[v_f]])
                            graph[node_i,node_f,v_i,v_f,h_i,h_f] = time
                            graph[node_f,node_i,v_f,v_i,h_f,h_i] = time
                        end

                    end
                end
            end
        end

        return graph
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

    function initialization(graph::Array{Float64,6})
        #Insert the first 3 locations randomly
        sequence::Vector{Int64} = shuffle(collect(1:3))
        #Vector corresponding to (heading angle, speed)
        configurations::Vector{Tuple{Int64, Int64}}


        return sequence
    end

end # module VsTsp