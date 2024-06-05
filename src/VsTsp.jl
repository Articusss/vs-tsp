module VsTsp
    import IterTools as itr
    import Plots as plt
    using Random
    include("AcceleratedDubins.jl")
    include("Visual.jl")
    include("Helper.jl")
    include("Vns.jl")
    include("FastReject.jl")

    struct VehicleParameters
        v_min::Float64
        v_max::Float64
        a_min::Float64
        a_max::Float64
        r_min::Float64
        r_max::Float64
    end

    function get_graph(instance_path::String, params::VehicleParameters, num_speeds::Int64 = 4, num_headings::Int64 = 8, radii_samples::Int64 = 8)
        points = Helper.read_tsp_file(instance_path)

        return compute_trajectories(points, params, num_speeds, num_headings, radii_samples)
    end

    function get_cessna172_params()
        return VehicleParameters(30., 67., -3., 2., 65.7, 264.2)
    end

    function vns_tsp(graph::Array{Float64, 6}, instance_path::String, params::VehicleParameters, speeds::Vector{Float64}, headings::Vector{Float64}, radii::Vector{Float64}, max_iterations::Int64 = 2000, num_windows::Int64 = 4, verbose::Bool = false)
        _, initial_seq, _ = cheapest_insertion(graph)
        vns_seq, time = variable_neighborhood_search_fast(graph, initial_seq, max_iterations, num_windows, verbose)

        points = Helper.read_tsp_file(instance_path)
        _, config = Helper.shortest_configuration_by_sequence(graph, vns_seq)
        path, _ = Helper.retrieve_path(points, config, params, speeds, headings, radii)

        return path, points, config
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

    function cheapest_insertion(graph::Array{Float64,6})
        #Cheapest insertion algorithm
        num_headings = size(graph,5)
        num_speeds = size(graph, 3)
        
        #TODO follow standard cheapest insertion and start with vertex with lowest time?
        #Connect first two points in the best configuration possible
        best_starting, best_ending, total_time = Helper.find_lowest_2_tour(graph, 1, 2, num_headings, num_speeds)
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

    #BASE VNS -> no fast reject
    function variable_neighborhood_search(graph::Array{Float64, 6}, initial_sequence::Vector{Int64}, max_iterations::Int64, verbose::Bool = false)
        best_time = Helper.shortest_time_by_sequence(graph, initial_sequence)
        best_sequence = deepcopy(initial_sequence)
        len = length(initial_sequence)

        for i in 1:max_iterations
            if verbose
                println((i, best_time))
            end
            #Shake
            local_sequence = rand([Vns.path_exchange, Vns.path_move])(deepcopy(best_sequence))
            local_time = Helper.shortest_time_by_sequence(graph, local_sequence)

            #Search
            for j in 1:len^2
                search = Vns.rand([Vns.one_point_move, Vns.open_point_exchange])(deepcopy(local_sequence))
                search_time = Helper.shortest_time_by_sequence(graph, search)

                if search_time < local_time
                    local_time = search_time
                    local_sequence = search
                end
            end

            if local_time < best_time
                best_time = local_time
                best_sequence = local_sequence
            end
        end


        return best_sequence, best_time
    end

    #VNS with fast reject
    function variable_neighborhood_search_fast(graph::Array{Float64, 6}, initial_sequence::Vector{Int64}, max_iterations::Int64, num_windows::Int64, verbose::Bool = false)
        best_time = Helper.shortest_time_by_sequence(graph, initial_sequence)
        best_sequence = deepcopy(initial_sequence)
        len = length(initial_sequence)
        stored = fill(-1., (len,len))

        for i in 1:max_iterations
            if verbose
                println((i, best_time))
            end
            #Shake
            local_sequence = rand([Vns.path_exchange, Vns.path_move])(deepcopy(best_sequence))
            local_costs::Vector{Float64} = fill(-1., num_windows)
            
            #Search
            for j in 1:len^2
                search = Vns.rand([Vns.one_point_move, Vns.open_point_exchange])(deepcopy(local_sequence))

                reject, costs = FastReject.fast_reject(graph, search, local_sequence, local_costs, num_windows, stored, best_time)
                if !reject
                    local_sequence = search
                    local_costs = costs
                end
            end

            local_time = Helper.shortest_time_by_sequence(graph, local_sequence)
            if local_time < best_time
                best_time = local_time
                best_sequence = local_sequence
            end
        end


        return best_sequence, best_time
    end
end # module VsTsp