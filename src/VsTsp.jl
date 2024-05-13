module VsTsp
    import PyPlot as plt
    import IterTools as itr
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

        compute_trajectories(points, params, num_speeds)
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
                        stop::Vector{Float64} = [locations[node_f][1], locations[node_f][1], headings[h_f]]
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

end # module VsTsp