module Helper
    import IterTools as itr

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

        for (starting_speed, ending_speed) in itr.product(1:num_speeds, 1:num_speeds)
            for (starting_heading, ending_heading) in itr.product(1:num_headings, 1:num_headings) 
                if graph[starting, ending, starting_speed, ending_speed, starting_heading, ending_heading] < best
                    best = graph[starting, ending, starting_speed, ending_speed, starting_heading, ending_heading]
                    best_starting = (starting, starting_speed, starting_heading)
                    best_ending = (ending, ending_speed, ending_heading)
                end
            end
        end

        return best_starting, best_ending
    end

    function find_lowest_edge_between(graph::Array{Float64}, starting::Tuple{Int64, Int64, Int64}, ending::Tuple{Int64, Int64, Int64}, middle::Int64, num_headings::Int64, num_speeds::Int64)
        best_time = Inf
        best_middle::Tuple{Int64, Int64, Int64} = (middle, -1, -1)

        #Starting and ending configurations are already defined, need to define middle configurations
        for speed in 1:num_speeds
            for heading in 1:num_headings
                dist = graph[starting[1], middle, starting[2], speed, starting[3], heading] + graph[middle, ending[1], speed, ending[2], heading, ending[3]]
                if dist < best_time
                    best_time = dist
                    best_middle = (middle, speed, heading)
                end
            end
        end

        best_middle, best_time
    end

end