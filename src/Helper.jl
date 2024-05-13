module Helper
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
end