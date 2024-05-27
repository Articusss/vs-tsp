module FastReject
    using ..Helper
    import IterTools as itr
    function window_cost(graph::Array{Float64, 6}, sequence::Vector{Int64}, window_size::Int64, stored::Array{Float64, 2})
        num_headings = size(graph,5)
        num_speeds = size(graph, 3)
        total_sum = 0.
        for pos_idx in eachindex(sequence)
            window_sum = 0.
            next_idx = pos_idx == length(sequence) ? 1 : pos_idx + 1

            pos = sequence[pos_idx]
            next = sequence[next_idx]
            
            #Use hash table
            if window_size + 1 == ndims(stored)
                #Check if first element is already calculated, if not, calculate it
                if stored[pos,next] == -1
                    _, _, window_sum = Helper.find_lowest_edge(graph, pos, next, num_headings, num_speeds)
                    #Store
                    stored[pos,next] = window_sum
                else
                    window_sum = stored[pos,next]
                end
            #Calculate manually using open loop forward search
            else
                #Holds best path of arbitrary position considering (speed, headingAngle)
                prev::Array{Tuple{Float64, Bool}, 2} = fill((0., false), (num_speeds, num_headings))
                curr::Array{Tuple{Float64, Bool}, 2} = fill((0., false), (num_speeds, num_headings))

                #Validates value, to remove necessity of creating new matrix everytime this runs
                valid = true
                for i in 1:window_size
                    for (prev_speed, curr_speed) in itr.product(1:num_speeds, 1:num_speeds)
                        for (prev_heading, curr_heading) in itr.product(1:num_headings, 1:num_headings)
                            #Not valid, assign first value for future comparisons
                            if curr[curr_speed, curr_heading][2] != valid
                                curr[curr_speed, curr_heading] = (prev[prev_speed,prev_heading][1] + graph[pos, next, prev_speed, curr_speed, prev_heading, curr_heading], valid)
                            else
                                #Valid, compare with previously calculated value
                                val = prev[prev_speed,prev_heading][1] + graph[pos, next, prev_speed, curr_speed, prev_heading, curr_heading]
                                if val < curr[curr_speed, curr_heading][1]
                                    curr[curr_speed, curr_heading] = (val, valid)
                                end
                            end
                        end
                    end
                    
                    #Swap, swap validity if necessary
                    curr, prev = prev, curr
                    valid = i % 2 == 0 ? !valid : valid

                    #Advance
                    next_idx = next_idx == length(sequence) ? 1 : next_idx + 1
                    pos = next
                    next = sequence[next_idx]
                end

                window_sum = minimum(prev)[1]
            end

            total_sum += window_sum
        end

        return (1. / window_size) * total_sum
    end

    function fast_reject(graph::Array{Float64}, candidate_sequence::Vector{Int64}, original_sequence::Vector{Int64}, num_windows::Int64, stored::Array{Float64, 2}, best_time::Float64, exact_comparison::Bool = false)
        max_window = num_windows
        for i in 1:num_windows
            if window_cost(graph, candidate_sequence, 2^(i-1), stored) > best_time
                max_window = i
                break
            end
        end

        for j in 1:max_window
            if exact_comparison && j == max_window
                return Helper.shortest_time_by_sequence(graph, candidate_sequence) > Helper.shortest_time_by_sequence(graph, original_sequence)
            end

            cand_cost = window_cost(graph, candidate_sequence, 2^(j-1), stored)
            og_cost = window_cost(graph, original_sequence, 2^(j-1), stored)
            if cand_cost > og_cost
                return true
            end
        end

        return false
    end
end