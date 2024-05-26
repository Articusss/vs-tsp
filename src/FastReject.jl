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
            #Calculate manually TODO fix 
            else
                for _ in 1:window_size
                    #Find smallest edge between curr_config and next
                    curr_config, time = Helper.find_lowest_edge(graph, curr_config, ending, num_headings, num_speeds)
                    window_sum += time

                    #Advance
                    next_idx = next_idx == length(sequence) ? 1 : next_idx + 1
                    pos = next
                    next = sequence[next_idx]
                end
            end

            total_sum += window_sum
        end

        return (1. / window_size) * total_sum
    end
end