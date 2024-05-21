module Vns
    using Random

    function path_move(sequence::Vector{Int64})
        len = length(sequence)
        if len <= 1
            return sequence
        else
            start_idx = rand(1:len)
            end_idx = rand(start_idx:len)

            segment = sequence[start_idx:end_idx]
            removed = vcat(sequence[1:start_idx-1], sequence[end_idx+1:len])
            insertion_point = rand(0:len - (end_idx - start_idx + 1))

            res = insertion_point == 0 ? vcat(segment, removed) : vcat(removed[1:insertion_point], segment, removed[insertion_point+1:length(removed)])

            return res
        end
    end

    function path_exchange(sequence::Vector{Int64})
        len = length(sequence)
        if len < 2
            return sequence
        end

        start1 = rand(1:len)
        start2 = rand(1:len)
        #Start2 must be different
        while start2 == start1
            start2 = rand(1:len)
        end
        #Start1 should be smaller
        if start2 < start1
            start1, start2 = start2, start1
        end

        end1 = rand(start1:start2-1)
        end2 = rand(start2:len)

        seg1 = sequence[start1:end1]
        seg2 = sequence[start2:end2]

        res = vcat(sequence[1:start1-1], seg2, sequence[end1 + 1:start2-1], seg1, sequence[end2+1:len])

        return res
    end

    function one_point_move(sequence::Vector{Int64})
        len = length(sequence)

        idx = rand(1:len)
        val = sequence[idx]
        deleteat!(sequence, idx)

        new_idx = rand(1:len)
        #Assure index is different
        while new_idx == idx
            new_idx = rand(1:len)
        end

        insert!(sequence, new_idx, val)

        return sequence
    end

    function open_point_exchange(sequence::Vector{Int64})
        len = length(sequence)

        idx1 = rand(1:len)
        idx2 = rand(1:len)

        while idx1 == idx2
            idx2 = rand(1:len)
        end

        sequence[idx1], sequence[idx2] = sequence[idx2], sequence[idx1]

        return sequence
    end
end