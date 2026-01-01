module CFSjul

greet() = print("Hello World!")


export iter_int
export iter_lie

function iter_int(utemp, dt, Ntrunc)

    num_input = size(utemp, 1)


    total_iterint = num_input*(1-num_input^Ntrunc)/(1-num_input)

    total_iterint = Int(total_iterint)

    Etemp = zeros(total_iterint, size(utemp, 2))

    ctrEtemp = zeros(Ntrunc+1)

    for i in 0:Ntrunc-1
        ctrEtemp[i+2] =  ctrEtemp[i+1] + num_input^(i+1)
    end

    #ctrEtemp[1] = 1

    sum_acc = cumsum(utemp, dims=2)*dt

    #Etemp[1:num_input,begin:end] = hcat(zeros(num_input, 1), sum_acc[:,begin:end-1])
    Etemp[1:num_input,:] = sum_acc

    for i in 1:Ntrunc-1
        
        start_prev_block = Int(ctrEtemp[i])
        end_prev_block = Int(ctrEtemp[i+1])
        end_current_block = Int(ctrEtemp[i+2])
        num_prev_block = end_prev_block - start_prev_block
        num_current_block = end_current_block - end_prev_block
        #print("i: ", i, " start_prev_block: ", start_prev_block, " end_prev_block: ", end_prev_block, " end_current_block: ", end_current_block, " num_prev_block: ", num_prev_block, " num_current_block: ", num_current_block)
        U_block = repeat(utemp, inner=(num_prev_block,1)) # inputs for current permutation
        #print("U_block size: ", size(U_block))
        prev_int_block = repeat(Etemp[start_prev_block+1:end_prev_block,:], outer=(num_input,1)) # block of previous permutations
        #print("prev_int_block size: ", size(prev_int_block))
        current_int_block = cumsum(U_block.*prev_int_block, dims = 2)*dt
        #print(size(current_int_block))
        #Etemp[end_prev_block+1:end_current_block,:] = hcat(zeros(num_current_block,1), current_int_block[:,begin:end-1])
        Etemp[end_prev_block+1:end_current_block,:] = current_int_block

    end

end



function iter_lie(h,vector_field,x,Ntrunc)
    num_vfield = size(vector_field, 2)

    # Compute total number of Lie derivatives
    total_lderiv = Int(num_vfield * (1 - num_vfield^Ntrunc) / (1 - num_vfield))

    # Allocate symbolic matrix for Lie derivatives
    Ltemp = Matrix{Num}(undef, total_lderiv, 1)

    # Track block boundaries
    ctrLtemp = zeros(Int, Ntrunc + 1)
    for i in 0:Ntrunc-1
        ctrLtemp[i+2] = ctrLtemp[i+1] + num_vfield^(i+1)
    end

    # First-order Lie derivatives
    for i in 1:num_vfield
        Ltemp[i, 1] = (Symbolics.jacobian([h], x) * vector_field[:, i])[1]
    end

    # Higher-order Lie derivatives
    for i in 1:Ntrunc-1
        start_prev_block = ctrLtemp[i] + 1
        end_prev_block = ctrLtemp[i+1]
        end_current_block = ctrLtemp[i+2]
        num_prev_block = end_prev_block - start_prev_block + 1

        # Track where to write new entries
        write_index = end_prev_block + 1

        for k in 1:num_vfield
            for j in 0:num_prev_block-1
                idx = start_prev_block + j
                Ltemp[write_index, 1] = (Symbolics.jacobian([Ltemp[idx, 1]], x) * vector_field[:, k])[1]
                write_index += 1
            end
        end
    end

end


end # module CFSjul
