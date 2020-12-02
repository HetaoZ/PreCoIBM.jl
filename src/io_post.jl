function save_time(frame, time, fname)
    if frame == 0
        open(fname*".txt", "w") do file
            writedlm(file, [frame time])
        end
    else
        open(fname*".txt", "a") do file
            writedlm(file, [frame time])
        end
    end
end