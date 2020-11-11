function outputfig!(f::Fluid, s::Structure; varname::String = "rho", axis::Int = 1, frame::Int = 1, figpath::String = "../outputfig/", plotdim::String = "2D", levels::Vector=[0,1,1])
    ioff()
    if plotdim in ("2D", "2D-x", "2D-y")
        x = [f.d[1]*(i-0.5-f.ng)+f.point1[1]  for  i=1:f.nmesh[1]+f.ng*2]
        y = [f.d[2]*(i-0.5-f.ng)+f.point1[2]  for  i=1:f.nmesh[2]+f.ng*2]
        if varname in ("rho","p","e")
            data = Array{Float64}(undef, Tuple(f.nmesh .+  2 * f.ng ))
            for k in eachindex(f.cells)
                data[k] = getfield(f.cells[k], Symbol(varname))
            end
        elseif varname == "u"
            data = Array{Float64}(undef, Tuple(f.nmesh .+  2 * f.ng ))
            for k in eachindex(f.cells)
                data[k] = f.cells[k].u[axis]
            end  
        else
            error("undef varname")          
        end  
        
        if plotdim == "2D-x"
            FVM.output_plot_func(x, data[:, floor(Int, size(data, 2)/2)], levels = levels)
        elseif plotdim == "2D-y"
            FVM.output_plot_func(y, data[floor(Int, size(data, 1)/2), :], levels = levels)  
        else        
            FVM.output_plot_func(x, y, data, levels = levels)
        end
    elseif plotdim == "1D"
        x = [f.d[1]*(i-0.5-f.ng)+f.point1[1]  for  i=1:f.nmesh[1]+f.ng*2]
        data = Array{Float64}(undef, length(x))
        for k in eachindex(f.cells)
            data[k] = getfield(f.cells[k], Symbol(varname))
        end
        FVM.output_plot_func(x, data, levels = levels)        
    else
        error("undef plotdim")
    end
    
    # if plotdim == "2D"
    #     for e in s.mesh.edge
    #         if typeof(s) == RigidStructure 
    #             RS.output_plot_func(e)
    #         else
    #             @warn("undef structure type")
    #         end
    #     end
    # end
    savefig(figpath*varname*"_"*string(frame+FRAME_BASE)*".png",dpi=100)
    close(1)
end

function output!(frame::Int, time::Float64; filepath::String = "../outputdata/")
    if frame == 0
        open(filepath*"time.txt","w") do f
            writedlm(f,Union{Float64,Int}[frame time])
        end
    else
        open(filepath*"time.txt","a") do f
            writedlm(f,Union{Float64,Int}[frame time])
        end
    end
end