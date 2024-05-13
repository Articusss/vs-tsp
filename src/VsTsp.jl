module VsTsp
    import PyPlot as plt
    include("AcceleratedDubins.jl")
    include("Visual.jl")

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

end # module VsTsp