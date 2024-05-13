module Visual
    import PyPlot as plt
    using ..AcceleratedDubins

    function plot_dubins_curve(path)
        confx, confy = AcceleratedDubins.sample_path(path)
        plt.plot(confx, confy)
    end

    function plot_speeds(path, max_min_speed, initial_final_speeds)
        times, speeds = AcceleratedDubins.speed_profile(path, max_min_speed, initial_final_speeds)
        plt.plot(times, speeds)
    end
end