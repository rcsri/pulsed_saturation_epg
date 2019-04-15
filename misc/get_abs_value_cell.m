function sim_norm = get_abs_value_cell(sim_signal)

if iscell(sim_signal)
    sim_norm = cell(size(sim_signal)); 
    for ind = 1:length(sim_signal)
        sim_norm{ind} = abs(sim_signal{ind});
    end
else
    sim_norm = abs(sim_signal);
end
