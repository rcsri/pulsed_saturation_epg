function y_norm = normalize_data_to_steady_state_of_sim(y, sim_norm, sequence)

if isfield(sequence,'fsat_vect')
    y_norm = zeros(size(y));
    for ind_scan = 1:length(sequence.b1_nom_vect)
        curr_sim_norm = sim_norm{ind_scan};        
        curr_mzss    = curr_sim_norm(sequence.index_norm_sim_vect{ind_scan});
        y_norm(:,ind_scan)  = curr_mzss * y (:, ind_scan) / y(sequence.index_norm_acq, ind_scan);
    end
else    
    mzss = sim_norm(sequence.index_norm_sim, :); 
    y_norm      = bsxfun(@times, mzss, bsxfun(@times, y (:, :), 1./y(sequence.index_norm_acq, :)));   
end
