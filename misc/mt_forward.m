function [signal_out] = mt_forward(x, sequence, current_pool, lstype)
dfa = 0;
if strcmp(current_pool,'mt')
    dfs = 0;
end
M0a = 1;
R1a = x(1);
R2a = x(2);
R1s = x(3);
R2s = x(4);
M0s = x(5);
Rs  = x(6);
if numel(x) >= 7
    fit_b0 = 1;
    b0_in = x(7:end);
else
    fit_b0 = 0;
end

%------------%
% BM sim
%------------%
p0 = bmsim_params_0CEST_1MTC(dfa, R1a, R2a, M0a, dfs, R1s, R2s, M0s, Rs);  % set up pools
n_scans = length(sequence.b1_nom_vect);
signal_out = zeros(length(sequence.fsat), n_scans);
for ind_scan = 1:n_scans
    
    sequence.b1_sat = sequence.b1_nom_vect(ind_scan);  % B1 amplitude (uT)
    
    num_sat = round(sequence.ts_nom_vect(ind_scan)/(1e3*sum(sequence.t_sat))); % e.g."900ms"-->num_sat = 900/10ms = 90. "450ms"--> num_sat=45
    sequence.num_sat = num_sat; % e.g."900ms"-->num_sat = 900/10ms = 90. "450ms"--> num_sat=45
    
    % delay at the end of sat-TR, different depending on the duration of sat (either 45 blocks or 90 blocks)
    dur_block = sum(sequence.t_sat) + sequence.tspoil_sat;
    n_blocks = 1; % number of blocks within each sat block. In MT sequence, each of the (90 or 45) sat blocks are not further divided into smaller blocks so n_blocks = 1
    tdel_end = (sequence.TR - sequence.num_sat*((n_blocks*dur_block)+sequence.tspoil_sat)) - sequence.tdel_pre;
    sequence.tdel_end = tdel_end;
    
    if fit_b0
        sequence.b0 = b0_in(ind_scan);
    else
        sequence.b0 = sequence.b0_vect(ind_scan);
    end
    sequence.b1_scale = sequence.b1_scale_vect(ind_scan);
    
    [~, Mxy_out] = bmsim_sequence_epg_ver2(p0, sequence, lstype);
        
    signal_out(:,ind_scan) = Mxy_out;
    
end
