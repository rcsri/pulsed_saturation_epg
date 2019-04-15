function [res, sim_signal, sim_norm_exclude_dummy, y_norm, R1W] = obj_qMT(x, y, sequence, fit_options, current_pool, t1obs, varargin)
% x: parameter values
% y: data
% t1obs: observed T1(s)
% sim_signal: full forward model including dummy TRs
% y_norm: data normalized by steady state signal without (unacquired) dummy TRs
% sim_norm_exclude_dummy: model, matched to data (excluding unacquired dummy TRs) 

%------------------------------------------%
% scaling factors to help with minimization 
%------------------------------------------%
if isempty(varargin) 
    x_in = x; 
    SCALE_Y = 1;
else
    scale_x = varargin{1};
    x_in = x./scale_x;
    SCALE_Y = varargin{2};
end

%-------------------------------%
% compute R1a (Henkelman 1993) 
%-------------------------------%
% Assume that first 6 parameters are: (R1A, R2A, R1B, R2B, M0B, RB)
R1obs = 1/t1obs; 
R1MTC   = x_in(3);
M0MTC   = x_in(5);
RMTC    = x_in(6);
x0 = 0.5; % (1/s)
lb = 0;
ub = Inf;
options = optimoptions('fmincon','Display','notify','Algorithm','active-set','MaxFunEvals',1e4);
R1W = fmincon(@(R1W) R1W_minus_R1obs_1MTC_0CEST(R1W,RMTC,M0MTC,R1MTC,R1obs),x0,[],[],[],[],lb,ub,[],options);
x_in(1) = R1W; 

%-----------------------------------------%
% x_in: current parameters 
% sim_signal: simulation signal 
% sim_norm: normalized simulation signal
%-----------------------------------------%
sim_signal = mt_forward(x_in, sequence, current_pool, fit_options.lstype); % [n freq sat, n powers]

if iscell(sim_signal)
    sim_norm = cellfun(@(x) x * SCALE_Y, sim_signal,'un',0); % 'un',0: shorthand way of writing 'UniformOutput'=false so that it outputs a cell rather than array
else
    sim_norm = SCALE_Y * get_abs_value_cell(sim_signal); 
end

%------------------------------------------------------------------%
% set steady state signal from data to that of sim, assumes M0z = 1
%------------------------------------------------------------------%
y_norm = normalize_data_to_steady_state_of_sim(y, sim_norm, sequence); 

%-------------------------------------------------------------------------------%
% select appropriate range for comparing sim and data (exclude scanner dummy TRs)
%-------------------------------------------------------------------------------%
if iscell(sim_norm)
    sim_norm_exclude_dummy = zeros(length(sim_norm{1})-sequence.n_dummy_vect(1),length(sim_norm));   
    for ind_scan = 1:length(sim_norm)
        curr_sim_norm = sim_norm{ind_scan};
        sim_norm_exclude_dummy(:,ind_scan) = curr_sim_norm(sequence.n_dummy_vect(ind_scan)+1:end); % note: no dummies but including 100k        
    end   
    
else
    sim_norm_exclude_dummy = sim_norm (sequence.n_dummy+1:end, :); % note: no dummies but including acquired 100k frequency offsets
end

%---------------------------------------------------------------%
% TFE shots: use every other simulated signal due to 2 TFE shots
%---------------------------------------------------------------%
if sequence.tfe.n_shots == 2
    sim_norm_exclude_dummy = sim_norm_exclude_dummy(1:2:end,:);
end

%----------------------------------------------------------%    
% include range: fit only frequency offset indices in m_vect 
%----------------------------------------------------------%
if isfield(fit_options, 'm_vect')
    m_range = fit_options.m_vect;
end
sim_norm_exclude_dummy = sim_norm_exclude_dummy(m_range, :);
y_norm = y_norm(m_range, :); 

%--------------------------------%
% W: weights for each freq offset
% res: residual
%--------------------------------%
if isfield(sequence, 'W') 
    if ~isempty(sequence.W)
        W = sequence.W; 
        W = W(:); 
        W = repmat(W, [1 size(y_norm,2)]);
        res = W(:) .* ( sim_norm_exclude_dummy(:) - y_norm(:) );    
    else
        res = (sim_norm_exclude_dummy(:) - y_norm(:));        
    end
else
    res = (sim_norm_exclude_dummy(:) - y_norm(:));
end