clear; clc;
close all; 
addpath(genpath('bmsim'));
addpath(genpath('epg'));
addpath(genpath('misc'));

%-----------------------------------------------------------%
% load simulated data and corresponding sequence parameters
%-----------------------------------------------------------%
load('y_mt.mat');        % simulated MT data
load('mt_sequence.mat'); % load sequence and other parameters 

%----------------------%
% fitting parameters
%----------------------%
pars.dt2a = 0.1000;
pars.max_iter = 1000;
pars.max_func_eval = 5000;
pars.alpha_ci = 0.6800;
pars.n_iter = 3;
pars.do_x_scaling = 1;
pars.optimizer_use_parallel = 1; 
pars.lstype = 'GS'; 

fit_options.reduce_dim_x        = 1;
fit_options.vary_t1a            = 0;
fit_options.vary_t2a            = 1;
fit_options.constrain_t2s       = 0;
fit_options.lstype              = 'GS'; % Gaussian lineshape for agar
fit_options.m_vect              = 2:32; % frequency offset indices include from y_mt

par_labels = {'T1a'    'T2a'    'T1s'    'T2s'    'M0s'    'Rs'};
lb_all = [12.5000    0.5882    1.0000    2.5000    0.1000         0];
ub_all = [12.5000   14.2857    1.0000   20.0000   50.0000    0.5000];
x0_all = [12.5000    1.4286    1.0000    6.6667    1.0000    0.1500];
scale_x = [10 0.1 1 0.0001 100 0.001];  % scale parameters for fitting 
scale_y = 1000;                         % scale data for fitting 
pool_in_mt = 'mt';
t1obs = 0.800; % observed T1 (s)

% oscillations
mt_sequence.W = ones(length(fit_options.m_vect),1);
indices_oscillation = 11:29; % which indices of fit_options.m_vect are oscillating
mt_sequence.W(11:29) = 1/numel(indices_oscillation); % lower weighting for samples that are closer together

%---------------------%
% objective function 
%---------------------%
F_full = @(x) obj_qMT (x, y_mt, mt_sequence, fit_options, pool_in_mt, t1obs, scale_x, scale_y);

if fit_options.reduce_dim_x % reduce number of free parameters    
    if fit_options.vary_t2a
       inds = [2, 4, 5, 6];   % vary: T2a, T2b, M0b, R
    else
       inds = [4, 5, 6];      % vary: T2b, M0b, R
    end     
    if fit_options.vary_t1a   % additionally vary: T1a
        inds = [1 inds];
    end
    par_labels_red = par_labels(inds);    
    lb_red       = lb_all(inds);
    ub_red       = ub_all(inds);
    x0_mt_red    = x0_all(inds);
    F = @(x) F_full( replace_param(x0_all, inds, x) );  
    scale_x_red = scale_x(inds); 
else
    F = F_full;
    x0_mt_red = x0_all;
    lb_red      = lb_all;
    ub_red      = ub_all;
    par_labels_red = par_labels;
    scale_x_red = scale_x;  
end

%-----------------------------%
% minimize
%-----------------------------%
[x_out_s, resnorm, resid, J] = run_lsqnonlin(F, x0_mt_red, lb_red, ub_red, pars);                  
x_out = x_out_s ./ scale_x_red;
if fit_options.reduce_dim_x
    xfit = x0_all ./ scale_x;
    for ind = 1:length(inds)
        p = inds(ind);
        xfit(p) = x_out(ind);
    end
else
    xfit = x_out;
end
[~, ~, ~, ~, R1W]  = obj_qMT (xfit, y_mt, mt_sequence, fit_options, pool_in_mt, t1obs, 1, scale_y); 
xfit(1) = R1W; 
disp(['xfit : ',get_str_x(xfit, par_labels)]);

%--------------------%
% forward and plot
%--------------------%
[~, ~, y_model, y_data_out] = obj_qMT (xfit, y_mt, mt_sequence, fit_options, pool_in_mt, t1obs);                    
freq = mt_sequence.facq(fit_options.m_vect); 

figure;
subplot(1,2,1); 
plot(freq, y_data_out, '.','MarkerSize',10); hold on; reset_color;
plot(freq, y_model, 'LineWidth',1); 
xlabel('Frequency offset (Hz)'); ylabel('MT signal');
subplot(1,2,2); 
ind_range = find(abs(freq) < 3000);
plot(freq(ind_range), y_data_out(ind_range,:), '.','MarkerSize',10); hold on; reset_color;
plot(freq(ind_range), y_model(ind_range,:), 'LineWidth',1);
xlabel('Frequency offset (Hz)'); ylabel('MT signal');
