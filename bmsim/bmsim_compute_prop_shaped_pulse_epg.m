function [expm_UAUt, inv_UAU] = bmsim_compute_prop_shaped_pulse_epg(p, f, w1, dt, U, inv_U, lstype)
% Compute sequence of matrix propagators
%
%   [expm_At, inv_A] = bmsim_compute_prop_shaped_pulse(p, f, w1, dt);
%
%   p = matrix of pool parameters
%   f = irradiation frequency
%   w1 = shaped pulse amplitude (rad/s)
%   dt = shaped pulse time steps (s)
%
% This function is short hand for
%   A = bmsim_mtx(p,f,w1);            % generate system matrix
%   [E,I] = bmsim_compute_prop(dt,A); % propagator for time step dt

num_components = bmsim_num_components(p);
num_pulses = numel(w1);

expm_UAUt = zeros(num_components,num_components,num_pulses);
inv_UAU   = zeros(num_components,num_components,num_pulses);
for ix = 1:num_pulses
    
    A = bmsim_mtx( p, f, w1(ix), lstype );
    UAU = U * A * inv_U;
    [expm_UAUt(:,:,ix), inv_UAU(:,:,ix)] = bmsim_compute_prop(dt(ix), UAU);
    
end
