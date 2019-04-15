function [expm_At, inv_A] = bmsim_compute_prop_shaped_pulse(p, f, w1, dt, lstype)
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

expm_At = zeros(num_components,num_components,num_pulses);
inv_A   = zeros(num_components,num_components,num_pulses);
for ix = 1:num_pulses
    
    A = bmsim_mtx( p, f, w1(ix) , lstype);
    [expm_At(:,:,ix), inv_A(:,:,ix)] = bmsim_compute_prop(dt(ix), A);
    
end
