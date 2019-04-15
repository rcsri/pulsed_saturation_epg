function M = bmsim_propagate_ver2( M0, expm_At, B )
%function M = bmsim_propagate( M0, expm_At, inv_A, C )

% Propagate solution for coupled differential equation
% Solution to dM/dt = A*M + C;
%
% M = exp(A*t)*M(t=0) + A^(-1) * (exp(A*t) - I) * C
%
% inputs:
%   M0:         initial condition
%   expm_At:    exp(A*t)
%   inv_A:      inv(A)
%   C:          vector
%
%   expm_At and inv_A are [n,n,num_intervals]

num_intervals = size(expm_At,3);

M = M0;

for ix = 1:num_intervals
    E = expm_At(:,:,ix);
    M = E*M + B(:,ix);
end

