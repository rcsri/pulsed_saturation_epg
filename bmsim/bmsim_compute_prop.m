function [expm_At, inv_A] = bmsim_compute_prop(t,A)
% oompute propagator for matrix differential equation
%
% M' = A*M + C

expm_At = expm(A*t);
inv_A  = inv(A);