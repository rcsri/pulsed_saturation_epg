function p = bmsim_params_0CEST_0MTC( ...
    dfa, R1a, R2a, M0a)
% Build parameter matrix for BM exchange simulation.
%
% 0CEST_0MTC:
%
% p = bmsim_params_0CEST_0MTC( ...
%     dfa, R1a, R2a, M0a);

%===================================================
%      df     R1      R2      M        R       MT
%===================================================
p = [  dfa    R1a     R2a     M0a      0       0; 
       ];