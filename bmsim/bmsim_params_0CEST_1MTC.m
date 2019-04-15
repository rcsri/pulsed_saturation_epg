function p = bmsim_params_0CEST_1MTC( ...
    dfa, R1a, R2a, M0a, dfs, R1s, R2s, M0s, Rs)
% Build parameter matrix for BM exchange simulation.
%
% 0CEST_1MTC:
%
% p = bmsim_params_0CEST_1MTC( ...
%     dfa, R1a, R2a, M0a, dfs, R1s, R2s, M0s, Rs);

%===================================================
%      df     R1      R2      M        R       MT
%===================================================
p = [  dfa    R1a     R2a     M0a      0       0;   % water pool  
       dfs    R1s     R2s     M0s      Rs      1;   % MT pool 
       ];