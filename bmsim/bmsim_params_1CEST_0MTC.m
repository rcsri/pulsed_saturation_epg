function p = bmsim_params_1CEST_0MTC( ...
    dfa, R1a, R2a, M0a, dfb, R1b, R2b, M0b, Rb)
% Build parameter matrix for BM exchange simulation.
%
% 0CEST_1MTC:
%
% p = bmsim_params_1CEST_0MTC( ...
%     dfa, R1a, R2a, M0a, dfb, R1b, R2b, M0b, Rb);

%===================================================
%      df     R1      R2      M        R       MT
%===================================================
p = [  dfa    R1a     R2a     M0a      0       0; 
       dfb    R1b     R2b     M0b      Rb      0;
       ];