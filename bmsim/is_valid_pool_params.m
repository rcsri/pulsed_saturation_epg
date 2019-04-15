function valid = is_valid_pool_params(params)
% return TRUE if valid

% Specify 1 row per compound, set MT flag to 1
% For each compound, specify
%   df: resonance frequency (Hz)
%   R1: longitudinal relaxation rate (Hz)
%   R2: transverse relaxation rate (Hz)
%   M:  concentration
%   R:  exchange rate
%   MT: if semisolid pool
%       df     R1      R2      M       R MT
%
% example:
%
% params = [  dfa    R1a     R2a     Ma   0   0  ; 
%             dfb    R1b     R2b     Mb   Rb  0  ; 
%             dfc    R1c     R2c     Mc   Rc  0  ;
%             dfs    R1s     R2s     Ms   Rd  1  ;
%          ];

valid = true;

if ( any(params(:,2) <= 0) || ...
     any(params(:,3) <= 0) || ...
     any(params(:,4) < 0)  || ...
     params(1,5) ~= 0 || ...
     any(params(1:end-1,6) == 1) )
 
    valid = false;
    
end
