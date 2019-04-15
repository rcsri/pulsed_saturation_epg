function [A] = bmsim_mtx(p0, f, w1, lstype)
% bmsim_mtx
%   build system matrix for bloch-mcconnell exchange equations
%
%   M' = A * M + C
%

% Specify 1 row per compound, set MT flag to 1
% For each compound, specify
%   df: resonance frequency (Hz)
%   R1: longitudinal relaxation rate (Hz)
%   R2: transverse relaxation rate (Hz)
%   M0:  concentration (M)
%   R:  exchange rate
%   MT: if semisolid pool
%

% exchange rate R is defined by first order steady state
%   k_ab  = R * M0b
%   k_ba  = R * M0a
%   k_ab = k_ba * f_b, where f_b = M0b / M0a

% 3 CEST, 1  MTC
%===================================================
%       df     R1      R2      M        R       MT
%===================================================
% p0 = [  dfa    R1a     R2a     M0a      0       0;
%         dfb    R1b     R2b     M0b      Rb      0;
%         dfc    R1c     R2c     M0c      Rc      0;
%         dfs    R1s     R2s     M0s      Rs      1;
%         ];

% validate input
% if ~is_valid_pool_params(p0)
%     disp('invalid pool parameters');
% end

% number of exchanging compounds
num_compounds = size(p0,1);
num_MT = p0(end,6);

% number of magnetization components
num_components = 3 * num_compounds - (num_MT * 2);
A = zeros(num_components, num_components);

M0a = p0(1,4);

A_relax     = 0 * A;
A_exchange  = 0 * A;
A_rf        = 0 * A;

for ixc = 1:num_compounds
    is_compound_MT = (p0(ixc,6) == 1);
    dfx = p0(ixc,1);
    R1x = p0(ixc,2);
    R2x = p0(ixc,3);
    M0x = p0(ixc,4);
    R   = p0(ixc,5);
    
    %=============================================
    % RELAXATION AND OFF-RESONANCE
    %=============================================
    if (~is_compound_MT)
        ix_mtx = (1:3) + 3*(ixc-1);
        A_relax(ix_mtx,ix_mtx) = ...
            [ -R2x               2*pi*(f-dfx)      0   ;
            -2*pi*(f-dfx)      -R2x              0   ;
            0                  0               -R1x ];
    else
        ix_mtx = 3*(ixc-1) + 1;
        
        if strcmp(lstype,'SL')
            g = lineshape_superlorentzian(f, 1./R2x);
        elseif strcmp(lstype,'GS')
            g = lineshape_gaussian(f, 1./R2x);
        elseif strcmp(lstype,'LZ')
            g = lineshape_lorentzian(f, 1./R2x);
        end
        
        Rrfb = pi * w1^2 * g;
        
        A_relax(ix_mtx,ix_mtx) = -R1x - Rrfb;
    end
    
    %=============================================
    % EXCHANGE TERMS
    %=============================================
    kxa = R * M0a;
    kax = R * M0x;
    
    if (~is_compound_MT)
        K_ax = [kax 0   0
            0   kax 0
            0   0   kax];
        K_xa = [kxa 0   0
            0   kxa 0
            0   0   kxa];
        
        ix_mtx_a = 1:3;
        ix_mtx_x = (1:3) + 3*(ixc-1);
        
        A_exchange(ix_mtx_a,ix_mtx_a) = A_exchange(ix_mtx_a,ix_mtx_a) - K_ax;
        A_exchange(ix_mtx_a,ix_mtx_x) = A_exchange(ix_mtx_a,ix_mtx_x) + K_xa;
        A_exchange(ix_mtx_x,ix_mtx_a) = A_exchange(ix_mtx_x,ix_mtx_a) + K_ax;
        A_exchange(ix_mtx_x,ix_mtx_x) = A_exchange(ix_mtx_x,ix_mtx_x) - K_xa;
    else
        ix_mtx = [3, 3*(ixc-1) + 1];
        
        % Note: Rrfb has already been taken into account above (line 80)
        A_exchange(ix_mtx,ix_mtx) = A_exchange(ix_mtx,ix_mtx) + ...
            [ -kax     kxa
            kax      -kxa ];
        
    end
    
    %=============================================
    % RF EXCITATION
    %=============================================
    if (~is_compound_MT)
        ix_mtx = (1:3) + 3*(ixc-1);
        A_rf(ix_mtx,ix_mtx) = ...
            [ 0 0 0
            0 0 w1
            0 -w1 0 ];
    end
end

A = A_relax + A_exchange + A_rf;
