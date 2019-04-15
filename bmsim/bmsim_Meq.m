function [Meq, C] = bmsim_Meq(p0)
% Extract equilibrium magnetization and inhomogeneous term

% number of exchanging compounds
% num_compounds = bmsim_num_compounds(p0);

% number of magnetization components
[num_components, num_compounds] = bmsim_num_components(p0);

Meq = zeros(num_components, 1);
C   = zeros(num_components, 1);

% Build equilibrium magnetization and C vector
for ixc = 1:num_compounds
     is_compound_MT = (p0(ixc,6) == 1);
     
     R1x = p0(ixc,2);
     M0x = p0(ixc,4);
     
     if (~is_compound_MT)    
         ix_mtx = 3 + 3*(ixc-1);
     else
         ix_mtx = 3*(ixc-1) + 1;
     end
     
     C(ix_mtx) = R1x .* M0x;
     Meq(ix_mtx) = M0x;
end
end