function [num_components, num_compounds] = bmsim_num_components(p)
% bmsim_num_components
%
%   extract number of compounds and components from pool parameter matrix

    % number of exchanging compounds
    num_compounds = bmsim_num_compounds(p);
    num_MT = p(end,6);

    % number of magnetization components
    num_components = 3 * num_compounds - (num_MT * 2);

end