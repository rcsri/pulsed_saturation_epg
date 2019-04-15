function MzA = bmsim_MzA(M)
% bmsim_MzA
%   extract MzA component from magnetization vector

MzA = M(3,:,:);
MzA = permute(MzA, [3 2 1]);