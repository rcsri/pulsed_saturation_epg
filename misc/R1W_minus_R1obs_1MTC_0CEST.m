function difference = R1W_minus_R1obs_1MTC_0CEST(R1W,RMTC,M0MTC,R1MTC,R1obs)
%   MzW             MzMTC      
Az = [
    -R1W-RMTC*M0MTC RMTC;        % dMzW/dt
    RMTC*M0MTC      -R1MTC-RMTC; % dMzMTC/dt
    ];
[~,eigval] = eig(Az);

% Not sure why negative sign is necessary, but it was necessary to
% reproduce Eq. [21] in Henkelman. Magn Reson Med (1993) 29(6): 759-766.
eigval = -eigval;

% for i = 1:5
%     disp(eigval(i,i) - R1obs)
% end
% The minimum uses eigenvalue 2.

difference = min(abs(diag(eigval) - R1obs));
