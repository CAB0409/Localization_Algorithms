% This program provides an example on how to call TDOALoc.m 
% and TDOALocCRLB.m. 
%
% This program will reproduce the new method result in Fig. 4 of Y. T. Chan 
% and K. C. Ho, "A Simple and Efficient Estimator for Hyperbolic Location,"
% IEEE Trans. Signal Processing, vol. 42, no. 8, pp. 1905-1915, Aug. 1994.
%
%
% Ming Sun, K. C. Ho     08-01-2009
%
%       Copyright (C) 2009
%       Computational Intelligence Signal Processing Laboratory
%       University of Missouri
%       Columbia, MO 65211, USA.
%       hod@missouri.edu
%

clc; close all; clear all;warning off;  % Program initialization.
L = 10e3;                               % Number of ensemble runs.            
randn('seed',0);                        % Initialize random number generator.

uo = [-50 250]';                        % True source position.

x=[0 -5 4 -2 7 -7 2 -4 3 1];            % True sensor position matrix.
y=[0  8 6  4 3  5 5  2 3 8];
S=[x; y];

M = size(S,2);                          % Number of sensors.
N = size(S,1);                          % Dimension of localization.

ro = sqrt(sum((uo*ones(1,M)-S).^2))';   % True source-sensor ranges
rdo = ro(2:end)-ro(1);

R = (eye(M-1)+ones(M-1))/2;             % covariance structure of TDOA

NsePwrVecdB=-60:4:-24; 

fprintf('Simulation in progress');
for nseIdx=1:length(NsePwrVecdB),       % vary measurement noise level
    fprintf('.');
    nsePwr = 10^(NsePwrVecdB(nseIdx)/10);
    Q = nsePwr * R;                     % Covariance matrix of TDOA (range difference) noise
    if nseIdx == 1
        Q_1 = Q;
    end
    %lowest crlb contains most likely location
    crlb(nseIdx)=trace(TDOALocCRLB(S,uo,Q));
    if nseIdx == 1
        nse_1 = nsePwr;
        crlb_1 = crlb(1);
        R_1 =R;
    end
    
    SimulationMSE = 0;
    for k = 1 : L,                      % Monte Carlo Simulation.
        rdNse = sqrt(nsePwr/2) * randn(M,1);
        rd = rdo + rdNse(2:end)-rdNse(1);   % Noisy source TDOAs (range differences).
   
        u = TDOALoc(S,rd,Q);
        if nseIdx == 1 && k == 1
            fprintf('\nsaving variables\n');
            save('loc_vars.mat', 'S', 'rd','Q', 'uo');
            S_1 =S;
            rd_1 = rd;
            Q_2 = Q;
            actual_location = u;
            sim_1  = norm(u-uo,2)^2;
            norm_1 = norm(u-uo,2)^2;
        end
        SimulationMSE = SimulationMSE + norm(u-uo,2)^2;
    end;

    mse(nseIdx) = SimulationMSE/L;
end;
fprintf('\n');
% Plot the results.
figure(1); plot(NsePwrVecdB/2,10*log10(mse),'xk','MarkerSize',8); hold on;
plot(NsePwrVecdB/2,10*log10(crlb),'k'); grid on; hold off;

xlabel('10 log(c\sigma)'); 
ylabel('10 log(MSE)');
legend('New Method','CRLB');
ylim([0 60]);

