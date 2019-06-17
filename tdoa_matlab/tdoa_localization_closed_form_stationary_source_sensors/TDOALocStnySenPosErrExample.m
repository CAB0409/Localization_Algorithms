% This program provides an example on how to call TDOALocStnySenPosErr.m 
% and TDOALocStnySenPosErrCRLB.m. 
%
% This program will reproduce the Fig. 10 in K. C. Ho, Xiaoning Lu, 
% L. Kovavisaruch, "Source Localization Using TDOA and FDOA Measurements 
% in the Presence of Receiver Location Errors: Analysis and Solution," IEEE
% Trans. on Signal Processing, Vol. 55, pp. 684-696, 2007.
% Please note that the paper has a typo on page 686 near the end of column
% 2 in specifying the sensor position noise power.  Rs should be
% Rs=sigma_s^2[1,1,1,2,2,2,10,10,10,40,40,40,20,20,20,3,3,3] and it is
% corrected in this code.
%
%
% Le Yang, K. C. Ho      08-01-2009
%                        09-05-2010, revised
%
%       Copyright (C) 2009
%       Computational Intelligence Signal Processing Laboratory
%       University of Missouri
%       Columbia, MO 65211, USA.
%       hod@missouri.edu
%

clc; clear all; warning off;	      % Program initialization.
L = 5e3;                              % Number of ensemble runs.            
randn('seed',0);                      % Initialize random number generator.

uo = [2000 2500 3000]';               % True source position.

x = [300,400,300,350,-100,200];       % True sensor position matrix.
y = [100 150 500 200 -100,-300];
z = [150 100 200 100 -100,-200];
so = [x; y; z];  

M = size(so,2);                       % Number of sensors.
N = size(so,1);                       % Dimension of localization.

ro = sqrt(sum((repmat(uo,1,M-1)-so(:,2:M)).^2,1))';
ro = ro - norm(uo-so(:,1));           % True source TDOAs (range differences).

R = (eye(M-1)+ones(M-1))/2;
J = diag([1,1,1,2,2,2,10,10,10,40,40,40,20,20,20,3,3,3]);
Q_alpha = 1e-4 * R;                   % Covariance matrix of source TDOA (range difference) noise. 
i = 1;                                % Program counter.
fprintf('Simulation in progress');

for sigma2_s = -60 : 5 : 0            % Vary sensor position noise power.
    fprintf('.');
    Q_beta = J * 10^(sigma2_s/10);    % Covariance matrix of sensor position noise.
    crlb(i) = trace(TDOALocStnySenPosErrCRLB(so,uo,Q_beta,Q_alpha));
    
    SimulationMSE = 0;
    chol_Q_alpha=chol(Q_alpha);
    chol_Q_beta=chol(Q_beta);
    
    for k = 1 : L                     % Monte Carlo Simulation.
        delta_s = chol_Q_beta'*randn(M*N,1);
        s = so+ reshape(delta_s,N,M); % Noisy sensor positions.
        
        delta_r = chol_Q_alpha'*randn(M-1,1);
        r = ro + delta_r;             % Noisy source TDOAs (range differences).
        u = TDOALocStnySenPosErr(s,r,Q_beta,Q_alpha);
        SimulationMSE = SimulationMSE + norm(u-uo,2)^2;
    end;

    mse(i) = SimulationMSE/L;
    i = i + 1;                        % Update program counter.
end;
fprintf('\n');

% Plot the results.
figure(1); plot(-60:5:0,10*log10(mse),'kx','MarkerSize',8); hold on;
plot(-60:5:0,10*log10(crlb),'k','LineWidth',1);
xlabel('10log(\sigma_s^2(m^2))'); grid on; box on;
ylabel('10log(position MSE(m^2))');
legend('Simulation MSE, proposed solution', ...
    'CRLB with receiver location error',2);       
hold off;