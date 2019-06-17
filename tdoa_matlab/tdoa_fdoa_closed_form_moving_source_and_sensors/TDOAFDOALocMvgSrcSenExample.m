% This program provides an example on how to call TDOAFDOALocMvgSrcSen.m 
% and TDOAFDOALocMvgSrcSenCRLB.m. 
%
% This program will reproduce the Fig. 6 in K. C. Ho and Wenwei Xu, 
% "An accurate algebraic solution for moving source location using TDOA and
% FDOA measurements," IEEE Trans. on Signal Processing, Vol. 52, 
% pp. 2453-2463, 2004.
% To generate Fig. 7 in the paper, the source location should be 
% uo = [600 650 550]' and not uo = [300 325 275]'.  The paper specifies a 
% wrong uo for Fig. 7.
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

clc; clear all; warning off;          % Program initialization.
L = 5e3;                              % Number of ensemble runs.            
randn('seed',0);                      % Initialize random number generator.

uo = [2000 2500 3000]';               % True source position.
u_doto = [-20 15 40]';                % True source velocity.

x = [300,400,300,350,-100];
y = [100 150 500 200 -100];
z = [150 100 200 100 -100];
s = [x; y; z];                        % True sensor position matrix.

M = size(s,2);                        % Number of sensors.
N = size(s,1);                        % Dimension of localization problem.

xd = [ 30 -30  10 10 -20];
yd = [-20  10 -20 20  10];
zd = [ 20  20  10 30  10];
s_dot = [xd; yd; zd];                 % True sensor velocity.

ro = sqrt(sum((repmat(uo,1,M-1)-s(:,2:M)).^2,1))'; %distance between source and receiver
r_doto=(repmat(uo,1,M-1)-s(:,2:M))'*(repmat(u_doto,1,M-1)-s_dot(:,2:M));

r_doto = diag(r_doto)./ro;
% True TDOAs (range differences).
ro = ro - norm(uo - s(:,1));                           
% True FDOAs (range rate differences).
r_doto = r_doto - (uo-s(:,1))'*(u_doto-s_dot(:,1))/norm(uo-s(:,1)); 

R = (eye(M-1)+ones(M-1))/2;           % Covariance matrix of measurements.
R = [R,zeros(M-1); zeros(M-1),0.1*R]; 
i = 1;                                % Program counter.
fprintf('Simulation in progress');

for sigma2_d = -30 : 2 : 0            % Vary measurement noise power.
    fprintf('.'); 
    
    Q_alpha = 10^(sigma2_d/10) * R;
    CRLB = TDOAFDOALocMvgSrcSenCRLB(s,s_dot,uo,u_doto,Q_alpha);
    
    crlbu(i) = trace(CRLB(1:N,1:N));
    crlbudot(i) = trace(CRLB(N+1:end,N+1:end));
    
    SimulationMSEu = 0; SimulationMSEudot = 0;
    chol_Q_alpha=chol(Q_alpha);
    for k = 1 : L                     % Monte Carlo simulation.
        % Noisy TDOA and FDOA measurements.
        Delta_r = chol_Q_alpha'*randn(2*(M-1),1);
        r = ro + Delta_r(1:M-1);
        r_dot = r_doto + Delta_r(M:end);
       
        [u,u_dot] = TDOAFDOALocMvgSrcSen(s,s_dot,r,r_dot,Q_alpha);
        SimulationMSEu = SimulationMSEu + norm(u-uo)^2;
        SimulationMSEudot = SimulationMSEudot + norm(u_dot-u_doto)^2;
    end;
    mseu(i) = SimulationMSEu/L;
    mseudot(i) = SimulationMSEudot/L;
    i = i + 1;                        % Update program counter.    
end;
fprintf('\n');

% Plot the result.
figure(1); subplot(211);
semilogy(-30:2:0,sqrt(mseu),'kx','MarkerSize',8); grid on; hold on;
semilogy(-30:2:0,sqrt(crlbu),'k','LineWidth',1); 
xlabel('10log(c^2\sigma_d^2(m^2))'); ylabel('Position RMSE(m)'); 
legend('Simulation RMSE, proposed solution','CRLB',2); %#ok<LEGINTPAR>
axis([-30,0,10,1000]);
hold off;

subplot(212);
semilogy(-30:2:0,sqrt(mseudot),'kx','MarkerSize',8); grid on; hold on;
semilogy(-30:2:0,sqrt(crlbudot),'k','LineWidth',1); 
xlabel('10log(c^2\sigma_d^2(m^2))'); ylabel('Velocity RMSE(m/s)'); 
legend('Simulation RMSE, proposed solution','CRLB',2); %#ok<LEGINTPAR>
axis([-30,0,10,1000]);
hold off;
