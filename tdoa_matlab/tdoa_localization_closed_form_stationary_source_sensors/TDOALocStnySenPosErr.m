function u = TDOALocStnySenPosErr(s,r,Q_beta,Q_alpha)
%
% This program realizes the algorithm for localizing a stationary source 
% using TDOAs in the presence of sensor position errors. The sensors are 
% stationary as well. The details of the algorithm development can be found
% in K. C. Ho, Xiaoning Lu, L. Kovavisaruch, "Source Localization Using 
% TDOA and FDOA Measurements in the Presence of Receiver Location Errors: 
% Analysis and Solution," IEEE Trans. on Signal Processing, Vol. 55, 
% pp. 684-696, 2007.
%
% Usage: u = TDOALocStnySenPosErr(s,r,Q_beta,Q_alpha).
% 
% Input Parameter List:
% s : 2xM or 3xM sensor position matrix.
%     M: number of sensors.
%     s(:,i) is the known position of the ith sensor.
% r =[r_21,r_31,...r_M1]': (M-1)x1 TDOA (range difference) measurement vector.
% Q_beta: 2Mx2M, or 3Mx3M covariance matrix of the sensor positions
%                                       reshape(s,2M,1) or reshape(s,3M,1).
% Q_alpha: (M-1)x(M-1) covariance matrix of r.
%
% The program returns a 2x1 or 3x1 source location estimate.
%
% Note: W1 is updated 3 times (RptCnt=3) in Stage-1, however in most
% cases updating W1 once (RptCnt=1) is sufficient.
%
% The program can be used for 2D(Dim=2) or 3D(Dim=3) localization.
%
%
% Le Yang, K. C. Ho      08-01-2009
%                        10-01-2010, revised
%
%       Copyright (C) 2009
%       Computational Intelligence Signal Processing Laboratory
%       University of Missouri
%       Columbia, MO 65211, USA.
%       hod@missouri.edu
%

RptCnt = 3;     % number of repetitions in Stage-1 to recompute W1

Qb = Q_beta; Qa = Q_alpha; 
M = size(s,2);                   % Number of Sensors.
N = size(s,1);                   % Dimension of the localization problem.

% Stage-1 processing.    
for j = 2 : M
    h1(j-1,1) = r(j-1)^2 - s(:,j)'*s(:,j) + s(:,1)'*s(:,1);
    G1(j-1,:) = -2 * [s(:,j)'-s(:,1)',r(j-1)];
end;
W1 = inv(Qa);
phi1 = inv(G1'*W1*G1)*G1'*W1*h1;
u = phi1(1:N); 

for iteration = 1 : max(1,RptCnt)            % Iterate to update W1.              
    for j = 2 : M
        b(j-1,1) = norm(u-s(:,j),2);
        D1(j-1,:)= [-(u-s(:,1))',zeros(1,(j-2)*N),(u-s(:,j))',zeros(1,(M-j)*N)];
    end;
    B1 = 2 * diag(b);
    D1 = 2 * D1;
    W1 = inv(B1*Qa*B1'+D1*Qb*D1');
    phi1 = inv(G1'*W1*G1)*G1'*W1*h1;
    u = phi1(1:N);
end;
cov_phi1 = inv(G1'*W1*G1);

% Stage-2 processing.
h2 = (phi1-[s(:,1);0]).^2;
G2 = [eye(N); ones(1,N)];
B2 = 2 * diag([u-s(:,1); norm(u-s(:,1))]);
D2 = [zeros(N,M*N);...
      (u - s(:,1))',zeros(1,(M-1)*N)];
D2 = 2* D2;
W2 = inv(B2*cov_phi1*B2'+D2*Qb*D2' + B2*cov_phi1*G1'*W1*D1*Qb*D2'...
                                   + D2*Qb*D1'*W1*G1*cov_phi1*B2');
phi2 = inv(G2'*W2*G2)*G2'*W2*h2;

% Return the source location estimate.
u = diag(sign(u-s(:,1)))*sqrt(phi2)+ s(:,1);