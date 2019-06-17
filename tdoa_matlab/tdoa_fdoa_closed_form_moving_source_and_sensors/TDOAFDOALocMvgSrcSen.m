function [u,u_dot] = TDOAFDOALocMvgSrcSen(s,s_dot,r,r_dot,Q_alpha)
%
% This program realizes the algorithm for localizing a moving source 
% using TDOAs and FDOAs. The sensors are moving as well. The details of the
% algorithm development can be found in K.C. Ho and Wenwei Xu, "An accurate
% algebraic solution for moving source location using TDOA and FDOA 
% measurements," IEEE Trans. on Signal Processing, Vol. 52, pp. 2453-2463, 
% 2004.
%
% Usage: [u,u_dot] = TDOAFDOALocMvgSrcSen(s,s_dot,r,r_dot,Q_alpha).
%
% Input Parameter List:
% s : 2xM or 3xM sensor position matrix.
%     M: number of sensors.
%     s(:,i) is the known position of the ith sensor.
% s_dot: 2xM or 3xM sensor velocity matrix.
%        s_dot(:,i) is the known velocity of the ith sensor.
% r = [r_21,r_31,...r_M1]': (M-1)x1 TDOA (range rate) measurement vector.
% r_dot = cx[f_21,f_31,...f_M1]': (M-1)x1 FDOA measurement vector.
%         c is the signal propagation speed.
% Q_alpha: 2(M-1)x2(M-1) covariance matrix of [r;r_dot].
%
% The program returns a 2x1 or 3x1 source location estimate u and a 2x1 or 
% 3x1 source velocity estimate u_dot.
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

RptCnt = 3;              % number of repetitions in Stage-1 to recompute W1

Qa = Q_alpha; 
M = size(s,2);           % Number of sensors.
N = size(s,1);           % Dimension of the localization problem.

% Stage-1 processing.
for j = 2 : M
    ht(j-1,1) = r(j-1)^2 - s(:,j)'*s(:,j) + s(:,1)'*s(:,1);
    hf(j-1,1) = 2*(r(j-1)*r_dot(j-1)-s_dot(:,j)'*s(:,j) + s_dot(:,1)'*s(:,1));

    Gt(j-1,:) = -2 * [(s(:,j)-s(:,1))',  r(j-1),zeros(1,N+1)];
    Gf(j-1,:) = -2 * [(s_dot(:,j)-s_dot(:,1))',r_dot(j-1),(s(:,j)-s(:,1))',r(j-1)];
end;

h1 = [ht; hf];
G1 = [Gt; Gf];
W1 = inv(Qa);
phi1 = inv(G1'*W1*G1)*G1'*W1*h1;

u = phi1(1:N);
u_dot= phi1(N+2:end-1);

for m = 1 : max(1,RptCnt)            % Iterate to update W1.
    for j = 2 : M
        b(j-1,1) = norm(u-s(:,j));
        b_dot(j-1,1)= (u-s(:,j))'*(u_dot-s_dot(:,j))/norm(u-s(:,j));
    end;
    B = 2 * diag(b);
    B_dot = 2 * diag(b_dot);
    B1 = [B,zeros(size(B));B_dot,B];
    W1 = inv(B1*Qa*B1');

    phi1 = inv(G1'*W1*G1)*G1'*W1*h1;
    u = phi1(1:N);
    u_dot = phi1(N+2:end-1);
end;
cov_phi1 = inv(G1'*W1*G1);

% Stage-2 processing.
h2 = [(u-s(:,1)).^2;                  phi1(N+1)^2;...
      (u-s(:,1)).*(u_dot-s_dot(:,1)); phi1(N+1)*phi1(end)];
G2 = [eye(N),    zeros(N);...
      ones(1,N), zeros(1,N);...
      zeros(N),  eye(N);...
      zeros(1,N),ones(1,N)];
for m = 1 : 2            % Iterate once to update W2.
    Bt =  diag([u-s(:,1);    norm(u-s(:,1))]);
    Btd = diag([u_dot-s_dot(:,1);(u-s(:,1))'*(u_dot-s_dot(:,1))/norm(u-s(:,1))]);
    B2 = [2*Bt,zeros(N+1); Btd, Bt];
    W2 = inv(B2*cov_phi1*B2');
    phi2 = inv(G2'*W2*G2)*G2'*W2*h2;

    u = diag(sign(u-s(:,1)))*sqrt(abs(phi2(1:N))) + s(:,1);
    u_dot = phi2(N+1:end)./(u-s(:,1)) + s_dot(:,1);
end;           