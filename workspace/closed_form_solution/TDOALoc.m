function [SourceLocation] = TDOALoc(SensorPositions, r, Q),
% SensorPositons: (Dim x M) matrix, each column is a sensor position and 
%                 first column is the reference sensor
%                 the sensors should not lie in one plane or line
% r:           a (M-1) x 1 vector of TDOA measurements times signal propagation speed
%                 M is the number of sensors and should be at least Dim+2
% Q:              the covariance matrix of the r vector
% SourceLocation: estimated source location
%
% Note: W1 is updated 3 times (RptCnt=3) in Stage-1, however in most
% cases updating W1 once (RptCnt=1) is sufficient.
%
% The program can be used for 2D(Dim=2) or 3D(Dim=3) localization
%
%
% Ming Sun, K. C. Ho     08-01-2009
%                        10-01-2010, revised
%
%       Copyright (C) 2009
%       Computational Intelligence Signal Processing Laboratory
%       University of Missouri
%       Columbia, MO 65211, USA.
%       hod@missouri.edu
%

RptCnt = 9;     % number of repetitions in Stage-1 to recompute W1

M = length(r) + 1;

if (M<size(SensorPositions,1)+2)
    fprintf('Number of sensors must be at least %d\n',size(SensorPositions,1)+2);
    SourceLocation=NaN;
    return;
end;

if (rank(SensorPositions) < size(SensorPositions,1))
    disp('The sensors should not lie in one plane or line!');
    SourceLocation=NaN;
    return;
end

S=SensorPositions;
R=sqrt(sum(S.^2))';

%=========== construct related vector and matrix ============
h1 = r.^2 - R(2:end).^2 + R(1)^2;
G1 = -2*[S(:,2:end)'-ones(M-1,1)*S(:,1)' ,  r];

%============= first stage ===================================  
B = eye(M-1);
W1 = inv(B*Q*B');
u1 = inv(G1'*W1*G1)*G1'*W1*h1;

for j = 1:max(1,RptCnt),
    ri_hat = sqrt(sum((S-u1(1:end-1)*ones(1,M)).^2));
    B = 2*diag(ri_hat(2:M));  
    W1 = inv(B*Q*B');
    u1 = inv(G1'*W1*G1)*G1'*W1*h1;
end

u1p = u1 - [S(:,1);0];

%========== second stage =====================================
h2 = u1p.^2;
G2 = [eye(length(u1p)-1);ones(1,length(u1p)-1)];
    
B2 = 2*diag(u1p);
W2 = inv(B2')*(G1'*W1*G1)*inv(B2);
u2 = inv(G2'*W2*G2)*G2'*W2*h2;

%=========== mapping ========================================
SourceLocation = sign(diag(u1p(1:length(u2))))*sqrt(abs(u2)) + S(:,1);
%============================================================

if u1(end) < 0 || min(u2) < 0
    SourceLocation = u1(1:length(u2));
end