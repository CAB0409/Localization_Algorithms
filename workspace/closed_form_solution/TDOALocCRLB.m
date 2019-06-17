function [CRLB] = TDOALocCRLB(SensorPositions,SourceLocation,Q)
% SensorPositons: (Dim x M) matrix, each column is a sensor position and 
%                 first column is the reference sensor
%                 M is the number of sensors and should be at least Dim+2
%                 the sensors should not lie in one plane or line
% SourceLocation: (Dim x 1) vector of the estimated source location
% Q:              the covariance matrix of the TDOA (range difference) vector
% CRLB:           (Dim x Dim) CRLB matrix of the estimated source localization
%
% The program can be used for 2D(Dim=2) or 3D(Dim=3) localization
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

M = length(Q) + 1;

if (M<size(SensorPositions,1)+2)
    fprintf('Number of sensors must be at least %d\n',size(SensorPositions,1)+2);
    return;
end;

if (rank(SensorPositions) < size(SensorPositions,1))
    disp('The sensors should not lie in one plane or line!');
    return;
end

S = SensorPositions; 
u = SourceLocation;

M = size(SensorPositions,2);
ro = sqrt(sum((u*ones(1,M)-S).^2));

d_u = (S(:,2:end)-u*ones(1,M-1))'./(ro(2:end)'*ones(1,size(S,1))) ...
      -ones(M-1,1)*((S(:,1)-u)'/ro(1));

J = d_u'*inv(Q)*d_u;    % FIM
CRLB = inv(J);

