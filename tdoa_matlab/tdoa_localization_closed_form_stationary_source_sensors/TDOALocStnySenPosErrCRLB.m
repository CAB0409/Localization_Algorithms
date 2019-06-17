function CRLB = TDOALocStnySenPosErrCRLB(so,uo,Q_beta,Q_alpha)
%
% This program finds the CRLB of the source position when the sensor
% positions are subject to random noises. The positioning parameter is
% TDOA. The derivation of the CRLB can be found in K. C. Ho, Xiaoning Lu, 
% L. Kovavisaruch, "Source Localization Using TDOA and FDOA Measurements 
% in the Presence of Receiver Location Errors: Analysis and Solution," IEEE
% Trans. on Signal Processing, Vol. 55, pp. 684-696, 2007.
%
% Usage: CRLB = TDOALocStnySenPosErrCRLB(so,uo,Q_beta,Q_alpha)
%
% Input Parameter List
% so: 2xM or 3xM true sensor position matrix.
%     M: number of sensors.
%     so(:,i) is the true position of the ith sensor.
% uo: 2x1 or 3x1 true source position.
% Q_beta: 2Mx2M, or 3Mx3M covariance matrix of the sensor position vector.
% Q_alpha: (M-1)x(M-1) covariance matrix of the TDOAs (range differences) from uo.
%
% The program returns the 2x2 or 3x3 CRLB of the source position.
%
%
% The program can be used for 2D(Dim=2) or 3D(Dim=3) localization.
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

M = size(so,2);       % Number of sensors.
N = size(so,1);       % Dimension of the localization problem.

for j = 2 : M         % Generate the partial derivatives.
    dr_du(j-1,:) =  (uo-so(:,j))'/norm(uo-so(:,j)) - (uo-so(:,1))'/norm(uo-so(:,1));
    dr_ds(j-1,:) = [(uo-so(:,1))'/norm(uo-so(:,1)), zeros(1,(j-2)*N), -(uo-so(:,j))'/norm(uo-so(:,j)), zeros(1,(M-j)*N)];
end;

X = dr_du'*inv(Q_alpha)*dr_du;
Y = dr_du'*inv(Q_alpha)*dr_ds;
Z = inv(Q_beta) + dr_ds'*inv(Q_alpha)*dr_ds;

CRLB = inv([X,Y;Y',Z]); 
CRLB = CRLB(1:N,1:N); % Return the CRLB of the source position.