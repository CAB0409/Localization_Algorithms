function CRLB = TDOAFDOALocMvgSrcSenCRLB(s,s_dot,uo,u_doto,Q_alpha)
%
% This program finds the CRLB of the source position and velocity. The 
% positioning parameters are TDOA and FDOA. The derivation of the CRLB can 
% be found in K. C. Ho and Wenwei Xu, "An accurate algebraic solution for 
% moving source location using TDOA and FDOA measurements," IEEE Trans. on 
% Signal Processing, Vol. 52, pp. 2453-2463, 2004.
%
% Usage: CRLB = TDOAFDOALocMvgSrcSenCRLB(s,s_dot,uo,u_doto,Q_alpha)
%
% Input Parameter List
% s : 2xM or 3xM sensor position matrix.
%     M: number of sensors.
%     s(:,i) is the position of the ith sensor.
% s_dot: 2xM or 3xM sensor velocity matrix.
%        s_dot(:,i) is the velocity of the ith sensor.
% uo: 2x1 or 3x1 true source position. 
% u_doto: 2x1 or 3x1 true source velocity.
% Q_alpha: 2(M-1)x2(M-1) covariance matrix of the TDOAs (range differences)
%                           and FDOAs (range rate differences).
%
% The program returns the 4x4 or 6x6 CRLB of [uo;u_doto].
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

M = size(s,2);           % Number of sensors.
N = size(s,1);           % Dimension of the localization problem.

for j = 2 : M            % Generate the partial derivatives.
    dr_du(j-1,:) = (uo-s(:,j))'/norm(uo-s(:,j))-(uo-s(:,1))'/norm(uo-s(:,1));
    dr_dudot(j-1,:) = zeros(1,N);
    
    drdot_du(j-1,:) = (u_doto-s_dot(:,j))'/norm(uo-s(:,j)) - (u_doto-s_dot(:,j))'*(uo-s(:,j))*(uo-s(:,j))'/norm(uo-s(:,j))^3; 
    drdot_du(j-1,:) = drdot_du(j-1,:) - (u_doto-s_dot(:,1))'/norm(uo-s(:,1)) + (u_doto-s_dot(:,1))'*(uo-s(:,1))*(uo-s(:,1))'/norm(uo-s(:,1))^3; 
    drdot_dudot(j-1,:) = dr_du(j-1,:);
end;

dq_dtheta = [dr_du,   dr_dudot;
             drdot_du,drdot_dudot];

X = dq_dtheta'*inv(Q_alpha)*dq_dtheta;
CRLB = inv(X);           % Return the CRLB.