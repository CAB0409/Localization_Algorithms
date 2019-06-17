
clear all;

load('loc_vars.mat');

r = rd;

RptCnt = 9;     % number of repetitions in Stage-1 to recompute W1

M = length(r) + 1;

%S=SensorPositions;
R=sqrt(sum(S.^2))';

%=========== construct related vector and matrix ============
h1 = r.^2 - R(2:end).^2 + R(1)^2;
G1 = -2*[S(:,2:end)'-ones(M-1,1)*S(:,1)' ,  r];

%============= first stage ===================================  
B = eye(M-1);
W1 = inv(B*Q*B');
pre_inv = B*Q*B';
pre_B = B;
Q_1 = Q;
W1_1 = W1;
u1 = inv(G1'*W1*G1)*G1'*W1*h1;
u1_1 = u1;

for j = 1:max(1,RptCnt),
    ri_hat = sqrt(sum((S-u1(1:end-1)*ones(1,M)).^2));
    B = 2*diag(ri_hat(2:M));  
    W1 = inv(B*Q*B');
    u1 = inv(G1'*W1*G1)*G1'*W1*h1;
    if j == 1
       ri_hat_1 = ri_hat;
       B_1 = B;
       W1_2 = W1;
       u1_2 = u1;
    end
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
