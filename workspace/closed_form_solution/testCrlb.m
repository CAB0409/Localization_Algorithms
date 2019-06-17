
clear all;
load('loc_vars.mat');

SensorPositions = S;
SourceLocation = uo;

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