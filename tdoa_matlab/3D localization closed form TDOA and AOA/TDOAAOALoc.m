% This program demonstrates the performance of the proposed method for
% source localization using hybrid TDOA and AOA measurements
%
% Reference: J. Yin, Q. Wan, S. Yang, and K. C. Ho, �A simple and accurate 
% TDOA-AOA localization method using two stations,� IEEE Signal Process. 
% Lett., vol. 23, pp. 144-148, Jan. 2016.
%
% This code generates Fig. 2 in the reference by setting VryNseTDOA equal 
% 1 and Fig. 3 by setting it to 0.
%
% Jihao Yin, K. C. Ho    06-28-2016
%
%       Copyright (C) 2016, Jihao Yin, Dominic Ho
%       All rights reserved
%
%       hod@missouri.edu

%clc
clear;

VryNseTDOA=0;               % performance vs tdoa noise 
VryNseAOA=1-VryNseTDOA;     % performance vs aoa noise 

% the source's actual position, unit:meter.
u=[1000;1000;1000];

% the number of station placements.
L1=10;

% the radius of the placement circle, unit:meter.
cr=300;

% TDOA error, unit:meter square.
sigmaTDOA=10^2;
tdoanf=[10^0.5,10,10^1.5,100,10^2.5,1000].^2;

% AOA error, unit:rad square.
sigmaAOA=(pi/360)^2;
aoanf=([0.5:0.5:2.5]/180*pi).^2;

% performance vs sigma_TDOA or vs sigma_AOA
if (VryNseTDOA) nf=tdoanf; else nf=aoanf; end;

% -- test results. --

% results of the proposed method.
eresult=zeros(length(nf),L1);
teresult=zeros(length(nf),1);

% results of Taylor series method.
eresult1=zeros(length(nf),L1);
teresult1=zeros(length(nf),1);

% theoretical mse results.
mseresult=zeros(length(nf),L1);
tmseresult=zeros(length(nf),1);

% CRB results.
CRBresult=zeros(length(nf),L1);
tCRBresult=zeros(length(nf),1);

% estimation bias square results.
sebias=zeros(length(nf),L1);
tsebias=zeros(length(nf),1);

% for the solution of Taylor series method by good initial guess.
sebias1=zeros(length(nf),L1);
tsebias1=zeros(length(nf),1);

% the number of the ensemble of Monte-Carlo
L=5000;

for il=1:L1,
    
    fprintf('.');
    
    % produce a random number in [0,2pi] with uniform distribution.
    sphai=unifrnd(0,2*pi);
    
    % the stations' positions (each row, each station), unit:meter.
    s0=cr*[cos(sphai),sin(sphai),0];
    
    S=[s0;-s0];
    
    [srow scol]=size(S);
    
    % the actual AOAs ( the first row is theta, the second row is phi), unit:rad.
    AOA=zeros(2,srow);
    
    for irow=1:srow,
        AOA(1,irow)=atan2(u(2)-S(irow,2),u(1)-S(irow,1));
        AOA(2,irow)=atan2(u(3)-S(irow,3),sqrt((u(1)-S(irow,1))^2+(u(2)-S(irow,2))^2));
    end
    
    % the actual TDOA, unit:meter.
    TDOA=zeros(srow-1,1);
    
    r1=norm(u-S(1,:).',2);
    
    for irow=1:srow-1,
        TDOA(irow)=norm(u-S(irow+1,:).',2)-r1;
    end
    
    % background noise preprocessing
    tdoanoise=randn(L,srow-1);
    for icol=1:srow-1,
        tdoanoise(:,icol)=tdoanoise(:,icol)-sum(tdoanoise(:,icol))/L;
    end
    
    
    aoanoise1=randn(L,srow);
    aoanoise2=randn(L,srow);
    for icol=1:srow,
        aoanoise1(:,icol)=aoanoise1(:,icol)-sum(aoanoise1(:,icol))/L;
        aoanoise2(:,icol)=aoanoise2(:,icol)-sum(aoanoise2(:,icol))/L;
    end
    
    for nfnum=1:length(nf),
        
        if (VryNseTDOA)
            % variance of TDOA error, unit:meter square.
            sigmaTDOA=tdoanf(nfnum);
        else
            % variance of AOA error, unit:rad square.
            sigmaAOA=aoanf(nfnum);
        end;
        
        % covariance matrix.
        Q=zeros(3*srow-1,3*srow-1);
        
        Q(srow:end,srow:end)=sigmaAOA*eye(2*srow);
        
        Q(1:srow-1,1:srow-1)=sigmaTDOA*(eye(srow-1)+ones(srow-1))/2;
        
        % the theoretical mse.
        
        % G here is G^T in the paper.
        G=zeros(3*srow-1,3);
        
        rowindex1=1;
        rowindex2=1;
        
        for irow=1:srow-1,
            
            G(rowindex1,:)=2*[cos(AOA(2,irow+1))*cos(AOA(1,irow+1))-cos(AOA(2,1))*cos(AOA(1,1)),cos(AOA(2,irow+1))*sin(AOA(1,irow+1))-cos(AOA(2,1))*sin(AOA(1,1)),sin(AOA(2,irow+1))-sin(AOA(2,1))];
            
            rowindex1=rowindex1+1;
        end
        
        for irow=1:srow,
            
            G(srow-1+rowindex2,:)=[sin(AOA(1,irow)),-cos(AOA(1,irow)),0];
            G(srow-1+rowindex2+1,:)=[sin(AOA(2,irow))*cos(AOA(1,irow)),sin(AOA(2,irow))*sin(AOA(1,irow)),-cos(AOA(2,irow))];
            
            rowindex2=rowindex2+2;
        end
        
        B=zeros(3,srow);
        
        for irow=1:srow,
            B(:,irow)=[cos(AOA(2,irow))*cos(AOA(1,irow));cos(AOA(2,irow))*sin(AOA(1,irow));sin(AOA(2,irow))];
        end
        
        Lm=zeros(3,2,srow);
        
        for irow=1:srow,
            Lm(:,:,irow)=[-cos(AOA(2,irow))*sin(AOA(1,irow)),-sin(AOA(2,irow))*cos(AOA(1,irow)); ...
                cos(AOA(2,irow))*cos(AOA(1,irow)),-sin(AOA(2,irow))*sin(AOA(1,irow)); ...
                0,cos(AOA(2,irow))];
        end
        
        Km=zeros(3,2,srow);
        
        for irow=1:srow,
            Km(:,:,irow)=[-sin(AOA(2,irow))*sin(AOA(1,irow)),cos(AOA(2,irow))*cos(AOA(1,irow)); ...
                sin(AOA(2,irow))*cos(AOA(1,irow)),cos(AOA(2,irow))*sin(AOA(1,irow)); ...
                0,sin(AOA(2,irow))];
        end
        
        Jm=zeros(3,srow);
        
        for irow=1:srow,
            Jm(:,irow)=[cos(AOA(1,irow));sin(AOA(1,irow));0];
        end
        
        T=zeros(3*srow-1,3*srow-1);
        
        index1=1;
        index2=1;
        
        for irow=1:srow-1,
            
            T(index1,index1)=-(B(:,irow+1)-B(:,1)).'*B(:,1);
            T(index1,srow:srow+1)=(2*u-TDOA(irow)*(B(:,irow+1)-B(:,1))-(S(1+irow,:).'+S(1,:).'-TDOA(irow)*B(:,1))).'*Lm(:,:,1);
            
            T(index1,(srow+2+2*(index1-1)):(srow+3+2*(index1-1)))=(S(1+irow,:).'+S(1,:).'-TDOA(irow)*B(:,1)-2*u).'*Lm(:,:,1+irow);
            
            index1=index1+1;
        end
        
        for irow=1:srow,
            
            T(srow-1+index2,srow-1+index2)=(S(irow,:)-u.')*Jm(:,irow);
            T(srow-1+index2+1,srow-1+index2:srow-1+index2+1)=(S(irow,:)-u.')*Km(:,:,irow);
            
            index2=index2+2;
        end
        
        
        mse=inv((inv(T)*G)'*inv(Q)*inv(T)*G);
        
        mseresult(nfnum,il)=trace(mse);
        
        % comput the CRB.
        
        D=zeros(3*srow-1,3);
        
        index2=1;
        
        for irow=1:srow-1,
            D(irow,:)=[(u(1)-S(irow+1,1))/norm(u-S(irow+1,:).',2)-(u(1)-S(1,1))/norm(u-S(1,:).',2),(u(2)-S(irow+1,2))/norm(u-S(irow+1,:).',2)-(u(2)-S(1,2))/norm(u-S(1,:).',2),(u(3)-S(irow+1,3))/norm(u-S(irow+1,:).',2)-(u(3)-S(1,3))/norm(u-S(1,:).',2)];
        end
        
        for irow=1:srow,
            
            D(srow-1+index2,:)=[ -(u(2)-S(irow,2))/((u(1)-S(irow,1))^2+(u(2)-S(irow,2))^2),(u(1)-S(irow,1))/((u(1)-S(irow,1))^2+(u(2)-S(irow,2))^2),0];
            D(srow-1+index2+1,:)=[-(u(1)-S(irow,1))*(u(3)-S(irow,3))/((u-S(irow,:).')'*(u-S(irow,:).')*sqrt((u(1)-S(irow,1))^2+(u(2)-S(irow,2))^2)),-(u(2)-S(irow,2))*(u(3)-S(irow,3))/((u-S(irow,:).')'*(u-S(irow,:).')*sqrt((u(1)-S(irow,1))^2+(u(2)-S(irow,2))^2)),sqrt((u(1)-S(irow,1))^2+(u(2)-S(irow,2))^2)/((u-S(irow,:).')'*(u-S(irow,:).'))];
            
            index2=index2+2;
        end
        
        FIM=D'*inv(Q)*D;
        
        CRB=inv(FIM);
        
        CRBresult(nfnum,il)=trace(CRB);
        
        % the ensemble of Monte-Carlo
        
        ebias=zeros(length(u),1);
        
        ebias1=zeros(length(u),1);
        
        for ensemnum=1:L,
            
            % the TDOA estimate.
            eTDOA=zeros(srow-1,1);
            
            for irow=1:srow-1,
                eTDOA(irow,1)=TDOA(irow,1)+sqrt(sigmaTDOA)*tdoanoise(ensemnum,irow);
                %                 eTDOA(irow,1)=TDOA(irow,1)+sqrt(sigmaTDOA)*randn(1);
            end
            
            % the AOA estimate.
            eAOA=zeros(2,srow);
            
            for irow=1:srow,
                eAOA(1,irow)=AOA(1,irow)+sqrt(sigmaAOA)*aoanoise1(ensemnum,irow);
                %                 eAOA(1,irow)=AOA(1,irow)+sqrt(sigmaAOA)*randn(1);
                eAOA(2,irow)=AOA(2,irow)+sqrt(sigmaAOA)*aoanoise2(ensemnum,irow);
                %                 eAOA(2,irow)=AOA(2,irow)+sqrt(sigmaAOA)*randn(1);
            end
            
            eG=zeros(3*srow-1,3);
            
            rowindex1=1;
            rowindex2=1;
            
            for irow=1:srow-1,
                
                eG(rowindex1,:)=2*[cos(eAOA(2,irow+1))*cos(eAOA(1,irow+1))-cos(eAOA(2,1))*cos(eAOA(1,1)),cos(eAOA(2,irow+1))*sin(eAOA(1,irow+1))-cos(eAOA(2,1))*sin(eAOA(1,1)),sin(eAOA(2,irow+1))-sin(eAOA(2,1))];
                
                rowindex1=rowindex1+1;
            end
            
            for irow=1:srow,
                
                eG(srow-1+rowindex2,:)=[sin(eAOA(1,irow)),-cos(eAOA(1,irow)),0];
                eG(srow-1+rowindex2+1,:)=[sin(eAOA(2,irow))*cos(eAOA(1,irow)),sin(eAOA(2,irow))*sin(eAOA(1,irow)),-cos(eAOA(2,irow))];
                
                rowindex2=rowindex2+2;
            end
            
            eB=zeros(3,srow);
            
            for irow=1:srow,
                eB(:,irow)=[cos(eAOA(2,irow))*cos(eAOA(1,irow));cos(eAOA(2,irow))*sin(eAOA(1,irow));sin(eAOA(2,irow))];
            end
            
            eh=zeros(3*srow-1,1);
            
            for irow=1:srow-1,
                eh(irow,1)=eG(irow,:)/2*(S(irow+1,:).'+S(1,:).'-eTDOA(irow)*eB(:,1));
            end
            
            index1=1;
            index2=1;
            
            for irow=srow:2:3*srow-1,
                eh(srow-1+index2,1)=eG(irow,:)*S(index1,:).';
                eh(srow-1+index2+1,1)=eG(irow+1,:)*S(index1,:).';
                
                index1=index1+1;
                index2=index2+2;
            end
            
            % initial estimate.
            eW=eye(3*srow-1);
            eu=inv(eG'*inv(eW)*eG)*eG'*inv(eW)*eh;
            
            eLm=zeros(3,2,srow);
            
            for irow=1:srow,
                eLm(:,:,irow)=[-cos(eAOA(2,irow))*sin(eAOA(1,irow)),-sin(eAOA(2,irow))*cos(eAOA(1,irow)); ...
                    cos(eAOA(2,irow))*cos(eAOA(1,irow)),-sin(eAOA(2,irow))*sin(eAOA(1,irow)); ...
                    0,cos(eAOA(2,irow))];
            end
            
            eT=zeros(3*srow-1,3*srow-1);
            
            for irow=1:srow-1,
                eT(irow,irow)=-(eB(:,irow+1)-eB(:,1)).'*eB(:,1);
            end
            
            for inum=1:2,
                % the estimate of the source's position.
                
                index1=1;
                
                for irow=1:srow-1,
                    
                    eT(index1,srow:srow+1)=norm(eu-S(1,:).',2)*eB(:,1+irow).'*eLm(:,:,1);
                    eT(index1,(srow+2+2*(index1-1)):(srow+3+2*(index1-1)))=-norm(eu-S(1+irow,:).',2)*eB(:,1).'*eLm(:,:,1+irow);
                    
                    index1=index1+1;
                end
                
                index1=1;
                
                for irow=srow:2:3*srow-1,
                    eT(irow,irow)=-norm(eu-S(index1,:).',2)*cos(eAOA(2,index1));
                    eT(irow+1,irow+1)=-norm(eu-S(index1,:).',2);
                    
                    index1=index1+1;
                end
                
                eW=eT*Q*eT';
                
                eu=inv(eG'*inv(eW)*eG)*eG'*inv(eW)*eh;
            end
            
            % for the improved solution of the proposed method.
            eresult(nfnum,il)=eresult(nfnum,il)+(eu-u)'*(eu-u);
            ebias=ebias+(eu-u);
            
            % %             % for the solution by good initial guess.
            % %             eu1=tdoaaoataylorseries(S,Q,eTDOA,eAOA,u,10^-5);
            % %             eresult1(nfnum,il)=eresult1(nfnum,il)+(eu1-u)'*(eu1-u);
            % %             ebias1=ebias1+(eu1-u);
            
        end
        
        
        % mse results of the proposed method.
        eresult(nfnum,il)=eresult(nfnum,il)/L;
        
        % %         % mse results of the Taylor-series Method.
        % %         eresult1(nfnum,il)=eresult1(nfnum,il)/L;
        
        % bias results of the proposed method.
        sebias(nfnum,il)=sum((ebias/L).^2);
        
        % %         % bias results of the Taylor-series Method.
        % %         sebias1(nfnum,il)=sum((ebias1/L).^2);
        
    end
    
end

for nfnum=1:length(nf),
    
    teresult(nfnum)=sum(eresult(nfnum,:))/L1;
    
    % %     teresult1(nfnum)=sum(eresult1(nfnum,:))/L1;
    % %
    tmseresult(nfnum)=sum(mseresult(nfnum,:))/L1;
    
    tCRBresult(nfnum)=sum(CRBresult(nfnum,:))/L1;
    
    tsebias(nfnum)=sum(sebias(nfnum,:))/L1;
    
    % %     tsebias1(nfnum)=sum(sebias1(nfnum,:))/L1;
end

fprintf('\n');

teresult=sqrt(teresult);

tmseresult=sqrt(tmseresult);

tCRBresult=sqrt(tCRBresult);

tsebias=sqrt(tsebias);

figure;

if (VryNseTDOA)
    tdoaErr=sqrt(tdoanf);
    loglog(tdoaErr,teresult,'ko',tdoaErr,tCRBresult,'k-',tdoaErr,tsebias,'k+','LineWidth',1.5)
    xlabel('\sigma_R_D (m)')
    ylabel('RMSE and bias (m)')
    legend('RMSE of the proposed method','Root CRB','Bias of the proposed method')
    xlim([10^0.5 1000]);
    
else
    aoaErr=sqrt(aoanf)/pi*180;
    semilogy(aoaErr,teresult,'ko',aoaErr,tCRBresult,'k-',aoaErr,tsebias,'k+','LineWidth',1.5)
    xlabel('\sigma_A_O_A (degrees)')
    ylabel('RMSE and bias (m)')
    legend('RMSE of the proposed method','Root CRB','Bias of the proposed method')
    xlim([min(aoaErr) max(aoaErr)]);
    ylim([0.8e-1 1000]);
end;