% Calculates and plots NSVD, CORCONDIA, MSE and Missing Data MSE.
%
%--------------------------------------------------------------------------
% INPUTS
%   X : Input 3-mode tensor
%   FRange : Range of CP-ranks
%   coresNum: maximum number of workers (requires Parallel Computing Toolbox)
%   nsvdIt : number of samples for each CP-rank
%   cpTol : Tolerance for the PARAFAC optimization algorithm
%   cpMaxIt : Maximum iterations for the PARAFAC optimization algorithm
%   cpAlgo : Type of optimization algorithm
%            0 -> PARAFAC of Tensor Toolbox
%            1 -> PARAFAC of N-way Toolbox
%            2 -> PARAFAC with missing data of N-way Toolbox
%                 (20% of data are randomly removed for all nsvdIt samples)
%   compRatio : -1 for no compression. Otherwise, the modes of the tensor
%               are scaled by compRatio
%   compStep : Number of compressions for each PARACAC sample

function NSVD_demo(X,FRange,coresNum,nsvdIt,cpTol,cpMaxIt,cpAlgo,compRatio,compStep)

if nargin<1
    % rank-5 SynthSmall tensor
    load SynthSmall.mat
    X = data;
end
if nargin<2
    FRange = 1:8;
end
if nargin<3
    coresNum = 0;
end
if nargin<4
    nsvdIt = 10;
end
if nargin<5
    cpTol = 1e-6;
end
if nargin<6
    cpMaxIt = 50000;
end
if nargin<7
    cpAlgo = 0;
end
if nargin<8
    compRatio = -1;
end
if nargin<9
    compStep = -1;
end



% Intermediate Calculations
fprintf("NSVD - CORCONDIA - MSE calculations\n")
[sing, corcond, error] = NSVD_Var(X,FRange,coresNum,nsvdIt,cpTol,cpMaxIt,cpAlgo,compRatio,compStep);
corcond(corcond<0)=0;

% fprintf("Missing Data MSE calculations\n")
% cpAlgo=2;
% [~,~,~,err_mis,~] = NSVD_Var(X,FRange,parallelMode,nsvdIt,cpTol,cpMaxIt,cpAlgo,compRatio,compStep);
% 


% Final Estimates
if compRatio<0
    for k = 1:nsvdIt-1
        for i = 1:length(FRange)
            NSVD_inter(k,i) = sum(log(var(filloutliers(sing{i}(1:k+1,:),'center'))./mean(filloutliers(sing{i}(1:k+1,:),'center'))));
            MSE_inter(k,i) = mean(filloutliers(error(1:k+1,i),'center'));
            CORCONDIA_inter(k,i) = mean(filloutliers(corcond(1:k+1,i),'center'));
%             MSE_miss_inter(k,i)=  mean(filloutliers(err_mis(1:k+1,i),'center'));
        end
    end
else
    for s = 1:nsvdIt/compStep
        k0=compStep*(s-1)+1;
        for k = k0:(k0-1)+compStep-1
            for i = 1:length(FRange)
                NSVD_inter(k-(k0-1),i,s) = sum(log(var(filloutliers(sing{i}(k0:k+1,:),'center'))./mean(filloutliers(sing{i}(k0:k+1,:),'center'))));
                MSE_inter(k-(k0-1),i,s) = mean(filloutliers(error(k0:k+1,i),'center'));
                CORCONDIA_inter(k-(k0-1),i,s) = mean(filloutliers(corcond(k0:k+1,i),'center'));
%                 MSE_miss_inter(k-(k0-1),i,s)=  mean(filloutliers(err_mis(k0:k+1,i),'center'));
            end
        end
    end
end






close all

% NSVD
NSVD_inter(abs(NSVD_inter)==inf)=NaN;
if compRatio>0
    for ind = 1: size(NSVD_inter,2)
        tmp = NSVD_inter(:,ind,:);
        tmp2(:,ind) = tmp(:);
    end
    NSVD_inter = tmp2;
end

pl4 = squeeze(nanmean(NSVD_inter,1))';
res = abs(quantile(NSVD_inter,[0.25 0.75])-[pl4';pl4']);

subplot(2,2,1)
errorbar(FRange,pl4,res(1,:),res(2,:));

axis tight
grid
xlabel('Number of Components')
ylabel('NSVD')
xticks(FRange)
xticklabels(FRange)




% MSE
subplot(2,2,2)
MSE_inter(abs(MSE_inter)==inf)=NaN;
if compRatio>0
    tmp2=[];
    for ind = 1: size(MSE_inter,2)
        tmp = MSE_inter(:,ind,:);
        tmp2(:,ind) = tmp(:);
    end
    MSE_inter = tmp2;
end

pll4 = squeeze(nanmean( MSE_inter,1))';
res = abs(quantile( MSE_inter,[0.25 0.75])-[pll4';pll4']);
errorbar(FRange,pll4,res(1,:),res(2,:))

axis tight
grid
xlabel('Number of Components')
ylabel('MSE')
xticks(FRange)
xticklabels(FRange)




% CORCONDIA
subplot(2,2,3)
CORCONDIA_inter(abs(CORCONDIA_inter)==inf)=NaN;
if compRatio>0
    tmp2=[];
    for ind = 1: size(NSVD_inter,2)
        tmp = CORCONDIA_inter(:,ind,:);
        tmp2(:,ind) = tmp(:);
    end
    CORCONDIA_inter = tmp2;
end

pll4 = squeeze(nanmean(CORCONDIA_inter,1))';
res = abs(quantile(CORCONDIA_inter,[0.25 0.75])-[pll4';pll4']);
errorbar(FRange,pll4,res(1,:),res(2,:));

axis tight
grid
xlabel('Number of Components')
ylabel('CORCONDIA')
xticks(FRange)
xticklabels(FRange)


% 
% % Missing Data MSE
% if compRatio<0
%     subplot(2,2,4)
%     
%     pll4 = squeeze(nanmean(MSE_miss_inter,1));
%     res = abs(quantile(MSE_miss_inter,[0.25 0.75])-[pll4;pll4]);
%     errorbar(FRange,pll4,res(1,:),res(2,:));
%     
%     axis tight
%     grid
%     xlabel('Number of Components')
%     ylabel('Missing Data MSE')
%     xticks(FRange)
%     xticklabels(FRange)
% end
