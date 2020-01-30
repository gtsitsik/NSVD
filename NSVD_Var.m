% Calculates the singular values required for the calculation of NSVD.
% Also calculates CORCONDIA, MSE and Missing Data MSE.
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
%
%--------------------------------------------------------------------------
% OUTPUTS
%   sing : singular values of the Khatri-Rao product of the PARAFAC factor
%          matrices.
%   corcond : CORCONDIA values
%   error : MSE values
%   missErrAll : Missing Data MSE values
%   cpIters : Iterations of the optimization algorithm

function [sing, corcond, error, missErrAll, cpIters] = NSVD_Var(X,FRange,coresNum,nsvdIt,cpTol,cpMaxIt,cpAlgo,compRatio,compStep)



if nargin<1
    % rank-5 SynthSmall tensor
    load SynthSmall.mat
    X = data;
    % X = artificial_data_generator([20 20 20],4,100,5,0);
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


if cpAlgo == 1
    Options(2)=2;
elseif cpAlgo == 0
    Options(2)=2;
elseif cpAlgo == 2
    X = double(X);
end



Options(1) = cpTol;
Options(5) = NaN;
Options(6) = cpMaxIt;

miss_ind_all = [];
XComp=[];
sing = cell(1,length(FRange));

if cpAlgo == 2
    for i = 1:nsvdIt
        tensorSize = prod(size(X));
        miss_ind_all(i,:) = logical([ones(1,ceil(tensorSize*0.8)) zeros(1,floor(tensorSize*0.2))]);
        miss_ind_all(i,:) = logical(miss_ind_all(i,randperm(tensorSize)));
    end
end


for i = 1:nsvdIt
    
    fprintf("\tSample %d\n", i);
    miss_ind_cur=0;
    if cpAlgo == 2
        miss_ind_cur = logical(miss_ind_all(i,:));
    end
    
    % Orthonormal Compression
    if(compRatio>=0 && (i-1)/compStep==ceil((i-1)/compStep))
        sz = size(X);
        U = randn(ceil(sz(1)*compRatio/100),sz(1));
        V = randn(ceil(sz(2)*compRatio/100),sz(2));
        W = randn(ceil(sz(3)*compRatio/100),sz(3));
        U = orth(U')';
        V = orth(V')';
        W = orth(W')';
        if strcmp(class(X),'ktensor')
            XComp = ttm(X,{U,V,W},[1 2 3]);
        else
            XComp = ttm(tensor(X),{U,V,W},[1 2 3]);
            XComp = XComp.data;
        end
    end
    
    parfor(k=1:length(FRange) , coresNum)
        
        singul = [];
        cp=[];
        F=FRange(k) ;
        
        if cpAlgo == 2
            X2 = NaN(size(X));
            X2(miss_ind_cur) = X(miss_ind_cur);
        else
            X2 = X;
        end
        if(compRatio>=0)
            X2 = XComp;
        end
        
        
        %***********************  PARAFAC ********************************
        if cpAlgo == 1
            evalc('[cp,cpIt] = parafac(X2,F,Options)');
            cpItersTmp(k) = cpIt;
        elseif cpAlgo == 0
            [cpp,~,out] = cp_als(tensor(X2),F,'tol',cpTol,'maxiters',cpMaxIt);
            cp{1}= cpp.U{1};
            cp{2}= cpp.U{2};
            cp{3}= cpp.U{3};
            cpIt = out.iters;
            cpItersTmp(k) = out.iters;
        elseif cpAlgo == 2
             evalc('[cp,cpIt] = parafac(X2,F,Options)');
            cpItersTmp(k) = cpIt;
        end
        
        cpAll{i,k}= cp;
        kr = khatrirao(cp{3},cp{2},cp{1});
        singul = svd(kr);
        
        
        %******************** Basic Calculations **************************
        sing{k}(i,:) = singul;
        decTensor = ktensor(ones(F,1),cp{1},cp{2},cp{3});
        error(i,k) = norm(full(X2)-full(decTensor))^2;
        corcond(i,k) = efficient_corcondia(tensor(X2), decTensor);
        if cpAlgo==2
            tmp_miss_err = X - double(ktensor(cpAll{i,k}));
            missErrAll(i,k) = norm(tmp_miss_err(logical(-(miss_ind_all(i,:)'-1))))^2;
        end
        
        
    end
    cpIters(i,:) = cpItersTmp;
end


end



