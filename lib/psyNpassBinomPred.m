function [PA, PC, PAlb, PAub] = psyNpassBinomPred(N,M,mu,cov,CIquant,numSim,numTrlPerSim,bPlot,bQuiet)

% function [PA, PC, PAlb, PAub] = psyNpassBinomPred(N,M,mu,cov,CIquant,numSim,numTrlPerSim,bPlot)
%
% example call: [PA, PC, PAlb, PAub] = psyNpassBinomPred(2,2,repmat([-2:0.2:2]',[1 2]),[1 0; 0 1],[0.05 0.95],10000,100,1)
%
% compute predictions of agreement with confidence intervals, for Gaussian
% distributed decision variable
%
% inputs:
%         N           : total number of passes
%         M           : M-agreement (e.g. if there are 8 passes, we might be
%                       interested in the probability of agreement on 6 of those
%                       passes. In that case, N = 8, M = 6.
%         mu          : means of decision variable: numRows = number of
%                       means to try, numCols = N
%         cov         : covariance matrix for difference signal
%         CIquant     : quantiles for confidence intervals
%         numSim      : number of simulations to generate
%         numTrlPerSim: number of trials per simulation
%
% outputs:
%          PA  : proportion agreement
%          PC  : proportion comparison chosen
%          PAlb: proportion agreement lower bounds
%          PAub: proportion agreement upper bounds
if ~exist('bQuiet','var') || isempty(bQuiet)
    bQuiet = 0;
end

for i = 1:length(mu) % FOR EACH MEAN
   % COMPUTE AGREEMENT AND CMP CHOSEN MONTE CARLO SIMULATIONS
   [PAtmp, PCtmp] = psyNpassMonteCarloSimCorr(N,M,mu(i,:),cov,numSim,numTrlPerSim);
   % STORE ACTUAL P(A) AND P(C)
   PC(i,:) = mean(PCtmp);
   PA(i,:) = mean(PAtmp);
   % STORE LOWER BOUNDS ON P(A)
   PAlb(i,:) = quantile(PAtmp,CIquant(1));
   PAub(i,:) = quantile(PAtmp,CIquant(2));
   if bQuiet == 0
       progressreport(i,10,length(mu));
   end
end

if bPlot==1
    figure;
    fill([PC; flipud(PC)],[PAlb; PAub],0.9.*[1 1 1],'EdgeColor','none'); hold on;
    plot(PC,PA,'--k','LineWidth',1.5); hold on;
    axis square;
    Fig.format('P(C)','P(A)');
    set(gca,'LineWidth',1.5);
end

end
