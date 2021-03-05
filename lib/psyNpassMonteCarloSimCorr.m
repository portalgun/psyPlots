function [PA, PC] = psyNpassMonteCarloSimCorr(N,M,mu,cov,numSim,numTrlPerSim)

% function [PA, PC] = psyNpassMonteCarloSimCorr(N,M,mu,cov,numSim,numTrlPerSim)
%
% example call: [PA, PC] = psyNpassMonteCarloSimCorr(2,2,repmat([-2:0.05:2]',[1 2]),[1 0; 0 1],10000,100)
%
% simulates 2AFC trials for which the observer's decision variable is
% Gaussian distributed, then computes proportion agreement and proportion
% cmp chosen
% 
% inputs: 
%         N           : total number of passes
%         M           : M-agreement (e.g. if there are 8 passes, we might be
%                       interested in the probability of agreement on 6 of those
%                       passes. In that case, N = 8, M = 6.
%         mu          : means of decision variable: numRows = number of
%                       means to try, numCols = N
%         cov         : covariance matrix for difference signal
%         numSim      : number of simulations to generate
%         numTrlPerSim: number of trials per simulation
%
% outputs: 
%          PA  : proportion agreement
%          PC  : proportion comparison chosen

% CHECK THAT M DOESN'T EXCEED N
if M>N
   error('psyNpassMonteCarloSimCorr0: M cannot exceed N'); 
end

% CHECK THAT mu DIMENSIONS ARE MATCHED TO N
if N ~=size(mu,2)
   error('psyNpassMonteCarloSimCorr0: make sure mu has correct number of entries!'); 
end

% CHECK THAT cov DIMENSIONS ARE MATCHED TO N
if N ~=size(cov,2)
   error('psyNpassMonteCarloSimCorr0: make sure sigma has correct number of entries!'); 
end

if sum(cov(boolean(1-eye(N))))==0 % IF THERE IS NO CORRELATION BETWEEN DECISION VARIABLES, CAN COMPUTE FAST
    % GENERATE DECISION VARIABLES
    Z = normrnd(mu(1),cov(1,1),[numTrlPerSim N numSim]);
    % COMPUTE PROPORTION CMP CHOSEN
    PC = sum(reshape(Z,[size(Z,1)*size(Z,2) size(Z,3)])>0)'./(numTrlPerSim*N);
    % COMPUTE AGREEMENT
    PA = squeeze(sum(sum(Z>=0,2)==M | sum(Z<0,2)==M)./numTrlPerSim);
else % IF THERE IS CORRELATION BETWEEN DECISION VARIABLES, LOOP (SLOWER)
    for i = 1:numSim % FOR EACH SIMULATION
       % GENERATE DECISION VARIABLES
       Z = mvnrnd(mu,cov,numTrlPerSim);
       % COMPUTE PC
       PC(i,:) = sum(Z(:)>0)./numel(Z);
       % COMPUTE M-AGREEMENT
       bZmoreThan0 = Z>=0;
       bZlessThan0 = Z<0;
       PA(i,:) = sum(sum(bZmoreThan0,2)==M | sum(bZlessThan0,2)==M)./numTrlPerSim;
    end
end

end