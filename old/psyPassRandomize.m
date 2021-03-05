function [DataNew] = psyPassRandomize(Data,rndSd)
% Data contains response values for a given comparison condition and standard
% Data [nTrl x nPass]

if ~exist('rndSd','var') || isempty(rndSd)
    rng('shuffle');
else
    setRndSd(rndSd);
end

ind=randperm(size(Data,2));
DataNew=Data(:,ind);
