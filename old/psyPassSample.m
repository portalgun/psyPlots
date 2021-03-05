function [DataSmp] = psyPassSample(Data,n,rndSd)

if ~exist('rndSd','var') || isempty(rndSd)
    rng('shuffle');
else
    setRndSd(rndSd);
end

ind=randperm(size(Data,1),n);
DataSmp=Data(ind,:);
