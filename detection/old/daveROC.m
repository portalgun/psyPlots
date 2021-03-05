function [PC,dp,idxM]=daveROC(DV_L,DV_R,X,minMax,nSmp,nIntrvl,bPlot)
% computes percent correct, d', and optimal criterion from p(LLR|A) and p(LLR|B), and plots ROC curve
%TODO
%   ADD INTERPOLATION in neg diag and best criterion
% ===============================================================================
    range=linspace(minMax(1),minMax(2),nSmp);
    tPr_ALL=zeros(length(range),1);
    fPr_ALL=zeros(length(range),1);
    %SWEEP CRITERION
    [PC,dp,idxM]=dave_AUC(fPr_ALL,tPr_ALL,nIntrvl,bPlot);
end
% ===============================================================================
