function [DV_R,DV_L,X] = daveSmp2LLRgauss(zA,zB,minMax,nBin,nIntrvl,bPlot)
%za?
%zb?
%minMax bin range
%nBins
%nIntrvl - 2AFC etc

    if exist('bPlot','var')~=1;   bPlot=0;   end
    if exist('nIntrvl','var')~=1; nIntrvl=2; end

    bins = linspace(minMax(1),minMax(2),nBin);

    mA=mean(zA);
    mB=mean(zB);

    varA=var(zA);
    varB=var(zB);

    %p(data is A | B data) / p(data is B | B data)
    DVb=log(mvnpdf(zB,mA,varA)./mvnpdf(zB,mB,varB));
    %p(data is A | A data) / p(data is B | A data)
    DVa=log(mvnpdf(zA,mA,varA)./mvnpdf(zA,mB,varB));
    hDV_A=hist(DVa, bins);
    hDV_B=hist(DVb, bins);

    DV_R = hDV_B./sum(hDV_B(:));
    DV_L = hDV_A./sum(hDV_A(:));
    X=linspace(minMax(1),minMax(2),length(DV_R));

    if bPlot==1
        [~,~,idxM] = daveROC(DV_L,DV_R,X,minMax,length(DV_L),nIntrvl,0);
        Fig.new()
            plot(DV_R); hold on
            plot(DV_L)
            scatter(X(idxM),0,'kX','LineWidth',2)
            legend({'$P(\Psi|a)$','$P(\Psi|b)$',},'Location','northwest','Interpreter','LaTeX')
            xlabel('Decision Variable')
            ylabel('Probability')
        hold off
    end
end
