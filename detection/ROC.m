classdef ROC < handle
properties
    range
    PC
    dp
    idxM
end
methods
    function obj=ROC(DV_L,DV_R,X,minMax,nSmp,nIntrvl)
        if size(DV_L,1)==1
            DV_L=DV_L';
            DV_R=DV_R';
        end
        A=[DV_L; DV_R]
        mini=min(A);
        maxi=max(A);
        n=numel(unique(A))

        if exist('X','var') && ~isempty(X)
            obj.X=X;
        else
            obj.X=linspace(mini,maxi,n);
        end
        if exist('nSmp','var') && ~isempty(nSmp)
            obj.nSmp=nSmp;
        else
            obj.nSmp=n;
        end
        if exist('minMax','var') && ~isempty(minMax)
            obj.minMax=minmax;
        else
            obj.minMax=[mini maxi];
        end

        obj.DV_L=DV_L;
        obj.DV_R=DV_R;
        obj.nIntervl=nIntrvl;
        obj.range=linspace(minMax(1),minMax(2),nSmp);
        obj.tPr=zeros(length(range),1);
        obj.fPr=zeros(length(range),1);

        obj.get_all_ratio();
        obj.AUC();
    end
    function obj=get_all_ratio(obj)
        for i = 1:length(range)
            criterion=obj.range(i);
            [obj.tPr(i),obj.fPr(i)]=obj.get_ratio(criterion);
        end
    end
    function [fPr,tPr]= get_ratio(obj,criterion)
        tP=sum(DV_R(X>=criterion));
        fP=sum(DV_L(X>=criterion));

        fN=sum(DV_R(X<criterion));
        tN=sum(DV_L(X<criterion));

        tPr=tP./(tP+fN);
        fPr=fP./(fP+tN);
    end
    function obj = AUC(obj)
        Rmax=-1;
        Mmin=Inf;
        for i = length(obj.tPr):-1:1
            %FIND BEST CRITERION
            R=obj.tPr(i)-obj.fPr(i);
            if R > Rmax
                Rmax=R;
                idxR=i;
            end
            %FIND CENTER OF NEGATIVE DIAGONAL
            M=sqrt((obj.tPr(i)-1).^2+obj.fPr(i).^2);
            if M < Mmin
                Mmin=M;
                idxM=i;
            end
        end
        maxTPr=obj.tPr(idxR);
        minFPr=obj.fPr(idxR);
        cTPr=obj.tPr(idxM);
        cFPr=obj.fPr(idxM);
        c1AFCx=[0 cFPr 1];
        c1AFCy=[0 cTPr 1];

        %2AFC vs 1AFC
        if obj.nIntrvl == 2
            obj.AUC=abs(trapz(obj.fPr,obj.tPr));
            dp=norminv(obj.AUC)*sqrt(2);
        elseif obj.nIntrvl == 1
            obj.AUC=abs(trapz(c1AFCx,c1AFCy));
            dp=norminv(obj.AUC)*2;
        end
        obj.PC=obj.AUC*100;

    end
    function []=plot(obj)
        Fig.new()
        plot(obj.fPr_ALL,obj.tPr_ALL,'LineWidth',2);hold on
        plot([0,1],[0,1],'k--');
        plot(c1AFCx,c1AFCy,'k');
        scatter(minFPr,maxTPr,'r')
        xlabel('False Alarm Rate (FP)')
        ylabel('Hit Rate (TP)')
        xlim([0,1]);
        ylim([0,1]);
        if obj.nIntrvl == 2
            text(.7,.3,['2AFC' newline  'PC = ' num2str(round(obj.PC,2))  newline 'd'' = ' num2str(round(dp,2))])
        elseif obj.nIntrvl == 1
            text(.7,.3,['1AFC' newline  'PC = ' num2str(round(obj.PC,2))  newline 'd'' = ' num2str(round(dp,2))])
        end
        axis square
        hold off
    end
end
end
