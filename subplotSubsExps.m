function [RC,r,c]=subplotSubsExps(fig_num,subs,subsAll,subPlotXSub)
    %subPlotXsub, which subsAll dim is on x axis

    %GET SUBPLOT SUBSCRIPTS
    nUnCol=sum(diff(sort(subsAll))~=0)+1;
    nF=unique(subsAll(:,2:end),'rows');
    nL=unique(subsAll(:,1:end-1),'rows');
    if ~isvar('subPlotXSub')
        subPlotXSub=1;
    end

    v=subPlotXSub;
    if v==1
        [~,i]=unique(subsAll(:,1));
        Rvals=subsAll(i,1);

        RC=[nUnCol(1) size(nF,1)];
        r=find(Rvals == subs(1));
        if size(nF,1)==1
            c=1;
        else
            c=findrow([subs(2:end)],nF);
        end
    elseif v==size(subsAll,2)
        [~,i]=unique(subsAll(:,end));
        Rvals=subsAll(i,end);

        RC=[nUnCol(end) size(nL,1)];
        r=find(Rvals == subs(end));
        if size(nF,1)==1
            c=1;
        else
            c=findrow(subs(1:end-1),nF);
        end
    else
        V=[1:v-1, v+1:size(subsAll,2)];
        nV=unique(subsAll(:,V),'rows');

        [~,i]=unique(subsAll(:,v));
        Rvals=subsAll(i,v);

        RC=[nUnCol(v) size(nV,1)];
        r=find(Rvals == subs(v));
        if size(nF,1)==1
            c=1;
        else
            c=findrow(sub(V),nF);
        end
    end
    figure(fig_num)
    subPlot(RC,r,c);
end
