function [yl] = expPsyThreshBlk(flagidx,range,shape,errorbars,fsLbl,fsTck,xtckSD,ytckSD,stdX,S,pcData,tFit,tCI,subStr,bYlabel,bXlabel)
%plot thresholds
    if ~exist('bYlabel','var')
        bYlabel=1;
    end
    if ~exist('bXlabel','var')
        bXlabel=1;
    end


    l=length(tFit);
    X=range;
    Y=tFit;
    %Confidence Intervals
    CIn=zeros(l,1);
    CIp=zeros(l,1);
    for i=1:l
        CIn(i)=tCI{i}(1);
        CIp(i)=tCI{i}(2);
    end

    if errorbars==1
        for i = 1:l
            plot([X(i),X(i)],[CIn(i),CIp(i)],'black','linewidth',2); hold on;
        end
    elseif errorbars==2
        plotfillederror(X,CIn',CIp',[.8,.8,.8],1); hold on;
        for i = 1:l
            if ismember(i,flagidx)
                plot([X(i),X(i)],[CIn(i),CIp(i)],'red','linewidth',2); hold on;
            else
                plot([X(i),X(i)],[CIn(i),CIp(i)],'black','linewidth',2); hold on;
            end
        end
    elseif errorbars==-1
        ;
    else
        plotfillederror(X,CIn',CIp',[.2,.2,.2],1); hold on;
    end
    %Plot means, and mark by color
    for i = 1:length(X)
        hold on
        if ismember(i,flagidx)
            plot(X(i),Y(i),'rs','linewidth',2,'MarkerFaceColor','w');
        else
            plot(X(i),Y(i),'ks','linewidth',2,'MarkerFaceColor','w');
        end
    end

    if length(X) > 10
        xtck=  X(1:2:end);
    else
        xtck=  X ;
    end
    xlim([X(1)-.5, X(end)+.5]);
    %[CIn,CIp]
    yl=[min(floor(CIn*2)/2)-.25, max(ceil(CIp*2)/2)+.25];

    ytck= floor(min(yl)):.5:ceil(max(yl));
    if length(ytck) > 9
        ytck= floor(min(yl)):2:ceil(max(yl));
    end
    if length(ytck) > 9
        ytck= floor(min(yl)):1:ceil(max(yl));
    end
    if exist('subStr','var')==1 && ~isempty(subStr)
        ylbl=({subStr;'Threshold (arcmin)'});
    elseif bYlabel==1
        ylbl=('Threshold (arcmin)');
    elseif bYlabel==0
        ylbl=('');
    end
    if bXlabel==1
        xlbl=('Block Number');
    else
        xlbl=('');
    end
    ttl= [ num2str(stdX) ' arcmins' ];

    axis square
    formatFigure(xlbl,ylbl,ttl,0,0,fsLbl,fsTck,xtck,ytck,xtckSD,ytckSD);
    hold on
end
