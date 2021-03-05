function [] = expPsyCurvePlotAll(fsLbl,fsTck,xtckSD,ytckSD,shape,colors,XS,fldNames,pcData,mFit,sFit,tFit,bFit,bSubPlot,xlbl,ylbl,ttl)
%plot all data points and curves


%COLORS
%colors=[[109,174,190];[106,200,168];[103,177,110];[160,200,106];[190,182,86]]; %pastel Green scheme
%colors=[[2,98,177];[9,187,173];[0,164,59];[43,200,9];[167,177,23]]; %Dark Green

% -------------------------------------------------------------------------------
nStds=length(fldNames);
%FIND MIN MAX X
firstP=min(min(XS.(fldNames{1}).cmpXunq),min(pcData{1}));
lastP =max(max(XS.(fldNames{end}).cmpXunq),max(pcData{end}));
r=abs(lastP-firstP);

%X's
X=zeros(size(tFit,1),1);
for i=1:nStds
    X(i)=XS.(fldNames{i}).stdX(1);
end
[X,~]=sort(X);

%FIND DISTANCE BETWEEN STDS
if length(X)==1
    xdist=r;
else
    xdist=abs(X(1)-X(2));
end

%extrapolate from fit
cmpExtr=cell(length(nStds)+4);
cFit=cell(length(nStds)+4);
for i = 1:nStds
    cmp=XS.(fldNames{i}).cmpXunq;
    cmpR=abs(cmp(end)-cmp(1));
    cmpExtr{i}=[ mFit{i}-cmpR*1.5; cmp; mFit{i}+cmpR*1.5 ];
    cFit{i}=(normcdf(cmpExtr{i},mFit{i},(sFit{i})));
end

% -------------------------------------------------------------------------------
hold off;

if bSubPlot == 0
    %SUBJECT MARKERS
    scatter(1000,1000,shape,'MarkerEdgeColor',[0,0,0],'MarkerFaceColor',[1,1,1]); hold on

    %DATA POINTS
    for i=1:nStds
        scatter(XS.(fldNames{i}).cmpXunq,pcData{i}, ...
            shape,'MarkerEdgeColor',colors(i,:),'MarkerFaceColor',[1,1,1],'linewidth',2); hold on;
    end
end

%FITS
for i=1:nStds
    datapoints=linspace(-20,0,201);

    if bSubPlot == 1
        subplot(nStds,1,i)
        %SUBJECT MARKERS
        scatter(1000,1000,shape,'MarkerEdgeColor',[0,0,0],'MarkerFaceColor',[1,1,1]); hold on

        %DATA POINTS
        scatter(XS.(fldNames{i}).cmpXunq,pcData{i}, ...
        shape,'MarkerEdgeColor','k','MarkerFaceColor',[1,1,1],'linewidth',2); hold on;
    end

    %CALCULATE CUMULATIVE GAUSSIAN FIT
    %[PC,~,~]=psyfuncgengauss(XS.(fldNames{i}).stdX,datapoints,mFit{i},sFit{i},bFit{i},[],0);
    %XXX CHANGE NUMBER OF INTERVALS TO BE A VARIABLE HANDLED
    [PC]=psyfitgengaussfunc(XS.(fldNames{i}).stdX(1),datapoints,mFit{i},sFit{i},bFit{i},1.36,2,0);
    hold on

    if bSubPlot == 0
        %PLOT FIT
        plot(datapoints,PC,'Color',colors(i,:),'linewidth',2,'HandleVisibility','off')
        hold on;
    elseif bSubPlot == 1
        %PLOT FIT
        plot(datapoints,PC,'Color','k','linewidth',2,'HandleVisibility','off')

        %FORMAT SUBPLOT
        xlbl=('Disparity (arcmin)');
        ylbl=('% Comparison Chosen');

        xtck=[X(i)-xdist X(i) X(i)+xdist];
        ytck=0:.2:1;

        %xlim([ firstP-r*.15, lastP+r*.15]);
        %xlim([min(X)-xdist*1.5,max(X)+xdist*1.5]);
        xlim([-3 3]+X(i))
        ylim([min(ytck),max(ytck)]);

        axis square
        ttl='All Blocks Combined';
        formatFigure(xlbl,ylbl,ttl,0,0,fsLbl,fsTck,xtck,ytck,xtckSD,ytckSD);
        hold off
    end
end

if bSubPlot == 0

    xtck=min(X)-xdist*2:xdist:max(X)+xdist*3;
    ytck=0:.2:1;

    %xlim([ firstP-r*.15, lastP+r*.15]);
    xlim([min(X)-xdist*1.5,max(X)+xdist*1.5]);
    ylim([min(ytck),max(ytck)]);

    if nStds==1
        pbaspect([1.5 1 1]);
    else
        pbaspect([2.5 1 1]);
    end
    formatFigure(xlbl,ylbl,ttl,0,0,fsLbl,fsTck,xtck,ytck,xtckSD,ytckSD);
    hold off
end
