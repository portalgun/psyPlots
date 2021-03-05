function [] = expPsyThreshSepPlot(color,shape,errorbars,fsLbl,fsTck,xtckSD,ytckSD,subjs,nStds,XSS,SpcData,StFit,StCI,fldNames,stdXunqAll)
if isempty(Plot)
	color='k';
end
if isempty(shape)
	shape='o';
end
%plot thresholds

%figure(nFn);

hold on

%FIND MIN MAX X
firstP=min(stdXunqAll);
lastP =max(stdXunqAll);
r=abs(lastP-firstP);


%X'S
X=zeros(length(StFit),1);
for i=1:nStds
    X(i)=XSS.(fldNames{i}).stdX(1);
end
[X,idx]=sort(X);

%FIND DISTANCE BETWEEN STDS
xdist=r;
%if length(X)==1
%    xdist=r;
%else
%    xdist=abs(X(1)-X(2));
%end

%Y'S
Y=StFit;
Y=Y(idx);
disp(subjs)
for i = 1:length(Y)
    disp(['   ' num2str(Y(i))]);
end

%CONFIDENCE INTERVALS
CIn=zeros(nStds,1);
CIp=zeros(nStds,1);
for i=1:nStds
    CIn(i)=StCI{i}(1);
    CIp(i)=StCI{i}(2);
end
CIn=CIn(idx);
CIp=CIp(idx);

% -------------------------------------------------------------------------------

%PLOT ERROR BARDS
if errorbars==1
    for i = 1:nStds
        plot([X(i),X(i)],[CIn(i),CIp(i)],'black','linewidth',2,'HandleVisibility','off'); hold on;
    end
elseif errorbars==2
    plotfillederror(X,CIn',CIp',color,0); hold on;
    for i = 1:nStds
        plot([X(i),X(i)],[CIn(i),CIp(i)],'black','linewidth',2,'HandleVisibility','off'); hold on;
    end
else
    plotfillederror(X,CIn',CIp',color,0); hold on;
end

prop=['k' shape];
plot(X,Y,prop,'linewidth',2,'MarkerFaceColor','w');

xtck=stdXunqAll;
ytck=[.3,1,3,10];
%r=abs(max(X)-min(X));
%xlim([ min(X)-r*.15,max(X)+r*.15]);
xlim([min(stdXunqAll)-.25*xdist,max(stdXunqAll)+.25*xdist]);
ylim([min(ytck),max(ytck)]);
xlbl=('Disparity (arcmin)');
ylbl=('Threshold (arcmin)');

axis square
formatFigure(xlbl,ylbl,[],0,1,fsLbl,fsTck,xtck,ytck,xtckSD,ytckSD);
hold off
