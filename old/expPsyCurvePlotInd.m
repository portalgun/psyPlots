function [] = expPsyCurvePlotInd(fsLbl,fsTck,xtckSD,ytckSD,shape,colors,XS,pcData,mFit,sFit,tFit,bFit,h,w,I,Z,ylbl)
%plot all data points and curves


%COLORS
%colors=[[109,174,190];[106,200,168];[103,177,110];[160,200,106];[190,182,86]]; %pastel Green scheme
%colors=[[2,98,177];[9,187,173];[0,164,59];[43,200,9];[167,177,23]]; %Dark Green

% -------------------------------------------------------------------------------
nBlks=length(mFit);
if exist('h','var')~=1 || exist('w','var')~=1
    h=floor(sqrt(length(mFit)));
    w=ceil(sqrt(length(mFit)));
    bMult=0;
else
    bMult=1;
end
if exist('I','var')~=1
    I=0;
end
if exist('Z','var')~=1
    Z=0;
end

%FIND MIN MAX X
firstP=min(XS.cmpXunq{1});
lastP =max(XS.cmpXunq{end});
r=abs(lastP-firstP);

%X's
X=zeros(size(tFit,1),1);
for i=1:nBlks
    X(i)=XS.stdX{1}(1);
end
[X,~]=sort(X);

%FIND DISTANCE BETWEEN STDS
if length(X)==1
    xdist=r;
else
    xdist=abs(X(1)-X(2));
end

%extrapolate from fit
cmpExtr=cell(length(nBlks)+4);
cFit=cell(length(nBlks)+4);
for i = 1:nBlks
    cmp=XS.cmpXunq{i};
    cmpR=abs(cmp(end)-cmp(1));
    cmpExtr{i}=[ mFit{i}-cmpR*1.5; cmp; mFit{i}+cmpR*1.5 ];
    cFit{i}=(normcdf(cmpExtr{i},mFit{i},(sFit{i})));
end

% -------------------------------------------------------------------------------

%SUBJECT MARKERS
hold on
scatter(1000,1000,shape,'MarkerEdgeColor',[0,0,0],'MarkerFaceColor',[1,1,1]); hold on

%DATA POINTS
for i=1:nBlks
    subplot(h,w,i+I)
    scatter(XS.cmpXunq{i},pcData{i}, ...
        shape,'MarkerEdgeColor',colors(i,:),'MarkerFaceColor',[1,1,1],'linewidth',2);
    hold on;
end

%FITS
datapoints=linspace(-20,0,200);
for i=1:nBlks
    [PC,~,~]=psyfuncgengauss(XS.stdX{i},datapoints,mFit{i},sFit{i},bFit{i},[],0);
    subplot(h,w,i+I)
    plot(datapoints,PC,'Color',colors(i,:),'linewidth',2,'HandleVisibility','off')
    hold on;

    if exist('ylbl','var')~=1 
        ylbl=('% Comparison Chosen');
    end

    xlbl=('Disparity (arcmin)');
    ttl=['Block ' num2str(i) ];
    ytck=0:.2:1;
    xtck=[min(XS.cmpXunq{1}),XS.stdX{1}(1),max(XS.cmpXunq{end})];
    pbaspect([1.5 1 1]);
    xlim([ firstP-r*.3, lastP+r*.3]);

    if i == 1 || mod(i,w)==1
        Ylbl=ylbl;
    else
        Ylbl='';
    end
    if Z==1 || I==0
        Ttl=ttl;
    else
        Ttl='';
    end
    if Z==h || bMult==0
        if i ==1 || mod(i,w)==1
            Xlbl=xlbl;
        else
            Xlbl='';
        end
    else
        Xlbl='';
    end
    formatFigure(Xlbl,Ylbl,Ttl,0,0,fsLbl,fsTck,xtck,ytck,xtckSD,ytckSD);
end



%xlim([min(X)-xdist*2.5,max(X)+xdist*1.5]);
%ylim([min(ytck),max(ytck)]);

hold off
