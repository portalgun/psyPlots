function [DP,rho,ncmp,sutitle,signal,Data,rMean,yMean] = expDoublePass(D,prjCodes,expIDs,SUBJ,bSave)
if ~iscell(SUBJ)
	SUBJ={SUBJ};
end
% -------------------------------------------------------------------------------
%SORT DATA
[DATA,DP,cmpXsrt,SS,E,signal,subj] = expDoubleSort(D,prjCodes,expIDs,SUBJ);
% -------------------------------------------------------------------------------
%FIT DATA AND PLOT
bPLOT=1;
bBOOTSTRAPPING=0;
%plot(DATA(:,1),DATA(:,2),'ko','markerface','w','markersize',8);
tSubplot=length(D.idxStd)*2;
cSubplot=0;
rho=zeros(length(D.idxStd),E.ncmp);
mu1=zeros(length(D.idxStd),E.ncmp);
mu2=zeros(length(D.idxStd),E.ncmp);
rhoCI=zeros(length(D.idxStd),E.ncmp,2);
mu1CI=zeros(length(D.idxStd),E.ncmp,2);
mu2CI=zeros(length(D.idxStd),E.ncmp,2);
cr1CI=zeros(length(D.idxStd),E.ncmp,2);
cr2CI=zeros(length(D.idxStd),E.ncmp,2);

%COLORS
colorInc=floor(length(colormap(spring))/E.ncmp);
coloridx=1:colorInc:E.ncmp*colorInc;
colors=colormap(spring);
colors=colors(coloridx,:);
%Purple to red
%colors(:,3)=.5;
%colors(:,1)=colors(:,2);
%colors(:,2)=0;
%BLACK TO GREY
colors(:,2)=colors(:,2)*.9;
colors(:,1)=colors(:,2);
colors(:,3)=colors(:,2);

PA=zeros(length(D.idxStd),E.ncmp);
PC=zeros(length(D.idxStd),E.ncmp);
figa=figure(nFn);
if bSave==1
	set(figa, 'Visible', 'off');
end
set(gcf,'Position',[10 10 640 300+270*tSubplot/2]);
try
    parpool;
end
tic

%D.idxStd
Data=zeros(length(D.idxStd),E.ncmp,E.trlPerLvl,length(prjCodes));
for s = 1:length(D.idxStd)
    cSubplot=cSubplot+1;
    data=squeeze(DATA(s,:,:));
	cmpXunq=unique(cmpXsrt(s,:,1));
    for c = 1:E.ncmp
        Dtmp=data(cmpXsrt(s,:,1)==cmpXunq(c),:); %trials comparison chosen
        Data(s,c,:,:)=Dtmp;
        %size(Data)

        if D.muFrmDp == 0
            if D.bCI==1
                [rhoCI(s,c,:),mu1CI(s,c,:),mu2CI(s,c,:),cr1CI(s,c,:),cr2CI(s,c,:),rhoDSTB,mu1DSTB,mu2DSTB] = psyFitDecisionVariableCorrBootstrap(D.modelType,Dtmp,D.rhoFix,D.mu1Fix,D.mu2Fix,D.cr1Fix,D.cr2Fix,CI,nBoot,bPLOT,cSubplot,tSubplot,colors(c,:),1);
                rho(s,c)=mean(rhoDSTB);
                mu1(s,c)=mean(mu1DSTB);
                mu2(s,c)=mean(mu2DSTB);
            else
                [rho(s,c),mu1(s,c),mu2(s,c)]=psyFitDecisionVariableCorr(D.modelType,Dtmp,D.rhoFix,D.mu1Fix,D.mu2Fix,D.cr1Fix,D.cr2Fix,bPLOT,bBOOTSTRAPPING,cSubplot,tSubplot,colors(c,:),1);
            end
        elseif D.muFrmDp == 1
            if D.bCI==1
                [rhoCI(s,c,:),mu1CI(s,c,:),mu2CI(s,c,:),cr1CI(s,c,:),cr2CI(s,c,:),rhoDSTB,mu1DSTB,mu2DSTB] = psyFitDecisionVariableCorrBootstrap(D.modelType,Dtmp,D.rhoFix,DP(s,c)/2,DP(s,c)/2,D.cr1Fix,D.cr2Fix,CI,nBoot,bPLOT,cSubplot,tSubplot,colors(c,:),1);
                rho(s,c)=mean(rhoDSTB);
                mu1(s,c)=mean(mu1DSTB);
                mu2(s,c)=mean(mu2DSTB);
            else
                [rho(s,c),mu1(s,c),mu2(s,c)]=psyFitDecisionVariableCorr(D.modelType,Dtmp,D.rhoFix,DP(s,c)/2,DP(s,c)/2,D.cr1Fix,D.cr2Fix,bPLOT,bBOOTSTRAPPING,cSubplot,tSubplot,colors(c,:),1);
            end
        elseif D.muFrmDp == 2
            if D.bCI==1
                [rhoCI(s,c,:),mu1CI(s,c,:),mu2CI(s,c,:),cr1CI(s,c,:),cr2CI(s,c,:),rhoDSTB,mu1DSTB,mu2DSTB] = psyFitDecisionVariableCorrBootstrap(D.modelType,Dtmp,D.rhoFix,DP(1,s,c)/2,DP(1,s,c)/2,D.cr1Fix,D.cr2Fix,CI,nBoot,bPLOT,cSubplot,tSubplot,colors(c,:),1);
                rho(s,c)=mean(rhoDSTB);
                mu1(s,c)=mean(mu1DSTB);
                mu2(s,c)=mean(mu2DSTB);
            else
                [rho(s,c),mu1(s,c),mu2(s,c)]=psyFitDecisionVariableCorr(D.modelType,Dtmp,D.rhoFix,DP(1,s,c)/2,DP(2,s,c)/2,D.cr1Fix,D.cr2Fix,bPLOT,bBOOTSTRAPPING,cSubplot,tSubplot,colors(c,:),1);
            end
        end

        if s == 1 && c == 1
            %disp([ 'Current time ' num2str(datetime(now)) ]);
            disp([ 'Estimated time remaining: ' num2str((toc*length(D.idxStd)*E.ncmp-toc)/60) ' min.' ]);
            progressreport(c+(s-1)*E.ncmp,1,length(D.idxStd)*E.ncmp)
        else
            progressreport(c+(s-1)*E.ncmp,5,length(D.idxStd)*E.ncmp)
        end
    end

    %plot binomial and data point
    hold off
    cSubplot=cSubplot+1;
    [PA(s,:),PC(s,:)]=psyNpassMagrEmpE(SS{s},length(prjCodes),1:length(prjCodes),.68*100,1,cSubplot,tSubplot,1);

    %fit data curve
    hold on
    med=median(rho(s,:));
    psyFitDecisionVariableCorrRspAgree(med,mu1(s,:),mu2(s,:),D.cr1Fix,D.cr2Fix,1,'-',1,1);%,symType)
    hold off

    if s == 1
        %disp([ 'Current time ' num2str(datetime(now)) ]);
        disp(['Estimated time remaining: ' num2str((toc*length(D.idxStd)-s*toc)/60) ' min.']);
    end
end

y=zeros(length(D.stdXo),E.ncmp);
yM=zeros(length(D.stdXo));
for s = 1:length(D.stdXo)
    for c = 1:E.ncmp
        rhoTmp=abs(real(rho(s,c)));
        rhoTmp2=abs(imag(rho(s,c)));
        if rhoTmp<rhoTmp2
            rhoTmp=rhoTmp2;
        end
        if rhoTmp==0
            y(s,c)=10;
        else
            y(s,c)=sqrt((1-rhoTmp)./(rhoTmp));
        end
    end
    yM(s,c)=geomean(y(s,:));
end


% -------------------------------------------------------------------------------

if bSave==1
    disp('TAGGED SAVING')
    cd(D.FigureFolder)
    LASTD=pwd;
end

%REMOVE ALL Xlabels other than last 2
for i = 1:tSubplot-2
    subplot(tSubplot/2,2,i);
    hold on
    xlabel('');
end

[eid,idx]=sort(expIDs);
eid=eid{1};
pjc=prjCodes(idx);
pjc=pjc{1};
titl=[pjc ' ' strrep(eid,'_','-') ' ' subj];

if D.muFrmDp==1
    strDP='-fDPtogether';
elseif D.muFrmDp==2
    strDP='-fDPseperate';
else
    strDP='';
end
if D.bSwap==1
    strSWP='-SWAP';
else
    strSWP='';
end

sutitle=[titl newline D.modelType strDP strSWP];
savetitl=[titl '__' D.modelType strDP strSWP];
savetitl=strrep(savetitl,' ','_');
suptitle(sutitle)

if bSave==1
    saveas(gcf,[savetitl '_' 'Agree'],'png')
    saveas(gcf,[savetitl '_' 'Agree'],'epsc')
    %saveas(gcf,[savetitl '_' 'Agree'],'fig')
	savefig(figa,[savetitl '_' 'Agree'])
end

disp([ num2str(toc/60) ' minutes elapsed']);

rMean=mean(rho(:));
yMean=geomean(y(:));
% -------------------------------------------------------------------------------
% PLOT RHO V. %comparison
figr=figure(nFn);
if bSave==1
	set(figr, 'Visible', 'off');
end
set(gcf,'Position',[10 10 640 600+400*tSubplot/2]);
for s = 1:length(D.idxStd)
    mea=mean(rho(s,:));
    med=median(rho(s,:));
    subplot(length(D.idxStd),1,s);
    hold on
    p1=plot([0 100],[mea mea],'k-');
    p2=plot([0 100],[med med],'k:');
    if D.bCI==1
        errorbar(PC(s,:),rho(s,:),rhoCI(s,:,1),rhoCI(s,:,2),'ko')
    else
        p3=scatter(100*PC(s,:),rho(s,:),'k','LineWidth',1);
    end
    legend([p1 p2],['\mu(\rho) ' num2str(round(mea,2))],['M(\rho) ' num2str(round(med,2))],'Location','southeast');
    if s==length(D.idxStd)
        formatFigure('% Comparison Chosen','\rho',D.stdXo(s));
    else
        formatFigure('','\rho',D.stdXo(s));
    end
    ylim([-.5,1]);
    axis square
    hold off
end
suptitle(sutitle)
if bSave==1
    saveas(gcf,[savetitl '_' 'Rho'],'png')
    saveas(gcf,[savetitl '_' 'Rho'],'epsc')
    %saveas(gcf,[savetitl '__' D.modelType strSWP strDP '_' 'Rho'],'fig')
	savefig(figr,[savetitl '_' 'Rho'])
end
% -------------------------------------------------------------------------------
%PLOT MU's
%XXX if D.muFrmDp==2 || D.mu1Fix~=D.mu2Fix || mu1~=mu2
    figm=figure(nFn);
	if bSave==1
		set(figm, 'Visible', 'off');
	end
    set(gcf,'Position',[10 10 400 50+400*tSubplot/2]);
    for s = 1:length(D.idxStd)
        subplot(length(D.idxStd),1,s);
        hold on
        plot([-100 100],[-100 100],'k:')
        scatter(mu1(s,:),mu2(s,:),'k','LineWidth',1);
        if s==length(D.idxStd)
            formatFigure(['\mu' num2str(char(8321))],['\mu' num2str(char(8322))],D.stdXo(s));
        else
            formatFigure('',['\mu' num2str(char(8322))],D.stdXo(s));
        end
        ylim([min(min(mu2))-.5,max(max(mu2))+.5]);
        xlim([min(min(mu1))-.5,max(max(mu1))+.5]);
        axis square
        hold off
    end
    suptitle(sutitle)
    if bSave==1
        saveas(gcf,[savetitl '_' 'Mu'],'png')
        saveas(gcf,[savetitl '_' 'Mu'],'epsc')
        %saveas(gcf,[savetitl '_' 'Mu'],'fig')
		savefig(figm,[savetitl '_' 'Mu'])
    end
%end
% -------------------------------------------------------------------------------
%HIST RHO
figh=figure(nFn);
if bSave==1
    set(figh, 'Visible', 'off');
end
binsz=10;
histogram(rho(:),binsz)
set(gca, 'Xscale', 'linear')
formatFigure('\rho','count',[sutitle newline '\mu(\rho) = ' num2str(mean(rho(:))) newline 'M(\rho) = ' num2str(median(rho(:))) ])
if bSave==1
    saveas(gcf,[savetitl '_' 'HistRho'],'png')
    saveas(gcf,[savetitl '_' 'HistRho'],'epsc')
    saveas(gcf,[savetitl '_' 'HistRho'],'svg')
    savefig(figh,[savetitl '_' 'HistRho'])
end
% -------------------------------------------------------------------------------
%HIST KAPPA
figh=figure(nFn);
if bSave==1
    set(figh, 'Visible', 'off');
end
binsz=10;
histogram(y(:),binsz)
[Xhst,~,~,XbinEnds]=histogramPlot(y(:),binsz,1,0);
XbinEnds=[XbinEnds(:,1); XbinEnds(end,2)];
set(gca, 'Xscale', 'log')
histogram('BinEdges',XbinEnds','BinCounts',Xhst),set(gca, 'Xscale', 'log')
formatFigure('\kappa','count',[sutitle newline '\mu(\kappa) = ' num2str(geomean(y(:))) newline 'M(\kappa) = ' num2str(median(y(:))) ])
if bSave==1
    saveas(gcf,[savetitl '_' 'HistKappa'],'png')
    saveas(gcf,[savetitl '_' 'HistKappa'],'epsc')
    saveas(gcf,[savetitl '_' 'HistKappa'],'svg')
    savefig(figh,[savetitl '_' 'HistKappa'])
end
% -------------------------------------------------------------------------------
%PLOT KAPPA
figk=figure(nFn);
if bSave==1
    set(figk, 'Visible', 'off');
end
plot(D.stdXo,yM,'k')
hold on
for s = 1:length(D.idxStd)
    scatter(D.stdXo(s)*ones(E.ncmp,1),y(s,:),'k');
end
formatFigure('','\kappa','')
set(gca, 'Yscale', 'log')
set(gca, 'Xscale', 'linear')
ylim([.1, 10])
xlim([-12,-3.5])
yticks([0 .1 .3 1 3 10])
yticklabels({'0', '.1', '.3', '1', '3','10'})
xlabel('disparity (arcmin)')
axis square
title(sutitle)
hold off
if bSave==1
    saveas(gcf,[savetitl '_' 'Kappa'],'png')
    saveas(gcf,[savetitl '_' 'Kappa'],'epsc')
    savefig(figk,[savetitl '_' 'Kappa'])
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%if bSave==1 && strcmp(info.localHostName,'BONICN')
%    cd LASTD;
%end
ncmp=E.ncmp;
