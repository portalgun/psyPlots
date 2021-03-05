function [] = expDoublePassPopulation(signal,D,prjCode,subj,DP,rho,y,yM,men,ncmp,sutitle,sI,sE,rhoMenSub,yMenSub,mCI)

for W=1:length(prjCode)
    if W==1
        ym=cell(length(prjCode),length(D{1}.idxStd));
    end

    %Kapp ALL
    figK=figure(nFn);
    if D{W}.bSave==1
        set(figK, 'Visible', 'off');
    end

    %STAGGER SUBJECT DATA POINTS ON X AXIS BY FRACTION OF SEPERATION BETWEEN POINTS
    stagger=floor(length(subj)/2);
    staggerR=((-1*stagger):1:stagger)*(D{1}.stdXo(2)-D{1}.stdXo(1))*.2;

    for Z = 1:length(subj)
    for s = 1:length(D{1}.idxStd)
        shape=markerDefs(subj{Z});
        scatter(D{W}.stdXo(s)*ones(ncmp,1)+staggerR(Z),y{W,Z}(s,:),['k' shape],'MarkerFaceColor',[1 1 1]);
        hold on
    end
    end

    %plot(D{W}.stdXo,med(W,:),'k')
    plot(D{W}.stdXo,men(W,:),'k--')
    formatFigure('','\kappa','')
    set(gca, 'Yscale', 'log')
    set(gca, 'Xscale', 'linear')
    ylim([.1, 10])
    xlim([-12.2,-2.75])
    yticks([0 .1 .3 1 3 10])
    yticklabels({'0', '.1', '.3', '1', '3','10'})
    xlabel('disparity (arcmin)')
    title(strrep(sutitle,subj{end},'ALL'));
    savetitl=strrep(sutitle,newline,'__');
    savetitl=strrep(savetitl,' ','_');
    hold off
    if D{W}.bSave==1
        saveas(gcf,[savetitl '_' 'KappaALL'],'png')
        saveas(gcf,[savetitl '_' 'KappaALL'],'epsc')
        %saveas(gcf,[savetitl '_' 'Agree'],'fig')
        savefig(figK,[savetitl '_' 'KappaSepALL'])
    end


    figM=figure(nFn);
    set(figM, 'Visible', 'on');
    scatter(D{W}.stdXo,men(W,:),'k')
    hold on
    plotfillederror(D{W}.stdXo,mCI(W,:,1),mCI(W,:,2),'k',0); hold on;
    hold on
    set(gca, 'Yscale', 'log')
    set(gca, 'Xscale', 'linear')
    ylim([.1, 10])
    xlim([-12.2,-2.75])
    yticks([0 .1 .3 1 3 10])
    yticklabels({'0', '.1', '.3', '1', '3','10'})
    saveas(gcf,[savetitl 'MENS'],'png');
    saveas(gcf,[savetitl 'MENS'],'svg');
    savefig(figM,[savetitl 'MENS']);

    %Kapp ALL
    figKN=figure(nFn);
    if D{W}.bSave==1
        set(figKN, 'Visible', 'off');
    end

    %STAGGER SUBJECT DATA POINTS ON X AXIS BY FRACTION OF SEPERATION BETWEEN POINTS
    %YE=zeros(length(subj),length(D{W}.idxStd),length(y{W,1}(s,:)));
    %YI=zeros(length(subj),length(D{W}.idxStd),length(y{W,1}(s,:)));
    for Z = 1:length(subj)
    for s = 1:length(D{W}.idxStd)
        shape=markerDefs(subj{Z});

        %Sig=squeeze(mean(signal{W,Z}(:,s,:)))';
        %dP=abs(squeeze(mean(     DP{W,Z}(:,s,:))))';
        %R=real(y{W,Z}(s,:)); %XXX
        %YE(Z,s,:)=Sig./(dP.*(R+1));
        %YI(Z,s,:)=squeeze(YE(Z,s,:)).*Sig';
        X=D{W}.stdXo(s)*ones(ncmp,1)+staggerR(Z);
        scatter(X,sI(W,s,Z,:),['k' shape],'MarkerFaceColor',[1 1 1]);
        hold on
        scatter(X,sE(W,s,Z,:),['r' shape],'MarkerFaceColor',[1 1 1]);
        hold on
    end
    end

    %YE(isnan(YE))=0;
    %YI(isnan(YI))=0;

    %menE=zeros(length(prjCode),length(D{W}.idxStd));
    %menI=zeros(length(prjCode),length(D{W}.idxStd));
    %for s = 1:length(D{1}.idxStd)
    %    menI(W,s)=mean(squeeze(mean(YI(:,s,:),1)));
    %    menE(W,s)=mean(squeeze(mean(YE(:,s,:),1)));
    %end

    % -------------------------------------------------------------------------------
    %Plot KAPA norm
    %plot(D{W}.stdXo,menI(W,:),'k--')
    %plot(D{W}.stdXo,menE(W,:),'r--')
    formatFigure('','\sigma','')
    set(gca, 'Yscale', 'linear')
    set(gca, 'Xscale', 'linear')
    %ylim([-.6, .3])
    xlim([-12.2,-2.75])
    %yticks([0 .1 .3 1 3 10])
    %yticklabels({'0', '.1', '.3', '1', '3','10'})
    xlabel('disparity (arcmin)')
    title(strrep(sutitle,subj{end},'ALL'));
    savetitl=strrep(sutitle,newline,'__');
    savetitl=strrep(savetitl,' ','_');
    hold off
    if D{W}.bSave==1
        saveas(gcf,[savetitl '_' 'NoiseSep'],'png')
        saveas(gcf,[savetitl '_' 'NoiseSep'],'epsc')
        %saveas(gcf,[savetitl '_' 'Agree'],'fig')
        savefig(figKN,[savetitl '_' 'NoiseSep'])
    end


    % -------------------------------------------------------------------------------
    %HISTOGRAM RHO
    figH=figure(nFn);
    if D{1}.bSave==1
        set(figH, 'Visible', 'off');
    end
    binsz=13;
    RHO=[];
    for Z=1:length(subj)
        RHO=[RHO; rho{W,Z}];
    end
    histogram(RHO(:),binsz)
    set(gca, 'Xscale', 'linear')
    xlim([-1,1])
    ylim([0,40])
    for Z = 1:length(subj)
        hold on
        shape=markerDefs(subj{Z});
        plot(rhoMenSub{W,Z},40,['k' shape]);
    end
	axis square
    formatFigure('\rho','count',strrep(sutitle,subj{end},'ALL'))
    if D{1}.bSave==1
        saveas(gcf,[savetitl '_' 'HistRhoALL'],'png')
        saveas(gcf,[savetitl '_' 'HistRhoALL'],'epsc')
        saveas(gcf,[savetitl '_' 'HistRhoALL'],'svg')
        %saveas(gcf,[savetitl '_' 'HistRho'],'fig')
        savefig(figH,[savetitl '_' 'HistRhoALL'])
    end
    hold off

    % -------------------------------------------------------------------------------
    %HISTOGRAM Kappa
    figH=figure(nFn);
    if D{1}.bSave==1
        set(figH, 'Visible', 'off');
    end
    binsz=12;
    K=[];
    for Z=1:length(subj)
        K=[K; y{W,Z}];
    end
    [Xhst,~,~,XbinEnds]=histogramPlot(K(:),binsz,1,0);
    XbinEnds=[XbinEnds(:,1); XbinEnds(end,2)];
    set(gca, 'Xscale', 'log')
    histogram('BinEdges',XbinEnds','BinCounts',Xhst),set(gca, 'Xscale', 'log')
    xlim([.1,10])
    ylim([0,40])
    for Z = 1:length(subj)
        hold on
        shape=markerDefs(subj{Z});
        plot(yMenSub{W,Z},40,['k' shape]);
    end
	axis square
    %ylim([0,max(Xhst)+.1*max(Xhst)])
    formatFigure('\kappa','count',strrep(sutitle,subj{end},'ALL'))
    if D{1}.bSave==1
        saveas(gcf,[savetitl '_' 'HistKappALL'],'png')
        saveas(gcf,[savetitl '_' 'HistKappALL'],'epsc')
        saveas(gcf,[savetitl '_' 'HistKappALL'],'svg')
        %saveas(gcf,[savetitl '_' 'HistRho'],'fig')
        savefig(figH,[savetitl '_' 'HistKappaALL'])
    end
    hold off

end
%[DP{W,Z},rho{W,Z},y{W,Z},yM{W,Z},ncmp,sutitle]=expDoublePass(D{W},prjcodes{W,:},expids{W,:},subj{Z},D{W}.bSave);

