classdef psyCurve < handle
properties
    %colors=[[18,133,158];[49,97,107];[45,209,148];[215,103,100];[158,18,73];[18,133,158]];
    %lcolors=colors ./ 255;
    mFix
    sFix
    bFix
    DPcrt
    nIntrvl
    nBoot
    CIsz
    prcntUse
    flipCmpChs

    fitOpts
    shape
    minFuncType

    nTrials
    nLvls
    nTrialsPerLvl

    nX % curve resolution
    lcolor
    color % line color
    CIalpha
    markersize
    LineWidth
    markerface
    Xname
    Xunits
    bPlotCI
    bExtendCurve
    bPlot
    bTitle

    % Fit
    mMU
    sMU
    bMU
    tMU

    mSE
    sSE
    bSE
    tSE

    mCI
    sCI
    bCI
    tCI

    mFit
    sFit
    bFit
    tFit
    PCfit
    DPfit

    mFitDstb
    sFitDstb
    bFitDstb
    tFitDstb

    Ymin
    Ymax
    Xmin
    Xmax
end
properties(Hidden=true)
    bP
    bTest
    i
   
    stdX
    cmpX
    RcmpChs

    stdXsel
    cmpXsel
    RcmpChsSel

    mFitSel
    bFitSel
    sFitSel
    tFitSel
    PCfitSel
    DPfitSel
end
methods
    function obj=psyCurve(stdX,cmpX,RcmpChs,Opts)
        bTest=0;
        if ~exist('stdX','var')
            bTest=1;
        else
            obj.stdX=stdX;
        end
        if ~exist('cmpX','var')
            bTest=1;
        else
            obj.cmpX=cmpX;
        end
        if ~exist('RcmpChs','var')
            bTest=1;
        else
            obj.RcmpChs=RcmpChs;
        end

        % XXX Make sure just 1 standard
        % XXX make sure std cmp Rchs are same size
        if ~exist('Opts','var')
            Opts=struct;
        end
        obj.parse_opts(Opts);


        if bTest
            obj.bTest=1;
        end

        if obj.bTest
            obj.gen_testX()
        end

        if numel(obj.stdX) == 1
            obj.stdX = obj.stdX*ones(numel(obj.cmpX),1);
        end

        % sizes
        if size(obj.stdX,2)    ~= 1
            obj.stdX    = obj.stdX(:);
        end
        if size(obj.cmpX,2)    ~= 1
            obj.cmpX    = obj.cmpX(:);
        end
        if size(obj.RcmpChs,2) ~= 1
            obj.RcmpChs = obj.RcmpChs(:);
        end
        %if obj.flipCmpChs
        %    obj.RcmpChs=~obj.RcmpChs;
        %end

        obj.nTrials=length(obj.stdX);
        obj.nLvls=numel(unique(obj.cmpX));
        obj.nTrialsPerLvl=obj.nTrials/obj.nLvls;

        % SET FMINCON OPTIONS
        %obj.minFuncType = 'fmincon'; % XXX
        if strcmp(obj.minFuncType,'fmincon')
            obj.fitOpts             = optimset('fmincon');
            obj.fitOpts.Algorithm   = 'active-set';
            obj.fitOpts.LargeScale  = 'off';
            obj.fitOpts.UseParallel = 'never';
            obj.fitOpts.Display     = 'none';
            obj.fitOpts.MaxIter     = 50;
        elseif strcmp(obj.minFuncType,'fminsearch')
            obj.fitOpts             = optimset('fminsearch');
            obj.fitOpts.UseParallel = 'never';
            obj.fitOpts.Display     = 'off';
            obj.fitOpts.MaxIter     = 50;
        end

        obj.run();
    end
    function obj=parse_opts(obj,Opts)
        if ~exist('Opts','var')
            Opts=struct;
        end
        p=inputParser();
        p.addParameter('mFix',[]);
        p.addParameter('sFix',[]);
        p.addParameter('bFix',1);
        p.addParameter('DPcrt',1); % XXX 1.36?
        p.addParameter('nIntrvl',2);
        p.addParameter('nBoot',1000);
        p.addParameter('CIsz',95);
        p.addParameter('prcntUse',100);
        p.addParameter('minFuncType','fmincon');
        p.addParameter('bTest',0);
        p.addParameter('bPlot',0);
        p.addParameter('bPlotCI',1);
        p.addParameter('LineWidth',2);
        p.addParameter('markersize',12);
        p.addParameter('markerface','w');
        p.addParameter('color','k');
        p.addParameter('lcolor','k');
        p.addParameter('CIalpha',.4);
        p.addParameter('shape','o');
        p.addParameter('Xname','X');
        p.addParameter('Xunits','');
        p.addParameter('nX',200);
        p.addParameter('bExtendCurve',0);
        p.addParameter('bTitle',1);
        p.addParameter('flipCmpChs',0)

        p=parseStruct(Opts,p);
        flds=fieldnames(p.Results);

        for i = 1:length(flds)
            fld=flds{i};
            obj.(fld)=p.Results.(fld);
        end
    end
    function obj=run(obj)
        obj.fit_boot();
        obj.fit_basic();
        if obj.bPlot==1
            obj.plot();
        end
    end
%% MAIN
    function obj=fit_boot(obj)
        obj.mFitDstb=zeros(obj.nBoot,1);
        obj.sFitDstb=zeros(obj.nBoot,1);
        obj.bFitDstb=zeros(obj.nBoot,1);
        obj.tFitDstb=zeros(obj.nBoot,1);
        for i = 1:obj.nBoot
            obj.i=i;

            obj.select_boot();
            obj.fit();
            obj.gen_gauss_sel(); %get T and PC & DP
            obj.pack_boot_ind();
        end
        obj.pack_boot();
    end
    function obj=fit_basic(obj)
    %NON BOOTSTRAPPED FIT
        obj.select_all();
        obj.fit();
        obj.PC_fit();
        obj.pack_all();
    end
%% SELECT
    function obj=select_all(obj)
        % SAMPLE ALL, NO REPLACMENT (AS OPPOSED TO BOOT)
        ind=~isnan(obj.stdX) & ~isnan(obj.cmpX) & ~isnan(obj.RcmpChs);
        obj.stdXsel=obj.stdX(ind);
        obj.cmpXsel=obj.cmpX(ind);
        obj.RcmpChsSel=obj.RcmpChs(ind);
    end
    function obj=select_boot(obj)
        ind=~isnan(obj.stdX) & ~isnan(obj.cmpX) & ~isnan(obj.RcmpChs);
        [~,indSmp] = datasample(obj.RcmpChs(ind),round(numel(obj.RcmpChs(ind)).*obj.prcntUse./100));
        % REFIT
        obj.stdXsel=obj.stdX(indSmp);
        obj.cmpXsel=obj.cmpX(indSmp);
        obj.RcmpChsSel=obj.RcmpChs(indSmp);
    end
%% PACK
    function obj=pack_boot_ind(obj)
        obj.mFitDstb(obj.i,:)=obj.mFitSel;
        obj.sFitDstb(obj.i,:)=obj.sFitSel;
        obj.bFitDstb(obj.i,:)=obj.bFitSel;
        obj.tFitDstb(obj.i,:)=obj.tFitSel;
    end
    function obj=pack_boot(obj)
        % CONFIDENCE INTERVAL: LO & HI BOUNDS
        CIlohi = 0.5*(1-obj.CIsz/100) + [0 obj.CIsz/100];
        % BOOSTRAPPED CONFIDENCE INTERVALS
        obj.mCI = quantile(obj.mFitDstb(~isnan(obj.mFitDstb)), CIlohi);
        obj.sCI = quantile(obj.sFitDstb(~isnan(obj.sFitDstb)), CIlohi);
        obj.bCI = quantile(obj.bFitDstb(~isnan(obj.bFitDstb)), CIlohi);
        obj.tCI = quantile(obj.tFitDstb(~isnan(obj.tFitDstb)), CIlohi);
        % BOOSTRAPPED STD ERR OF STATISTIC
        obj.mSE   = std(obj.mFitDstb(~isnan(obj.mFitDstb)));
        obj.sSE   = std(obj.sFitDstb(~isnan(obj.sFitDstb)));
        obj.bSE   = std(obj.bFitDstb(~isnan(obj.bFitDstb)));
        obj.tSE   = std(obj.tFitDstb(~isnan(obj.tFitDstb)));
        % BOOSTRAPPED MEAN
        obj.mMU = mean(obj.mFitDstb(~isnan(obj.mFitDstb)));
        obj.sMU = mean(obj.sFitDstb(~isnan(obj.sFitDstb)));
        obj.bMU = mean(obj.bFitDstb(~isnan(obj.bFitDstb)));
        obj.tMU = mean(obj.tFitDstb(~isnan(obj.tFitDstb)));
    end
    function obj=pack_all(obj)
        obj.mFit=obj.mFitSel;
        obj.sFit=obj.sFitSel;
        obj.bFit=obj.bFitSel;

        obj.tFit=obj.tFitSel;
        obj.PCfit=obj.PCfitSel;
        obj.DPfit=obj.DPfitSel;
    end
%% FIT
    function obj = fit(obj)
        pLB     = [min(obj.cmpXsel(:))*1.2 0.01.*(max(obj.cmpXsel(:))-min(obj.cmpXsel(:))) 0.35];
        pUB     = [max(obj.cmpXsel(:))*1.2 10.0.*(max(obj.cmpXsel(:))-min(obj.cmpXsel(:))) 3.00];

        % SET INITIAL PARAMETER VALUES
        m0  = obj.mFix;
        s0  = obj.sFix;
        b0  = obj.bFix;
        if isempty(m0)
            m0 = mean([min(obj.cmpXsel(:)) max(obj.cmpXsel(:))]);
            m0 = m0  + .1.*randn;
        end
        if isempty(s0)
            v=abs(obj.cmpXsel);
            mm=[min(v) max(v)];
            s0 = diff(mm)./6;
            s0 = s0  + .1.*s0.*randn;
        end
        if isempty(b0)
            b0 = 1;
            b0 = b0  + .1.*b0.*randn;
        end

        p0 = [m0 s0 b0];


        negLL = @(p) ...
              -(sum(log(     obj.gen_gauss(obj.cmpXsel(obj.RcmpChsSel==1),p(1),p(2),p(3)) )) + ...
                sum(log( 1 - obj.gen_gauss(obj.cmpXsel(obj.RcmpChsSel==0),p(1),p(2),p(3)) )) );

        % MINIMIZE NEGATIVE LOG-LIKELIHOOD
        switch obj.minFuncType
        case 'fmincon'
            [pFit,negLL] = fmincon(negLL,p0,[],[],[],[],pLB,pUB,[],obj.fitOpts);
        case 'fminsearch'
            [pFit,negLL] = fminsearch(neg_LL,p0,obj.fitOpts);
        end

        % FINAL FIT PARAMETERS
        if isempty(obj.mFix); obj.mFitSel = pFit(1); else; obj.mFitSel = obj.mFix; end
        if isempty(obj.sFix); obj.sFitSel = pFit(2); else; obj.sFitSel = obj.sFix; end
        if isempty(obj.bFix); obj.bFitSel = pFit(3); else; obj.bFitSel = obj.bFix; end
    end
    function obj=gen_gauss_sel(obj)
        [obj.PCfitSel,obj.tFitSel,obj.DPfitSel] = obj.gen_gauss(obj.cmpXsel,obj.mFitSel,obj.sFitSel,obj.bFitSel);
    end
    function obj=PC_fit(obj)
        cmpXUnq = unique(obj.cmpXsel);
        [obj.PCfitSel,obj.tFitSel,obj.DPfitSel] = obj.gen_gauss(cmpXUnq,obj.mFitSel,obj.sFitSel,obj.bFitSel);
    end
    function [PC,T,DP] = gen_gauss(obj,cmpX,mFit,sFit,bFit)
        DP = sign(cmpX-mFit).*(abs(cmpX-mFit)./sFit).^bFit; % abs() prevent complex numbers w. some betas
        PC = normcdf( 0.5.*sqrt(obj.nIntrvl).*DP,0,1); % sign() reinstates sign of Xcmp-mFit
        T  = sFit.*obj.DPcrt.^(1/bFit);
    end
%%
    function [PC,N1,N,N0]=get_PC(obj)
        % STANDARD VALUES
        XstdUnq = unique(obj.stdXsel);

        % LOOP OVER STANDARDS
        for s = 1:length(XstdUnq)
            % INDICES FOR EACH STANDARD
            indS = obj.stdXsel == XstdUnq(s);
            % COMPARISON VALUES
            XcmpUnq(:,s) = unique(obj.cmpXsel(indS));
            % LOOP OVER COMPARISONS
            for c = 1:length(XcmpUnq(:,s))
                % INDICES IN STD / CMP CONDITION
                indCnd  = obj.stdXsel ==XstdUnq(s) & obj.cmpXsel==XcmpUnq(c,s);
                % TOTAL NUMBER OF TRIALS IN CONDITION
                N(c,s)  = sum( indCnd );
                % TOTAL NUMBER OF CMP CHOSEN IN CONDITION
                N1(c,s) = sum( obj.RcmpChsSel(indCnd)==1 );
                % TOTAL NUMBER OF STD CHOSEN IN CONDITION
                N0(c,s) = sum( obj.RcmpChsSel(indCnd)==0 );
            end
        end

        PC = N1./N;
    end
%% PLOT
    function [] = plot_all(obj)
        figure(nFn);
        hold off
        subPlot([2,2],1,1)
        obj.plot_boot_curve();
        obj.plot_boot_params(2);
        subPlot([2 2],2,1)
        obj.plot_boot_DP();
    end
    function [] = plot(obj)
        figure(nFn)
        hold off
        obj.plot_boot_curve();
    end
    function [] = plot_params(obj)
        figure(nFn)
        hold off
        obj.plot_boot_params();
    end
    function [] = plot_DP(obj)
        figure(nFn)
        hold off
        obj.plot_boot_DP();
    end
%%%
end
methods(Hidden = true)
    function [] = plot_boot_params(obj,c)
        if ~exist('c','var') || isempty(c)
            c=1;
        end

        obj.get_bP();
        sSz=[sum(obj.bP)+1,c];

        subPlot(sSz,1,c);
        obj.plot_boot_T_p();

        i=1;
        if obj.bP(i)
            subPlot(sSz,i+1,c);
            obj.plot_boot_mu_p();
        end

        i=2;
        if obj.bP(i)
            subPlot(sSz,i+1,c);
            obj.plot_boot_sigma_p();
        end

        i=3;
        if obj.bP(i)
            subPlot(sSz,i+1,c);
            obj.plot_boot_bet_p();
        end


    end

    function obj = plot_boot_curve(obj,sym)
        if ~exist('sym','var') || isempty(sym)
            sym=obj.shape;
        end

        % UNIQUE COMPARISON VALUES
        Xdat = transpose(unique(obj.cmpX));
        % PROPORTION CMP CHOSEN
        for i = 1:length(Xdat)
            ind = obj.cmpX(:) == Xdat(i);
            Ydat(i) = mean(obj.RcmpChs(ind));
        end

        a=min(obj.cmpX);
        b=max(obj.cmpX);
        if obj.bExtendCurve
            rang=b-a;
            a=a-rang*2;
            b=b+rang*2;
        end
        if isempty(obj.nX)
            obj.nX=200;
        end

        Xfit = linspace(a,b,obj.nX);

        % PSYCHOMETRIC FUNCTION FIT MEAN
        Yfit = obj.gen_gauss(Xfit,obj.mFit,obj.sFit,obj.bFit);
        if  obj.bExtendCurve
            ind=Yfit < .999 & Yfit > .001;
            Yfit=Yfit(ind);
            Xfit=Xfit(ind);
        end

        obj.get_bP();

        if obj.bPlotCI
            if obj.bP(1)
                mlh(1)=obj.mCI(1);
                mlh(2)=obj.mCI(2);
            else
                mlh(1)=obj.mFix;
                mlh(2)=obj.mFix;
            end
            if obj.bP(2)
                slh(1)=obj.sCI(1);
                slh(2)=obj.sCI(2);
            else
                slh(1)=obj.sFix;
                slh(2)=obj.sFix;
            end
            if obj.bP(2)
                blh(1)=obj.bCI(1);
                blh(2)=obj.bCI(2);
            else
                blh(1)=obj.bFix;
                blh(2)=obj.bFix;
            end

            c=distribute(mlh,slh,blh);
            c=unique(c,'rows');
            Y=zeros(size(c,2),length(Xfit));
            for i = 1:size(c,1)
                ind=c(i,:);
                Y(i,:)  = obj.gen_gauss(Xfit,ind(1),ind(2),ind(3));
            end
            Ymin=min(Y,[],1);
            Ymax=max(Y,[],1);
            obj.Ymin=min(Yfit);
            obj.Ymax=max(Yfit);
            obj.Xmin=min(Xfit);
            obj.Xmax=max(Xfit);

            %plot(Xfit,Ymin,'r'); hold on
            %plot(Xfit,Ymax,'r');
            %HERE
            patch([Xfit fliplr(Xfit)], [Ymin fliplr(Ymax)], obj.color,'FaceAlpha',obj.CIalpha,'EdgeColor','none'); hold on
        end

        plot(Xfit,Yfit,'color',obj.lcolor,'LineWidth',obj.LineWidth); hold on;


        % Data
        plot(Xdat,Ydat,sym,'color',obj.lcolor, 'markerface',obj.markerface,'markersize',obj.markersize,'LineWidth',obj.LineWidth);
    end
    function obj=format_boot_curve(obj)

        if obj.bTitle
            titl=[

                '\mu='    num2str(obj.mFit,'%2.2f') ...
                ', \sigma=' num2str(obj.sFit,'%2.2f') ...
                ', \beta='  num2str(obj.bFit,'%2.2f') ...
            ];
        else
            titl=[];
        end
        formatFigure(obj.Xname,'Proportion Cmp Chosen',titl)
        axis square
    end
    function obj = plot_boot_DP(obj)
        cmpXunq = unique(obj.cmpX)';
        mm=[min(cmpXunq) max(cmpXunq)];

        c=polyfit(cmpXunq',obj.DPfit,1);
        f=@(x) c(1)*x + c(2);
        x=linspace(mm(1),mm(2),100);
        y=f(x);
        plot(x,y,obj.lcolor,'LineWidth',obj.LineWidth); hold on

        plot(cmpXunq,obj.DPfit,[obj.shape obj.lcolor],'markersize',obj.markersize,'LineWidth',obj.LineWidth,'markerface',obj.markerface);

        titl=[ 'T='        num2str(obj.tFit,'%2.2f') ', N=' num2str(numel(obj.RcmpChs)) ];
        formatFigure(obj.Xname,'d''',titl);
        axis square
        hold off
    end
    function [] = plot_boot_T(obj)
        plot(obj.stdX,obj.tMU,[obj.shape obj.lcolor],'LineWidth',obj.LineWidth); hold on
        hold off
    end
    function [] = errorbar_boot_T(obj)
        %? tCI need to be subtracted?
        errorbar(obj.stdX,obj.tMU,obj.tCI,[obj.shape obj.lcolor],'LineWidth',obj.LineWidth); hold on
        hold off
    end
    function [] = plot_boot_mu_p(obj)
        obj.format_param_plot('\mu',obj.mMU,obj.mCI,obj.mFitDstb);
    end
    function [] = plot_boot_sigma_p(obj)
        obj.format_param_plot('\sigma',obj.sMU,obj.sCI,obj.sFitDstb);
    end
    function [] = plot_boot_bet_p(obj)
        obj.format_param_plot('\beta',obj.bMU,obj.bCI,obj.bFitDstb);
    end
    function [] = plot_boot_T_p(obj)
        obj.format_param_plot('T',obj.tMU,obj.tCI,obj.tFitDstb);
    end
    function [] = format_param_plot(obj,name,mu,CI,dst)
        [H,B] = hist(dst,21);
        plot(B,H,'color',obj.lcolor,'LineWidth',obj.LineWidth);
        formatFigure(name,'Num Samples',['m=' num2str(mu,  '%2.2f') ', ' num2str(obj.CIsz) '%=[' num2str(CI(1), '%2.2f') ',' num2str(CI(2),'%2.2f')  ']']);
        lim=[min(B) max(B)];
        l=max(abs(mu-lim));
        xlim([mu-l mu+l])
        %Text(.1,.9,{[num2str(obj.prcntUse) '% Data Used']},'ratio',18);
        axis square
        hold off
    end

    %% IND PLOT
    function [] = plot_ind_fit(obj,PCdta)
        % PLOT FIT (IN HI-RES)
        XcmpPlt = linspace(min(obj.cmpXsel),max(obj.cmpXsel),obj.nX);
        [PCplt,T]=obj.gen_gauss(XcmpPlt,obj.mFitSel,obj.sFitSel,obj.bFitSel); hold on;

        cmpXUnq = unique(obj.cmpX)';

        plot(XcmpPlt,PCplt,'color',obj.lcolor,'linewidth',1.5); hold on

        % RAW DATA- COMPUTE PERCENT COMPARISON CHOSEN
        [PCdta] = obj.get_PC();

        plot(cmpXUnq,PCdta,obj.shape,'color',obj.lcolor,'Linewidth',obj.LineWidth,'markersize',obj.markersize,'markerface',obj.markerface);

        % WRITE STUFF TO SCREEN
        writeText(1-.1,.1,{['n=' num2str(numel(obj.RcmpChsSel))]},'ratio',18,'right')
        formatFigure('','',['T=' num2str(T,'%.2f') ': \mu=' num2str(obj.mFitSel,'%1.2f') ',\sigma=' num2str(obj.sFitSel,'%1.2f') ',\beta=' num2str(obj.bFitSel,'%1.2f')]);
        xlim([min(obj.cmpXsel) max(obj.cmpXsel)]+[-.1 .1]); ylim([0 1])
        axis square
        hold off
    end
    function [] = plot_gen_gauss(obj,T)
        % XXX
        plot(cmpXsel,PC,'color',color,'linewidth',2); hold on
        formatFigure('','',['T=' num2str(T,'%.2f') ': \mu=' num2str(obj.mFitSel,'%1.2f') ',\sigma=' num2str(obj.sFitSel,'%1.2f') ',\beta=' num2str(obj.bFitSel,'%1.2f') ',nIntrvl=' num2str(obj.nIntrvl)]);
        xlim(minmax(obj.cmpXsel)+[-.1 .1]);
        ylim([0 1])
        axis square

        % WRITE PARAMETER VALUES AND N SMP TO SCREEN
        writeText(.075,.900,{['d''= ( | x - \mu| / \sigma )^{\beta}' ]},'ratio',20)
        writeText(.075,.775,{['d''_{crit}=' num2str(obj.DPcrt,'%.2f')]},'ratio',20)

        color='k'

        % THRESHOLD LINES
        plot(obj.mFitSel*[1 1],        [0                  0.5],'color',obj.lcolor,'linewidth',1);
        plot((obj.mFitSel+T)*[1 1],    [0     normcdf(0.5.*sqrt(obj.nIntrvl).*obj.DPcrt)],'color',obj.lcolor,'linewidth',1);
        plot([min(xlim) obj.mFitSel],  [0.5                0.5],'color',obj.lcolor,'linewidth',1);
        plot([min(xlim) obj.mFitSel+T],[normcdf(0.5.*sqrt(obj.nIntrvl).*obj.DPcrt)*[1 1]],'color',obj.lcolor,'linewidth',1);
        hold off
    end
%% Misc
    function obj=get_bP(obj)
        obj.bP=[isempty(obj.mFix) isempty(obj.mFix)  isempty(obj.bFix)];
    end

    function obj=get_testX(obj)
            if isempty(obj.cmpX) && isemtpy(obj.RcmpChs) && (isempty(obj.stdX) || numel(obj.stdX)==1)
                nTrlPerLvl = 100;
                lvls       = -2:.5:2;
            elseif ~isempty(obj.cmpX)
                lvls=unique(obj.cmpX(:));
                nTrlsPerLvl=numel(lvls)/size(obj.cmpX,1);
            elseif ~isemtpy(obj.RcmpChs) && isempty(obj.stdX)
                nTrls=size(obj.RcmpChs(:));
            elseif numel(obj.stdX > 1)
                nTrls=size(obj.stdX(:));
                ctr=obj.stdX(1);
            elseif numel(obj.stdX)==1 && ~ismpty(obj.RcmpChs)
                ctr=obj.stdX(1);
                nTrls=size(obj.RcmpChs(:));
            elseif numel(obj.stdX)==1
                ctr=obj.stdX(1);
            end

            if exist('ctr','var') && exist('nTrls','var')
                f=divisor(nTrls);
                d=abs(f-9);
                ind=find(d==min(d),1,'last')
                nTrlsPerLvl=f(ind);

                nLvls=nTrls/nTrlsPerLevel;

                r=2; % XXX
                R=ctr+r;
                L=ctr-r;
                lvls=linspace(R,L,nLvls)
            elseif exist('ctr','var')
                nTrlPerLvl = 100;
                nLvls=9;

                r=2; % XXX
                R=ctr+r;
                L=ctr-r;
                lvls=linspace(R,L,nLvls)
            elseif exist('nTrls','var')
                f=divisor(nTrls);
                d=abs(f-9);
                ind=find(d==min(d),1,'last')
                nTrlsPerLvl=f(ind);

                nLvls=nTrls/nTrlsPerLevel;
            end

            if isempty(obj.cmpX)
                obj.cmpX = repmat(lvls,nTrlPerLvl,1);
                obj.cmpX = obj.cmpX(:);
            end
            if isempty(obj.stdX)
                obj.stdX = 0;
            end

            if isempty(obj.stdX)
                obj.RcmpChs = binornd(1,repmat(normcdf(unique(obj.cmpX)'),nTrlPerLvl,1));
            end
    end
end
end
