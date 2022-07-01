classdef psyCurve < handle
properties
    %colors=[[18,133,158];[49,97,107];[45,209,148];[215,103,100];[158,18,73];[18,133,158]];
    %lcolors=colors ./ 255;
    mFix
    sFix
    bMFix
    bFix
    DPcrt
    nIntrvl
    bBoot
    nBoot
    bBest
    nBest
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
    negLL

    mFitDstb
    sFitDstb
    bFitDstb
    tFitDstb
    negLLDstb

    Ymin
    Ymax
    Xmin
    Xmax
end
properties(Hidden=true)
    titl
    bP
    bTest
    i

    stdX
    stds
    cmpX
    RcmpChs

end
methods
    function obj=psyCurve(stdX,cmpX,RcmpChs,varargin)
        bTest=0;
        if nargin < 1 || isempty(cmpX)
            bTest=1;
        else
            obj.stdX=stdX;
        end
        if nargin < 2 || isempty(cmpX)
            bTest=1;
        else
            obj.cmpX=cmpX;
        end
        if nargin < 3 || isempty(RcmpChs)
            bTest=1;
        else
            obj.RcmpChs=RcmpChs;
        end
        if nargin < 4 || isempty(varargin)
            Opts=struct;
        elseif length(varargin)==1 && isstruct(varargin{1})
            Opts=varargin{1};
        else
            Opts=struct(varargin{:});
        end
        obj.parse_opts(Opts);

        if obj.bTest
            obj.gen_testX();
        end
        % XXX make sure std cmp Rchs are same size

        [n,m]=size(obj.stdX);
        if n < m && n > 1 && m > 1
            obj.stdX=obj.stdX';
        end
        obj.stds=obj.stdX(1,:);
        obj.stdX=psyCurve.rmUniformCols(obj.stdX);
        if ~Set.isUniform(obj.stdX,1)
            error('stdX contains multiple values');
        end
        obj.cmpX=psyCurve.rmUniformCols(obj.cmpX);
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
        if obj.flipCmpChs
            obj.RcmpChs=~obj.RcmpChs;
        end

        obj.nTrials=length(obj.stdX);
        obj.nLvls=numel(unique(obj.cmpX));
        obj.nTrialsPerLvl=obj.nTrials/obj.nLvls;

        % SET FMINCON OPTIONS
        %obj.minFuncType = 'fmincon'; % XXX
        if strcmp(obj.minFuncType,'fmincon')
            obj.fitOpts             = optimset('fmincon');
            obj.fitOpts.MaxPCGIter=[];
            obj.fitOpts.Algorithm   = 'active-set';
            obj.fitOpts.LargeScale  = 'off';
            obj.fitOpts.UseParallel = 'never';
            obj.fitOpts.Display     = 'none';
            obj.fitOpts.MaxIter     = 500;
        elseif strcmp(obj.minFuncType,'fminsearch')
            obj.fitOpts             = optimset('fminsearch');
            obj.fitOpts.UseParallel = 'never';
            obj.fitOpts.Display     = 'off';
            obj.fitOpts.MaxIter     = 500;
        end

        obj.run();
    end
    function obj=parse_opts(obj,Opts)
        if ~exist('Opts','var')
            Opts=struct;
        end
        P={...
            'mFix',[],'Num.is';
            'sFix',[],'Num.is';
            'bFix',1,'Num.is';
            'bMFix',[],'Num.is';
            'DPcrt',1,'Num.is'; % XXX 1.36?
            'nIntrvl',2,'Num.isInt';
            'nBoot',[],'Num.isInt';
            'bBoot',[],'Num.is';
            'CIsz',68,'Num.is';
            'prcntUse',100,'Num.is';
            'minFuncType','fmincon','ischar';
            'flipCmpChs',0,'isbinary';
            ...
            'bBest',1,'isbinary';
            'nBest',[],'Num.is';
            'bPlot',0,'isbinary';
            'bPlotCI',[],'isbinary';
            'LineWidth',2,'Num.is';
            'markersize',12,'Num.is';
            'markerface','w','@(x) true';
            'color','k','@(x) true';
            'lcolor','k','@(x) true';
            'CIalpha',.4,'Num.is';
            'shape','o','@(x) true';
            'Xname','X','ischar';
            'Xunits','','ischar_e';
            'nX',200,'Num.is';
            'bExtendCurve',0,'isbinary';
            'bTitle',1,'isbinary';
            ...
            'bTest',0,'isbinary';
        };

        Args.parse(obj,P,Opts);
        if ~isempty(obj.bMFix) && obj.bMFix && isempty(obj.mFix)
            obj.mFix=unique(obj.stdX(:,1));
        end
        if isempty(obj.bBoot)
            if isempty(obj.nBoot) && isempty(obj.bPlotCI)
                obj.bBoot=false;
            elseif isempty(obj.nBoot)
                obj.bBoot=obj.bPlotCI;
            else
                obj.bBoot = obj.nBoot > 1;
            end
        end
        if isempty(obj.nBoot) && obj.bBoot
            obj.nBoot=1000;
        end

        if isempty(obj.bBest)
            if isempty(obj.nBest)
                obj.bBest=false;
            else
                obj.bBest = obj.nBest > 1;
            end
        end
        if isempty(obj.nBoot) && obj.bBest
            obj.nBest=50;
        end

        if isempty(obj.bPlotCI)
            obj.bPlotCI=obj.bBoot;
        end


    end
    function obj=run(obj,bRefit)
        if nargin < 2
            bRefit=false;
        end
        if bRefit
            lMFit=obj.mFit;
            lTFit=obj.tFit;
            lSFit=obj.sFit;
            lBFit=obj.bFit;
            lNegLL=obj.negLL;
        end
        if obj.bBest
            obj.fit_best();
        else
            obj.fit_basic();
        end
        if bRefit && obj.negLL < lNegLL
            obj.mFit  =lMFit;
            obj.tFit  =lTFit;
            obj.sFit  =lSFit;
            obj.bFit  =lBFit;
            obj.negLL =lNegLL;
        end

        if obj.bBoot
            obj.fit_boot();
        end
        if obj.bPlot==1
            obj.plot();
        end
    end
%- MAIN
    function obj=fit_basic(obj)
    %NON BOOTSTRAPPED FIT
        [RcmpChs,cmpX,stdX]=obj.get_data();
        S = obj.init_sel(1);

        S = obj.fit(S,cmpX,RcmpChs,1);

        obj.pack_dist(S);
        obj.select_from_dist(1);
    end
    function obj=fit_best(obj)
        [RcmpChs,cmpX,stdX]=obj.get_data();
        S=obj.init_sel(obj.nBest);

        for i = 1:obj.nBest
            S = obj.fit(S,cmpX,RcmpChs,i);
            S=  obj.gen_gauss_sel(S,i,cmpX);
        end
        obj.pack_dist(S);
        obj.select_best();
    end
    function obj=fit_boot(obj)
        [RcmpChs,cmpX,stdX]=obj.get_data();
        S=obj.init_sel(obj.nBoot);

        nSmp=round(numel(RcmpChs).*obj.prcntUse./100);
        for i = 1:obj.nBoot
            obj.i=i;

            % Sample
            [stdXsel,cmpXsel,RcmpChsSel]=obj.get_sample_data(nSmp);
            S = obj.fit(S,cmpXsel,RcmpChsSel,i);
            S=  obj.gen_gauss_sel(S,i,cmpXsel);
        end
        obj.pack_dist(S);
        obj.pack_boot();
        if ~obj.bBest
            obj.select_mean();
        end
    end
%% INIT
    function S=init_sel(obj,n)
        S=struct('mFit',[],'sFit',[],'bFit',[],'negLL',[],'PCfit',[],'tFit',[],'DPfit',[]);
        S=repmat(S,n,1);
    end
%% GET DATA
    function [RcmpChs,cmpX,stdX]=get_data(obj)
        nnanind=~isnan(obj.stdX) & ~isnan(obj.cmpX) & ~isnan(obj.RcmpChs);
        RcmpChs=logical(obj.RcmpChs(nnanind));
        RcmpChs=RcmpChs(:);
        cmpX=obj.cmpX(:);
        stdX=obj.stdX(:);
    end
    function [stdXsel,cmpXsel,RcmpChsSel]=get_sample_data(obj,nSmp)
        indSmp = randi(nSmp,nSmp,1);
        stdXsel=obj.stdX(indSmp);
        cmpXsel=obj.cmpX(indSmp);
        RcmpChsSel=obj.RcmpChs(indSmp);
    end
%% PACK
    function obj=pack_dist(obj,S)
        obj.mFitDstb=vertcat(S.mFit);
        obj.sFitDstb=vertcat(S.sFit);
        obj.bFitDstb=vertcat(S.bFit);
        obj.tFitDstb=vertcat(S.tFit);
        obj.negLLDstb=vertcat(S.negLL);
        %obj.PCfitDstb=vertcat(S.PCfit); % XXX
        %obj.DPfitDstb=vertcat(S.DPfit); % XXX

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
%% SELECT OUTPUT
    function select_err(obj)
        [Xdat,Ydat]=obj.get_dataXY();
        Ydat=Vec.row(Ydat);
        ctr=ceil(length(Xdat/2));
        %Ydat(ctr)=[];
        %Xdat(ctr)=[];

        E2=zeros(obj.nBoot,1);
        for i = 1:length(obj.nBoot)
            Y = obj.gen_gauss(Xdat,obj.mFitDstb(i,:),obj.sFitDstb(i,:),obj.bFitDstb(i,:));
            E2(i)=mean((Y-Ydat).^2);
        end
        ind=find(E2==min(E2),1,'first');

        obj.select_from_dist(ind);
    end
    function select_best(obj)
        ind=find(obj.negLLDstb==min(obj.negLLDstb),1,'first');

        obj.select_from_dist(ind);
    end
    function select_mean(obj)
        obj.mFit = obj.mMU;
        obj.sFit = obj.sMU;
        obj.bFit = obj.bMU;

        [RcmpChs,cmpX,stdX]=obj.get_data();
        p0 = [obj.mFit obj.sFit obj.bFit];
        obj.negLL=negLLFunc(p,cmpX,RcmpChs,obj.DPcrt,obj.nIntrvl);

        obj.get_secondary_stats(cmpX);
        obj.tFit = obj.tMU;
    end
    function select_from_dist(obj,ind)
        obj.mFit   = obj.mFitDstb(ind,:);
        obj.sFit   = obj.sFitDstb(ind,:);
        obj.bFit   = obj.bFitDstb(ind,:);
        obj.tFit   = obj.tFitDstb(ind,:);
        obj.negLL  = obj.negLLDstb(ind,:);
        obj.get_secondary_stats();
    end
%% FIT
    function S = fit(obj,S,cmpX,RcmpChs,i)


        % SET LOWER AND UPPER BOUNDS ON PARAMETERS
        pLB     = [2.0.*min(cmpX-mean(cmpX))+mean(cmpX) 0.02.*(max(cmpX)-min(cmpX)) 0.35];
        pUB     = [2.0.*max(cmpX-mean(cmpX))+mean(cmpX) 2.00.*(max(cmpX)-min(cmpX)) 3.00];


        % SET INITIAL PARAMETER VALUES
        m0  = obj.mFix;
        s0  = obj.sFix;
        b0  = obj.bFix;
        if isempty(m0);
            m0 = mean([min(cmpX), max(cmpX)]);
            m0 = m0  + .1.*randn;
        end
        if isempty(s0);
            s0 = diff(Num.minMax(abs(cmpX)))./6;
            s0 = s0  + .1.*s0.*randn;
        end
        if isempty(b0);
            b0 = 1;
            b0 = b0  + .1.*b0.*randn;
        end
        %s0

        p0 = [m0 s0 b0];

        fun=@(p) psyCurve.negLLFunc(p,cmpX,RcmpChs,obj.DPcrt,obj.nIntrvl,obj.mFix,obj.sFix,obj.bFix);
        % MINIMIZE NEGATIVE LOG-LIKELIHOOD
        switch obj.minFuncType
        case 'fmincon'
            [pFit,S(i).negLL] = fmincon(fun,p0,[],[],[],[],pLB,pUB,[],obj.fitOpts);
        case 'fminsearch'
            [pFit,S(i).negLL] = fminsearch(neg_LL,p0,obj.fitOpts);
        otherwise
            error()
        end

        % FINAL FIT PARAMETERS
        S(i).mFit=pFit(1);
        S(i).sFit=pFit(2);
        S(i).bFit=pFit(3);
    end
    function plot_negLL(obj)
        [RcmpChs,cmpX,stdX]=obj.get_data();
        m=unique(stdX);
        s=.001:.001:.02;
        b=1;

        % SET INITIAL PARAMETER VALUES
        nl=zeros(length(s),1);
        m0  = obj.mFix;
        s0  = obj.sFix;
        b0  = obj.bFix;
        if isempty(m0); m0 = mean([min(cmpX) max(cmpX)]); m0 = m0  + .1.*randn; end
        if isempty(s0); s0 = diff(Num.minMax(abs(cmpX)))./6;        s0 = s0  + .1.*s0.*randn; end
        if isempty(b0); b0 = 1;                                 b0 = b0  + .1.*b0.*randn; end
        p0 = [m0 s0 b0];

        fun=@(p) psyCurve.negLLFunc(p,cmpX,RcmpChs,obj.DPcrt,obj.nIntrvl,obj.mFix,obj.sFix,obj.bFix);
        for i = 1:length(s)
            nl(i)=fun([p0(1),s(i),p0(3)]);
        end
        plot(s,nl);
    end
    function S=gen_gauss_sel(obj,S,i,cmpX)
        [S(i).PCfit,S(i).tFit,S(i).DPfit] = obj.gen_gauss(cmpX,S(i).mFit,S(i).sFit,S(i).bFit,obj.DPcrt,obj.nIntrvl);
    end
    function get_secondary_stats(obj,cmpXsel)
        if nargin < 2
            [~,cmpXsel,~]=obj.get_data();
        end
        [obj.PCfit,obj.tFit,obj.DPfit] = obj.gen_gauss(unique(cmpXsel),obj.mFit,obj.sFit,obj.bFit,obj.DPcrt,obj.nIntrvl);
    end
%%
    function [PC,N1,N,N0]=get_PC(obj,RcmpChsSel,cmpXsel,stdXsel)
        % STANDARD VALUES
        if nargin < 2
            [RcmpChsSel,cmpXsel,stdXsel]=obj.get_data();
        end
        XstdUnq = unique(stdXsel);

        % LOOP OVER STANDARDS
        for s = 1:length(XstdUnq)
            % INDICES FOR EACH STANDARD
            indS = stdXsel == XstdUnq(s);
            % COMPARISON VALUES
            XcmpUnq(:,s) = unique(cmpXsel(indS));
            % LOOP OVER COMPARISONS
            for c = 1:length(XcmpUnq(:,s))
                % INDICES IN STD / CMP CONDITION
                indCnd  = stdXsel ==XstdUnq(s) & cmpXsel==XcmpUnq(c,s);
                % TOTAL NUMBER OF TRIALS IN CONDITION
                N(c,s)  = sum( indCnd );
                % TOTAL NUMBER OF CMP CHOSEN IN CONDITION
                N1(c,s) = sum( RcmpChsSel(indCnd)==1 );
                % TOTAL NUMBER OF STD CHOSEN IN CONDITION
                N0(c,s) = sum( RcmpChsSel(indCnd)==0 );
            end
        end

        PC = N1./N;
    end
%% PLOT
    function Plot(obj,meas,units,mult,xfrmt);
        if nargin < 2
            meas=[];
        end
        if nargin < 3
            units=[];
        end
        if nargin < 4
            mult=[];
        end
        if nargin < 5
            xfrmt=[];
        end
        if nargin < 6
            stdNames=[];
        end
        obj.plot();

        obj.Format(meas,units,mult,xfrmt);
        hold off;
    end
    function Format(obj,meas,units,mult,xfrmt)
        if nargin < 2
            meas=1;
        end
        if nargin < 3
            units=1;
        end
        if nargin < 4
            mult=1;
        end
        if nargin < 5
            xfrmt=[];
        end
        if nargin < 6
            stdNames=[];
        end

        obj.xylim();
        obj.ylabel();
        obj.rylabel(mult);
        obj.cmplabel(meas,units);
        obj.cmpticks(mult,xfrmt);

        Axis.format;
    end
    function [] = plot_all(obj)
        Fig.new();
        hold off;
        subPlot([2,2],1,1);
        obj.plot_boot_curve();
        obj.plot_boot_params(2);
        subPlot([2 2],2,1);
        obj.plot_boot_DP();
    end
    function [] = plot(obj)
        obj.plot_boot_curve();
        hold on;
        obj.plot_data();
    end
    function [] = plot_params(obj)
        Fig.new();
        hold off;
        obj.plot_boot_params();
    end
    function [] = plot_DP(obj)
        Fig.new();
        hold off;
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

    function [Xdat,Ydat]=get_dataXY(obj)
        Xdat = transpose(unique(obj.cmpX));
        % PROPORTION CMP CHOSEN
        Ydat=zeros(length(Xdat),1);
        for i = 1:length(Xdat)
            ind = obj.cmpX(:) == Xdat(i);
            Ydat(i) = mean(obj.RcmpChs(ind));
        end
    end
    function plot_data(obj,sym,color)
        if ~exist('sym','var') || isempty(sym)
            sym=obj.shape;
        end
        if ~exist('color','var') || isempty(color)
            color=obj.lcolor;
        end
        [Xdat,Ydat]=obj.get_dataXY();
        % UNIQUE COMPARISON VALUES
        plot(Xdat,Ydat,sym,'color',color, 'markerface',obj.markerface,'markersize',obj.markersize,'LineWidth',obj.LineWidth);
    end
    function xylim(obj,Xdat)
        [Xdat,~]=obj.get_dataXY();
        d=max(Xdat) - min(Xdat);
        m=d*0.05;
        xlim([min(Xdat)-m max(Xdat)+m]);
        ylim([0 1]);
    end
    function ylabel(obj)
        ylabel('Prop. Cmp Chosen');
        yt=num2cell(yticks);
        yt=cellfun(@(x) sprintf('%.1f',x),yt,'UniformOutput',false);
        yticklabels(yt);;
    end
    function rylabel(obj,mult)
        if nargin < 2 || isempty(mult)
            mult=1;
        end
        text=obj.statsText(mult,newline);

        yyaxis right;
        ylabel(text);
        yticks('');
        yticklabels('');
        set(gca,'ycolor','k');
        set(get(gca,'ylabel'),'rotation',0,'VerticalAlignment','middle','HorizontalAlignment','left') ;
        yyaxis left;

    end
    function strs=stdstr(obj,mult, names, bUnit)
        if nargin < 2 || isempty(mult)
            mult=1;
        end
        if nargin < 3
            names=[];
        end
        if nargin < 4
            bUnit=true;
        end
        n=length(obj.stds);
        strs=cell(n,1);
        vals=Vec.row(mult).*obj.stds;
        for i = 1:n
            if isempty(names)
                if n==1
                    I='X';
                else
                    I=['X' num2str(i)];
                end
            elseif length(names)==n
                I=names{i};
            elseif iscell(names)
                I=names{i};
            else
                I=names;
            end
            if bUnit
                I=[I ' '];
            else
                I='';
            end
            strs{i}=[I Num.toStr(vals(i))];
        end
    end
    function cmplabel(obj,meas,units)
        if nargin < 1 || isempty(meas)
            meas='x';
        end
        if nargin < 2 || isempty(units)
            units='a.u.';
        end
        meas(1)=upper(meas(1));
        units=['(' units ')'];
        xlabel([ meas ' ' units ]);
    end
    function cmpticks(obj,mult,xfrmt,Xdat)
        if nargin < 2 || isempty(mult)
            mult=1;
        end
        if nargin < 3
            xfrmt=[];
        end
        if nargin < 4 || isempty(Xdat)
            [Xdat,~]=obj.get_dataXY();
        end
        xticks(Xdat);
        if mult==1 && isempty(xfrmt)
            return
        end
        lbls=xticklabels;
        if isempty(xfrmt)
            spl=strsplit(lbls{1},'.');
            n=numel(spl{2});
            xfrmt=['%.' num2str(n) 'f'];
        end
        lbls=cellfun(@(x) num2str(str2double(x)*mult,xfrmt),lbls,'UniformOutput',false);
        xticklabels(lbls);
    end
    function obj = plot_boot_curve(obj)


        a=min(obj.cmpX);
        b=max(obj.cmpX);
        if obj.bExtendCurve
            rang=b-a;
            a=a-rang*obj.bExtendCurve;
            b=b+rang*obj.bExtendCurve;
        end
        if isempty(obj.nX)
            obj.nX=200;
        end

        Xfit = linspace(a,b,obj.nX);

        % PSYCHOMETRIC FUNCTION FIT MEAN
        Yfit = obj.gen_gauss(Xfit,obj.mFit,obj.sFit,obj.bFit,obj.DPcrt,obj.nIntrvl);
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

            c=Set.distribute(mlh,slh,blh);
            c=unique(c,'rows');
            Y=zeros(size(c,2),length(Xfit));
            for i = 1:size(c,1)
                ind=c(i,:);
                % m s b
                Y(i,:)  = obj.gen_gauss(Xfit,ind(1),ind(2),ind(3),obj.DPcrt,obj.nIntrvl);
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
            patch([Xfit fliplr(Xfit)], [Ymin fliplr(Ymax)], obj.color,'FaceAlpha',obj.CIalpha,'EdgeColor','none'); hold on;
        end

        hold on
        plot(Xfit,Yfit,'color',obj.lcolor,'LineWidth',obj.LineWidth); hold on;


        % Data

    end
    function obj=format_boot_curve(obj,titl,units)
        if nargin < 2 || isempty(titl)
            titl='';
        else
            titl=[titl newline];
        end
        if nargin < 4 || isempty(units)
            units=obj.Xname;
        end
        if strcmp(units,'none')
            units='';
        end

        if obj.bTitle
        else
            titl=[];
        end
        %Axis.format(units,'Proportion Cmp Chosen',titl);
        Axis.format(units,'Prop. Cmp Chosen',titl);
        axis square;
    end

    function out=statsText(obj,mult,sep)
        if nargin < 2 || isempty(mult)
            mult=1;
        end
        if nargin < 3 || isempty(sep)
            %sep=newline;
            sep=', ';
        end
        out=[ ...
            'n=' num2str(size(obj.stdX,1)) sep ...
            '\mu='    num2str(mult*obj.mFit,'%2.2f') sep...
            '\sigma=' num2str(mult*obj.sFit,'%2.2f') sep...
            'T=' num2str(mult*obj.tFit,'%2.2f') sep...
            '\beta='  num2str(obj.bFit,'%2.2f') sep ...
        ];
    end
    function obj = plot_boot_DP(obj)
        cmpXunq = transpose(unique(obj.cmpX));
        mm=[min(cmpXunq) max(cmpXunq)];

        c=polyfit(transpose(cmpXunq),obj.DPfit,1);
        f=@(x) c(1)*x + c(2);
        x=linspace(mm(1),mm(2),100);
        y=f(x);
        plot(x,y,obj.lcolor,'LineWidth',obj.LineWidth); hold on;

        plot(cmpXunq,obj.DPfit,[obj.shape obj.lcolor],'markersize',obj.markersize,'LineWidth',obj.LineWidth,'markerface',obj.markerface);

        titl=[ 'T='        num2str(obj.tFit,'%2.2f') ', N=' num2str(numel(obj.RcmpChs)) ];
        Axis.format(obj.Xname,'d''',titl);
        axis square;
        hold off;
    end
    function [] = plot_boot_T(obj)
        plot(obj.stdX,obj.tMU,[obj.shape obj.lcolor],'LineWidth',obj.LineWidth); hold on;
        hold off;
    end
    function [] = errorbar_boot_T(obj)
        %? tCI need to be subtracted?
        errorbar(obj.stdX,obj.tMU,obj.tCI,[obj.shape obj.lcolor],'LineWidth',obj.LineWidth); hold on;
        hold off;
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
        Axis.format(name,'Num Samples',['m=' num2str(mu,  '%2.2f') ', ' num2str(obj.CIsz) '%=[' num2str(CI(1), '%2.2f') ',' num2str(CI(2),'%2.2f')  ']']);
        lim=[min(B) max(B)];
        l=max(abs(mu-lim));
        xlim([mu-l mu+l]);
        %Text(.1,.9,{[num2str(obj.prcntUse) '% Data Used']},'ratio',18);
        axis square;
        hold off;
    end

    %% IND PLOT
    function [] = plot_ind_fit(obj,PCdta)
        % PLOT FIT (IN HI-RES)
        [RcmpChsSel,cmpXsel,stdXsel]=obj.get_data();
        XcmpPlt = linspace(min(cmpXsel),max(cmpXsel),obj.nX);
        [PCplt,T]=obj.gen_gauss(XcmpPlt,obj.mFit,obj.sFit,obj.bFit,obj.nIntrvl); hold on;

        cmpXUnq = unique(obj.cmpX)';

        plot(XcmpPlt,PCplt,'color',obj.lcolor,'linewidth',1.5); hold on

        % RAW DATA- COMPUTE PERCENT COMPARISON CHOSEN
        [PCdta] = obj.get_PC();

        plot(cmpXUnq,PCdta,obj.shape,'color',obj.lcolor,'Linewidth',obj.LineWidth,'markersize',obj.markersize,'markerface',obj.markerface);

        % WRITE STUFF TO SCREEN
        Axis.writeText(1-.1,.1,{['n=' num2str(numel(RcmpChsSel))]},'ratio',18,'right');
        Axis.format('','',['T=' num2str(T,'%.2f') ': \mu=' num2str(obj.mFit,'%1.2f') ',\sigma=' num2str(obj.sFit,'%1.2f') ',\beta=' num2str(obj.bFit,'%1.2f')]);
        xlim([min(obj.cmpXsel) max(cmpXsel)]+[-.1 .1]); ylim([0 1]);
        axis square;
        hold off;
    end
    function [] = plot_gen_gauss(obj,T)
        % XXX
        [RcmpChsSel,cmpXsel,stdXsel]=obj.get_data();
        plot(cmpXSel,PC,'color',color,'linewidth',2); hold on;
        Axis.format('','',['T=' num2str(T,'%.2f') ': \mu=' num2str(obj.mFit,'%1.2f') ',\sigma=' num2str(obj.sFit,'%1.2f') ',\beta=' num2str(obj.bFit,'%1.2f') ',nIntrvl=' num2str(obj.nIntrvl)]);
        xlim(minmax(cmpXsel)+[-.1 .1]);
        ylim([0 1]);
        axis square;

        % WRITE PARAMETER VALUES AND N SMP TO SCREEN
        Axis.writeText(.075,.900,{['d''= ( | x - \mu| / \sigma )^{\beta}' ]},'ratio',20);
        Axis.writeText(.075,.775,{['d''_{crit}=' num2str(obj.DPcrt,'%.2f')]},'ratio',20);

        color='k';

        % THRESHOLD LINES
        plot(obj.mFit*[1 1],        [0                  0.5],'color',obj.lcolor,'linewidth',1);
        plot((obj.mFit+T)*[1 1],    [0     normcdf(0.5.*sqrt(obj.nIntrvl).*obj.DPcrt)],'color',obj.lcolor,'linewidth',1);
        plot([min(xlim) obj.mFit],  [0.5                0.5],'color',obj.lcolor,'linewidth',1);
        plot([min(xlim) obj.mFit+T],[normcdf(0.5.*sqrt(obj.nIntrvl).*obj.DPcrt)*[1 1]],'color',obj.lcolor,'linewidth',1);
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
                ind=find(d==min(d),1,'last');
                nTrlsPerLvl=f(ind);

                nLvls=nTrls/nTrlsPerLevel;

                r=2; % XXX
                R=ctr+r;
                L=ctr-r;
                lvls=linspace(R,L,nLvls);
            elseif exist('ctr','var')
                nTrlPerLvl = 100;
                nLvls=9;

                r=2; % XXX
                R=ctr+r;
                L=ctr-r;
                lvls=linspace(R,L,nLvls);
            elseif exist('nTrls','var')
                f=divisor(nTrls);
                d=abs(f-9);
                ind=find(d==min(d),1,'last');
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
methods(Static)
    function [PC,T,DP] = gen_gauss(cmpX,mFit,sFit,bFit,DPcrt,nIntrvl)

        DP = sign(cmpX-mFit).*(abs(cmpX-mFit)./sFit).^bFit; % abs() prevent complex numbers w. some betas
        PC = normcdf( 0.5.*sqrt(nIntrvl).*DP,0,1); % sign() reinstates sign of Xcmp-mFit
        T  = sFit.*DPcrt.^(1./bFit);
    end
    function out=negLLFunc(p,cmpX,RcmpChs,DPcrt,nIntrvl,mFix,sFix,bFix)
        if  nargin >= 6    && ~isempty(mFix)    p(1) = mFix; end
        if  nargin >= 7    && ~isempty(sFix)    p(2) = sFix; end
        if  nargin >= 8    && ~isempty(bFix)    p(3) = bFix; end

        out=-(sum(log(     psyCurve.gen_gauss(cmpX(RcmpChs==1),p(1),p(2),p(3), DPcrt, nIntrvl) )) + ...
            sum(log( 1 - psyCurve.gen_gauss(cmpX(RcmpChs==0),p(1),p(2),p(3), DPcrt, nIntrvl) )) );
    end
    function X=rmUniformCols(X)
        [n,m]=size(X);
        if n==1
            return
        end
        ind=false(1,m);
        for i = 1:m
            ind(i)=all(X(:,i)==X(1,i));
        end
        if all(ind)
            ind(1)=false;
        end
        X(:,ind)=[];
    end
end
end

