classdef EPsyCurve < handle
properties
    name
    lvlInds

    U

    partition
    S
    B

    % ALL
    SP
    DATA
    CURVE
    DVCORR
    C
    T
    MU
    X

    RHO
    MU1
    MU2
    CR1
    CR2
    CI

    STDX
    CMPX

    TCI
    bBoot
    ulvls
    lvls
    sz

    % cur
    num
    f
    i
    Sp
    Curve
    DVCorr
    c
    ulvl
    ulvlInds
    lvl

    nCmp

    bLog

    alias
    subj
    moude
    passes
    blks

end
% F
% R
% C
% X
methods(Static)
    function fname=get_fname(alias,subj,moude,passes,name)
        if nargin < 5 || isempty(name)
            name='';
        else
            name=[name '_'];
        end
        if passes==1
            pss='';
        else
            pss=['_' strrep(Num.toStr(passes),',','-')];
        end

        fname=['media/EP_' alias '_' name subj '_' num2str(moude) pss '.mat'];
    end
    function EP=load(alias,subj,moude,passes,name)
        fname=EPsyCurve.get_fname(alias,subj,moude,passes,name);
        S=load(fname,'EP');
        EP=S.EP;
    end
end
methods
    function obj=EPsyCurve(alias,subj,moude,passes,uNames,lvlInds,num,name)
        obj.bLog=true;
        obj.U=EUnits(uNames{:});
        if nargin < 2
            subj=':';
        end
        if nargin < 3
            moude=1;
        end
        if nargin < 4 || isempty(passes)
            passes=1;
        end
        if nargin < 5
            uNames={'X'};
        end
        if nargin < 6
            lvlInds=[];
        end
        if nargin < 7
            num=[];
        end
        if nargin < 8
            name='';
        end
        obj.passes=passes;
        obj.alias=alias;
        obj.num=num;
        obj.subj=subj;
        obj.moude=moude;
        obj.lvlInds=lvlInds;
        obj.name=name;
        obj.get_data();
    end
%- FITTING
    function Opts=parse(obj,varargin)
        Opts=struct();
        if length(varargin)==1 && isstruct(varargin{1})
            opts=varargin{1};
        else
            opts=struct(varargin{:});
        end
        if isfield(opts,'bBlk')
            Opts.bBlk=opts.bBlk;
            opts=rmfield(opts,'bBlk');
        else
            Opts.bBlk=false;
        end
        Opts.opts=opts;
        if Opts.bBlk
            obj.selDataPartition('B');
        else
            obj.selDataPartition('S');
        end
    end
    function getCurves(obj,varargin)
        opts=obj.parse(varargin{:});

        obj.T=zeros(obj.sz);
        obj.X=nan(obj.sz);
        obj.TCI=nan([2 obj.sz]);
        obj.bBoot=false(obj.sz);

        n=length(obj.lvls);
        p=Pr(n,1,'Fitting curves');
        for j = 1:size(obj.lvls,1)
            p.u();
            obj.get_curve(j,opts);
        end
        p.u();
    end
    function getDVCorr(obj,varargin)
        opts=obj.parse(varargin{:});

        obj.RHO=nan([obj.sz,obj.nCmp]);
        obj.MU1=zeros([obj.sz,obj.nCmp]);
        obj.MU2=zeros([obj.sz,obj.nCmp]);
        obj.CR1=zeros([obj.sz,obj.nCmp]);
        obj.CR2=zeros([obj.sz,obj.nCmp]);
        obj.CI=nan([2 obj.sz obj.nCmp]);

        n=length(obj.lvls);
        p=Pr(n,1,'Fitting Correlations');
        for j = 1:n
            %p.u();
            obj.get_DVcorr(j,opts);
        end
    end
%- PLOT
    function plotMagr(obj,varargin)
        if nargin < 2 || isempty(bSameRC)
            bSameRC=[false false];
        end
        if nargin < 3 || isempty(bSameDt)
            bSameDt=bSameRC;
        end
        bSameRC=[false false];
        bSameDt=bSameRC;
        bFlip=true;

        fn=obj.initFig('Magr');

        obj.init_subplot(fn,bSameRC,bSameDt,bFlip,varargin{:});

        for j = 1:length(obj.lvls)
            obj.sel(fn,j,2);
            obj.plot_magr();
        end
        obj.sel(fn);
        obj.label_curve_RC();
    end
    function plotCorr(obj,varargin)
        if nargin < 2 || isempty(bSameRC)
            bSameRC=[false false];
        end
        if nargin < 3 || isempty(bSameDt)
            bSameDt=bSameRC;
        end
        bSameRC=[false false];
        bSameDt=bSameRC;
        bFlip=true;

        fn=obj.initFig('Corr');

        % SUBPLOT
        obj.init_subplot(fn,bSameRC,bSameDt,bFlip,varargin{:});

        for j = 1:size(obj.lvls,1)
            [exitflag]=obj.sel(fn,j,2);
            %if exitflag || isempty(obj.DVCorr)
            %    continue
            %end
            obj.plot_corr();
        end
        obj.sel(fn);
        obj.label_curve_RC();
    end
    function histRho(obj,bRatio,varargin)
        if nargin < 2 || isempty(bRatio)
            bRatio=false;
        end
        RHO=obj.RHO(:);
        if bRatio;
            %RHO=sqrt((1-RHO)./RHO);
            RHO=sqrt((1-RHO)./RHO);
        end
        [name,units,mult]=obj.U{1,'name','units','mult'}{:};
        fn=obj.initFig(['RhoHist' num2str(double(bRatio))]);

        bSameRC=[true true];
        bSameDt=bSameRC;
        bFlip=true;
        obj.init_subplot(fn,bSameRC,bSameDt,bFlip,varargin{:});

        obj.sel(fn,1,2);
        Hist.plot(RHO,8,'bLog',bRatio);

        if bRatio
            xlabel('\sigma_i / \sigma_e');
            xlim([.1 10]);
            xticks([.1 .3 1 3 10]);
            a=gca;
            set(a,'XScale','log');
            XTickLabel = get(a,'XTick');
set(a,'XTickLabel',num2str(XTickLabel'))

        else
            xlabel('Between Pass Correlation');
            xticks([-1:.5:1]);
            xlim([-1 1]);
        end
        ylabel('Count');
        axis square;
        Axis.format();

    end
    function scatterRho(obj,bRatio,bMean)
        if nargin < 2 || isempty(bRatio)
            bRatio=false;
        end
        if nargin < 3 || isempty(bMean)
            bMean=false;
        end
        fn=obj.initFig(['RhoScatter' num2str(double(bRatio)) num2str(double(bMean))]);
        dim=2;
        o=size(obj.RHO,3);
        sz=[obj.sz o];

        bSameRC=[true true];
        bSameDt=bSameRC;
        bFlip=true;

        obj.init_subplot(fn,bSameRC,bSameRC,bFlip);

        colors=['k','r','y'];
        nseen=0;
        seen=[];
        RHO=obj.RHO;
        if bRatio;
            %RHO=sqrt((1-RHO)./RHO);
            RHO=sqrt((1-RHO)./RHO);
        end
        if bMean
            Y=zeros(1,length(obj.lvls));
            X=zeros(1,length(obj.lvls));
            IND=zeros(1,length(obj.lvls));
        end
        [name,units,mult]=obj.U{1,'name','units','mult'}{:};
        seen=[];
        for j = 1:length(obj.lvls)
            obj.sel(fn,j,2);
            if isempty(obj.DVCorr)
                continue
            end

            n=obj.ulvl(1);
            m=obj.ulvl(2);
            y=squeeze(RHO(n,m,:));
            if bMean
                y=mean(y);
            end
            x=squeeze(repmat(obj.STDX(j,1),1,numel(y)));
            if all(isnan(y)) || all(y==0)
                continue
            end
            lvls=obj.ulvls(j,:);
            if ismember(lvls(dim),seen)
                ind=find(ismember(seen,lvls(dim)));
            else
                nseen=nseen+1;
                seen=[seen; lvls(dim)];
                ind=nseen;
            end
            color=colors(ind);
            hold on;
            if bMean
                Y(j)=y;
                X(j)=x;
                IND(j)=ind;
            else
                plot(x*mult,y,'o','Color',color,'MarkerFaceColor',color);
            end
        end
        if bMean
            uinds=unique(IND);
            for i = 1:length(uinds)
                u=uinds(i);
                if u==0
                    continue
                end
                ind=IND==u;
                color=colors(u);
                plot(X(ind)*mult,Y(ind),'o-','Color',color,'MarkerFaceColor',color);
            end
        end
        obj.sel(fn);

        [name2,units,mult]=obj.U{2,'name','units','mult'}{:};
        [name, units,mult]=obj.U{1,'name','units','mult'}{:};
        if bRatio
            ylabel('\sigma_i / \sigma_e');
            ylim([.1 10]);
            yticks([.1 .3 1 3 10]);
            a=gca;
            set(a,'YScale','log');
            YTickLabel = get(a,'YTick');
set(a,'YTickLabel',num2str(YTickLabel'))

        else
            ylabel('\rho');
            ylim([0 1]);
        end
        xlabel([name ' (' units ')'] );

        n=obj.STDX(seen,dim);
        lgnd=cellfun(@(x) num2str(x), num2cell(n),'UniformOutput',false);
        legend(strcat(name2,{' '},lgnd),'Location','best');
        axis square;
        Axis.format();
        obj.label_curve_RC(false);
    end
    function plotRho(obj,varargin)

        if nargin < 2 || isempty(bSameRC)
            bSameRC=[true true];
        end
        if nargin < 3 || isempty(bSameDt)
            bSameDt=bSameRC;
        end
        bSameRC=[false false];
        bSameDt=bSameRC;
        bFlip=true;

        fn=obj.initFig('Rho');

        % SUBPLOT
        obj.init_subplot(fn,bSameRC,bSameDt,bFlip,varargin{:});

        dim=1;
        o=size(obj.RHO,3);
        sz=[obj.sz o];
        bGd=false(sz(dim),o);
        for j = 1:length(lvls)
            obj.sel(fn,j,2);
            %obj.plot_corr();
            hold on;
            bGd(j,:)=obj.plot_rho(dim,j);
        end
        obj.sel(fn);

        if dim==1
            thing=obj.U{2,'name'}{:};
        else
            thing=obj.U{1,'name'}{:};
        end

        %if dim==2
        %    bGd=permute(bGd,[2 1 3]);
        %end
        inds=find(bGd);
        [n,m]=find(bGd);
        [~,i]=sortrows([n,m]);
        inds=inds(i);
        stds=fliplr(sortrows(fliplr(obj.STDX(inds,:))));
        lgndT=cell(size(stds));
        for i = 1:size(stds,2)
            [frmt,mult,abbrev]=obj.U{i,'frmt','mult','nabbrev'}{:};
            lgndT(:,i)=arrayfun(@(x) sprintf(frmt,x) ,stds(:,i)*mult,'UniformOutput',false);
        end
        lgnd=cell(size(lgndT,1),1);
        for i = 1:size(lgnd,1)
            lgnd{i}=strjoin(lgndT(i,:),' ');
        end

        [name,units]=obj.U{1,'name','units'}{:};
        xlabel([name ' (' units ')']);
        ylabel('\rho');
        legend(lgnd{:});
        %obj.label_curve_RC();
    end
    function plot(obj,bSameRC,bSameDt,varargin)
        if nargin < 2 || isempty(bSameRC)
            bSameRC=[false false];
        end
        if nargin < 3 || isempty(bSameDt)
            bSameDt=bSameRC;
        end
        bFlip=true;

        fn=obj.initFig('Curves');

        % SUBPLOT
        obj.init_subplot(fn,bSameRC,bSameDt,bFlip,varargin{:});

        for j = 1:length(obj.lvls)
            obj.sel(fn,j);
            if isempty(obj.Curve)
                continue
            end
            hold on;
            obj.plot_curve();
        end
        obj.sel(fn);
        obj.label_curve_RC();
    end
    function plotT(obj,varargin)
        [bFlip,varargin,bSuccess]=Args.getPair('bFlip',varargin{:});

        fn=obj.initFig('Thresholds');

        bSameRC=false(1,2);
        bSameRC(1)=obj.partition=='S';
        if ~bSuccess
            bFlip=bSameRC(1);
        end
        bFlip=1;
        bSameDt=bSameRC;

        N=obj.init_subplot(fn,bSameRC,bSameDt,bFlip,varargin{:});
        n=numel(unique(N));
        for j = 1:n
            exitflag=obj.sel(fn,j);
            if exitflag || isempty(obj.Curve)
                continue
            end
            if obj.bBoot(j)
                obj.plot_thresh(j,obj.bBoot(j));
                hold on;
            end
            obj.plot_thresh(j,false);
            %obj.plot_thresh(j,bFlip);
            Axis.format;
            obj.label_curve();
        end
        obj.sel(fn);
        obj.label_curve_RC(false);
    end
    function save(obj)
        f=obj.f;
        Sp=obj.Sp;
        SP=obj.SP;

        obj.f=[];
        obj.Sp=[];
        obj.SP=[];
        fname=EPsyCurve.get_fname(obj.alias,obj.subj,obj.moude,obj.passes,obj.name);
        EP=obj;
        save(fname,'EP');

        obj.f=f;
        obj.Sp=Sp;
        obj.SP=SP;
    end
    function plotTT(obj,varargin)
        fn=obj.initFig('Thresholds Same');
        bFlip=true;
        bSameRC=[true true];
        bSameDt=[true false];
        obj.init_subplot(fn,bSameRC,bSameDt,bFlip,varargin{:});
        n=numel(obj.lvls);

        if ~isempty(obj.TCI)
            for j = 1:n
                exitflag=obj.sel(fn,j);
                if exitflag || isempty(obj.Curve)
                    continue
                end
                obj.plot_thresh(j,obj.bBoot(j));
            end
        end
        obj.resetSeen();
        for j = 1:n
            exitflag=obj.sel(fn,j);
            if exitflag || isempty(obj.Curve)
                continue
            end
            hold on;
            obj.plot_thresh(j,false);
            Axis.format;
            obj.label_curve();
            obj.lim(true);
        end
        obj.sel(fn);
        obj.ylabel();
        %lims=Axis.getLims('data');

        obj.label_curve_RC();
    end
    function fn=initFig(obj,parent)
        name=[obj.getName() ' ' parent];
        obj.f=findobj('Name',name);
        if ~isempty(obj.f)
            %close(obj.f);
            fn=get(obj.f,'Number');
            clf(obj.f,'reset');
            obj.f.Name=name;
            figure(obj.f)
        else
            obj.f=figure('Name',name);
            fn=get(obj.f,'Number');
        end
    end
    function refit(obj,n,m,val)
        obj.CURVE{n,m}.run;
        obj.plot;
    end
    function setM(obj,n,m,val)
        obj.CURVE{n,m}.mFit=val/60;
        obj.plot;
    end
    function setT(obj,n,m,val)
        obj.CURVE{n,m}.tFit=val/60;
        obj.CURVE{n,m}.sFit=val/60;
        obj.T(n,m)=val/60;
        obj.plot;
    end
    function name=getName(obj)
        name=[obj.alias ' ' obj.subj ' ' num2str(obj.moude) '_' num2str(obj.num)];
    end
    function [n,m]=getSubs(obj,j,~)
        %if nargin < 2 || isempty(bFlip)
            bFlip=false;
        %bFlipend
        lvls=obj.ulvls(j,:);
        if bFlip
            n=lvls(2);
            m=lvls(1);
        else
            n=lvls(1);
            m=lvls(2);
        end
        %if bFlip
        %    [m,n]=ind2sub(obj.sz,obj.ulvlInds(j));
        %else
        %    [n,m]=ind2sub(obj.sz,obj.ulvlInds(j));
        %end
    end
    function exitflag=sel(obj,fn,ii,typ)
        exitflag=false;
        if nargin < 3
            ii=[];
        end
        if nargin < 4 || isempty(typ)
            typ=1;
        end

        % FIG NUM
        if fn~=obj.f.Number
            if ~isempty(obj.c)
                obj.C{obj.f.Number}=obj.c;
            end
            obj.f=findobj('Number',fn);
            figure(obj.f)
            obj.Sp=obj.SP{fn};
        end
        if isempty(obj.c)
            obj.c=obj.C{fn};
        end

        obj.i=ii;
        if isempty(ii)
            exitflag=true;
            return
        end

        %[n,m]=ind2sub(obj.sz,obj.ulvlInds(ii));
        if strcmp(obj.partition,'B')
            o=obj.B.blks(ii);
        else
            o=1;
        end
        %if ~obj.bFlip
        %    nn=n;
        %    n=m;
        %    m=nn;
        %end
        % LVL NUM
        obj.lvl=obj.lvls(ii);
        obj.ulvl=obj.ulvls(ii,:);
        n=obj.ulvl(1);
        m=obj.ulvl(2);
        if typ==1
            obj.Curve=obj.CURVE{n,m,o};
        else
            obj.DVCorr=obj.DVCORR{n,m,o};
        end

        % SUBPLOT RC
        if all(obj.c.bSameRC)
            s=[1 1];
        elseif obj.c.bSameRC(1)
            s=[1 obj.ulvl(2)];
        elseif obj.c.bSameRC(2)
            s=[obj.ulvl(2) 1];
        else
            s=obj.ulvl;
        end
        if obj.c.bFlip
            s=fliplr(s);
        end
        obj.Sp.sel(s);

        % DATA RC
        if all(obj.c.bSameDt)
            lvl=1;
        elseif obj.c.bSameDt(1)
            lvl=obj.ulvl(2);
        elseif obj.c.bSameDt(2)
            lvl=obj.ulvl(1);
        else
            lvl=obj.ulvl;
        end
        if isempty(obj.c.bSeen)
            hold off;
        else
            hold on;
        end
        bSeen=ismember(lvl,obj.c.bSeen,'rows');
        if ~bSeen
            obj.c.bSeen(end+1,:)=lvl;
        else
            exitflag=1;
        end
    end
    function s=subsIndex(obj,j)
        [n,m]=obj.getSubs(j,false);
        s.type='()';
        if all(obj.c.bSameDt)
            s.subs={':'};
        elseif obj.c.bSameDt(1)
            s.subs={':',m};
        elseif obj.c.bSameDt(2)
            s.subs={n,':'};
        else
            s.subs={n,m};
        end
    end
    function [x,t,ci]=indexThresh(obj,j)
        [mult]=obj.U{1,'mult'};

        s=obj.subsIndex(j);
        x=subsref(obj.X,s)*mult;
        [x,ind]=sort(x);

        if obj.partition=='B';
            s.subs(end+1)={':'};
        end

        t=subsref(obj.T,s)*mult;
        t=t(ind,:,:);
        nn=isnan(t);
        t(nn)=[];

        if nargout < 2 || ~obj.bBoot(j)
            ci=[];
            return
        end
        s.subs=[':' s.subs];
        tCI=subsref(obj.TCI,s);
        tCI=tCI(:,ind,:)*mult;

        ci=cell(2,1);
        ci{1}=tCI(1,:,:);
        ci{2}=tCI(2,:,:);
        ci{1}(nn)=[];
        ci{2}(nn)=[];

    end
    function str=get_units(obj);
        units=obj.U{1,'units'}{:};
        str=['(' units ')'];
    end
    function colo=get_color(obj,j)
        if obj.partition=='B';
            g=1;
            colo=[0 0 0];
            return
        end

        G=linspace(1,0,obj.sz(1));
        g=G(j);
        colo=[1 g 0];
    end
    function lgnd=legend(obj,dim)
        lgnd=strcat( obj.U{dim,'meas'}{:} ,{' '},obj.c.clabels );
        legend(lgnd,'Location','northwest');
    end
    function ylabel(obj)
         ylabel(sprintf('Threshold (%s)',obj.U{1,'units'}{:}));
    end
    function lim(obj,bLogY)
        set(gca,'YScale','log');

        %d=maxi-mini;
        %m=d*0.20;
        %ylim([mini-m maxi+m]);

        h=gca;
        if bLogY
            %hold on
            set(h, 'YScale', 'log');
        end
        if obj.partition=='B'
            xlim([1 6]);
            ylim([.3 10]);
        else
            yt=[.3 1 3]; % XXX
            ylim(Num.minMax(yt));
            yticks(yt);
            if bLogY
                Axis.yticksLog(h);
            end
        end

        return
        d=max(x) - min(x);
        m=d*0.05;
        xlim([min(x)-m max(x)+m]);

        % XTICKS
        lbls=arrayfun(@(w) num2str(w,frmt),x,'UniformOutput',false);
        xticks(x);
        xticklabels(lbls);

        % YTICKS
        % ylim([10^-1 10]);
        % yticks([.1 1 10]);
        % yticklabels([.1 1 10]);
        yticklabels([.3 1 3]);

    end
    function xlabel(obj)
        [name,units,mult]=obj.U{1,'name','units','mult'}{:};
        str=sprintf('%s (%s)',name,obj.U{1,'units'}{:});
        xlabel(str);
    end
    function plot_thresh(obj,j,bBoot)
        colo=obj.get_color(j);
        [x,t,tCI]=indexThresh(obj,j);
        t(t==0)=nan;

        if bBoot
            Plot.interv(x,tCI{1},tCI{2},colo,'EdgeColor','none','FaceAlpha',0.3); hold on
        else
            if numel(x)==1 && obj.partition=='B'
                x=1:numel(t);
            end

            plot(x, squeeze(t),'-o','Color',colo,'MarkerEdgeColor','k','LineWidth',2,'MarkerSize',10,'MarkerFaceColor',colo);
        end

    end
    function plot_curve(obj)
        %obj.Curve.select_err();
        obj.Curve.Plot(obj.U{1,'meas','units','mult','frmt'}{:});
        obj.label_curve();
    end
    function plot_magr(obj)
        if ~isempty(obj.DVCorr)
            obj.DVCorr.plot_magr();
            %obj.label_rho();
        else
            Axis.format();
            axis square;
        end
        title('');
        xlabel(sprintf('Proportion Cmp Chosen'));
        ylabel(sprintf('Proportion Argreement'));
        obj.label_curve();
    end
    function bGd=plot_rho(obj,dim,j)
        colors=['r','k','b','m','k'];
        markers=['o','<','o','<'];
        ind=1;
        [name,mult,frmt,units]=obj.U{ind,'name','mult','frmt','units'}{:};
        %if nargin < 2 || isempty(rho)
            %[n,m]=ind2sub(obj.sz,obj.ulvlInds(obj.i));
        %end
        s=repmat({':'},1,3);
        s{dim}=j;
        rho=squeeze(subsref(obj.RHO,substruct('()',s)));

        o=size(obj.RHO,3);
        stds=obj.STDX;
        ustds=unique(stds(:,ind),'rows','stable');

        ulvls=obj.ulvls;
        if dim==1
            %ulvls=fliplr(sortrows(fliplr(ulvls)));
            %ulvls=fliplr(sortrows(fliplr(ulvls)));
        end
        mind=ulvls(j,1);
        %cind

        cmpx=squeeze(subsref(obj.CMPX*mult,substruct('()',s)));

        %x=repmat(,1,o);
        o=size(obj.RHO,3);
        s=repmat({':'},1,2);
        bGd=zeros(1,o);
        for i = 1:size(rho,1)
            s{1}=i;
            x=subsref(cmpx,substruct('()',s));
            y=subsref(rho,substruct('()',s));
            if all(y==0) || all(isnan(y))
                continue
            end
            bGd(i)=true;
            color=colors(i);
            plot(x,y,[markers(j) '-'],'Color',color,'MarkerFaceColor',color,'LineWidth',2);
            hold on;
        end
    end
    function label_rho(obj)
        yy=ylim;
        y=yy(1);
        xx=xlim;
        x=xx(end);
        [n,m]=ind2sub(obj.sz,obj.ulvlInds(obj.i));
        str=sprintf('\\rho=%.2f',obj.RHO(n,m,1));
        %ax=obj.Sp.text{obj.i};
        if isempty(ax)
            obj.Sp.text{obj.i}=text(x,y,str,'FontSize',14,'HorizontalAlignment','right','VerticalAlignment','bottom');
        else
            set(ax,'String',str);
        end
    end
    function plot_corr(obj)
        mult=obj.U{1,'mult'};
        if ~isempty(obj.DVCorr)
            obj.DVCorr.plot_ellipse_all_same();
            %obj.label_rho();
            %n=yticks;
            %yticklabels(n/mult);
            %xticklabels(n/mult);
        else
            Axis.format();
            n=[-4:2:4];
            l=[n(1) n(end)];
            ylim(l);
            xlim(l);
            yticks(n);
            xticks(n);
            yticklabels(n);
            xticklabels(n);
            %yticklabels(n)
            %xticklabels(n)
            axis square;
        end
        str='Decision Variable: Pass %d';
        %if obj.bFlip
        %    ps=flipud(obj.passes);
        %else
            ps=obj.passes;
        %end
        xlabel(sprintf(str,ps(1)));
        ylabel(sprintf(str,ps(2)));
        obj.label_curve();
    end
    function label_curve(obj)
        str=obj.std_str([],false);

        %if obj.ulvl(2)==1 || obj.c.bSameRC(1)
            obj.c.clabels{obj.ulvl(2)}=str{2};
        %end
        %if obj.ulvl(1)==1 || obj.c.bSameRC(2)
            obj.c.rlabels{obj.ulvl(1)}=str{1};
        %end
        obj.xlabel();
        obj.ylabel();
        obj.Sp.label('','',true);

        % TICKS
        %if obj.bFlip
        %    obj.Sp.xticks();
        %end
        obj.Sp.yticks();
    end
    function label_curve_RC(obj,bLegend)
        nam=obj.U{'nabbrev'};
        if nargin < 2 || isempty(bLegend)
            if obj.c.bSameDt(1)
                bLegend=true;
            else
                bLegend=false;
            end
        end

        % RClabel
        if ~obj.c.bSameRC(1)
            in=obj.c.rlabels;
            if ~obj.c.bFlip
                args={'right',  ['Std. ' nam{1}]};
                try
                    obj.Sp.rlabel(in,args{:});
                end
            else
                args={'top',  ['Std. ' nam{1}]};
                try
                    obj.Sp.clabel(in,args{:});
                end
            end
        end
        if ~obj.c.bSameRC(2)
            in=obj.c.clabels;
            if ~obj.c.bFlip
                args={'top',  ['Std. ' nam{2}]};
                try
                    obj.Sp.clabel(in,args{:});
                end
            else
                args={'right',  ['Std. ' nam{2}]};
                try
                    obj.Sp.rlabel(in,args{:});
                end
            end
        end

        % LEGEND
        % XXX correct?
        %
        if bLegend && obj.c.bSameDt(1)
            obj.legend(2);
        elseif bLegend && obj.cSameDt(2)
            obj.legend(1);
        end

        % SUPTITLE
        if obj.moude==1
            str='';
        elseif obj.moude==2
            str=' pilot';
        elseif obj.moude==3;
            str=' train';
        end
        if obj.passes==1
            ps='';
        else
            ps=[' ' strrep(Num.toStr(obj.passes),',',' ')];
        end
        titl=[obj.alias ' ' obj.subj str ps ];
        %titl=['{\fbox' titl '}];
        obj.Sp.suptitle(titl,22,'bottomleft');
    end
    function strs=std_str(obj,typ,bUnit)
        X=obj.STDX(obj.i,:);
        U=obj.U;

        if nargin < 2
            bUnit=false;
        end
        if isempty(typ)
            [mult,names]=U{'mult','name'};
        else
            [mult,names]=U{typ,'mult','name'};
        end

        vals=Vec.row(mult).*X;

        n=length(X);
        strs=cell(n,1);
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
    function N=init_subplot(obj,fn,bSameRC,bSameDt,bFlip,varargin)
        siz=obj.sz;
        if all(bSameRC)
            siz=[1,1];
            oMargin=[2 1 2 2];
            iMargin=[0 0 1 3];
            oUnits='char';
            iUnits='char';
        elseif bSameRC(1)
            siz=[1,siz(2)];
            oMargin=[2 1 2 2];
            iMargin=[0 0 1 3];
            oUnits='char';
            iUnits='char';
            N=obj.ulvls(:,2);
        elseif bSameRC(2)
            siz=[siz(1),1];
            oMargin=[2 1 2 2];
            iMargin=[0 0 1 3];
            oUnits='char';
            iUnits='char';
            N=obj.ulvls(:,2);
        else
            oMargin=[2 1 2 3];
            iMargin=[0 10 1 4];
            oUnits='char';
            iUnits='char';
            %N=1:siz(1);
            if obj.partition=='B'
                N=1:numel(obj.B.blks);
            else
                N=1:prod(siz);
            end
        end
        if bFlip
            Siz=fliplr(siz);
        else
            Siz=siz;
        end

        obj.c=[];
        obj.SP{fn}=[];
        obj.C{fn}=[];
        obj.C{fn}.bFlip=bFlip;
        obj.SP{fn}=SubPlot(Siz,'iMargin',iMargin,'oMargin',oMargin,'oUnits',oUnits,'iUnits',iUnits,varargin{:});
        obj.C{fn}.rlabels=cell(siz(1),1);
        obj.C{fn}.clabels=cell(siz(2),1);
        obj.C{fn}.bSameRC=bSameRC;
        obj.C{fn}.bSameDt=bSameDt;
        obj.C{fn}.bBoot=[];
        obj.Sp=obj.SP{fn};
        obj.c=obj.C{fn};
        obj.resetSeen();
    end
    function resetSeen(obj)
        if any(obj.c.bSameDt)
            n=1;
        else
            n=2;
        end
        obj.c.bSeen=zeros(0,n);
    end
%- DATA
    function get_data(obj)
        E=ETable.get(obj.alias);

        % TESTING
        %obj.moude=1;
        %obj.lvlInds=[1];

        sels={'subj', obj.subj, 'mode', obj.moude, 'status', 1, 'pass',obj.passes};
        if ~isempty(obj.lvlInds)
            sels=[sels, 'lvlInd', obj.lvlInds ];
        end
        if ~isempty(obj.blks)
            sels=[sels, 'blk', obj.blks ];
        end
        obj.DATA=E.loadRawData(sels{:});

        nc=numel(unique(obj.DATA{'cmpX'},'rows'));
        ns=numel(unique(obj.DATA{'stdX'},'rows'));
        obj.nCmp=nc/ns;


        gets={'lvlInd','stdX'};
        Stds=obj.DATA.unique_rows(gets{:});
        [obj.S]=stdFun(Stds,false);

        gets={'lvlInd','stdX','blk'};
        [StdsBlk,inds,ic]=obj.DATA.unique_rows(gets{:});
        dta=obj.DATA(inds);
        [obj.B]=stdFun(StdsBlk,true);
        %[dates,ind]=sort(,'ascend');

        obj.sz=max(obj.S.ulvls,[],1);
        function S=stdFun(Stds,bBlk)
            % NOTE NO SORTING GOING ON HERE

            S=struct();

            % LVLS
            if bBlk
                blk=Stds{'blk'};
                S.blks=blk;
            end
            lvls=Stds{'lvlInd'};
            stds=Stds{'stdX'};
            dates=Stds{'dates'};
            stdsAll=Set.distribute(unique(stds(:,1)),unique(stds(:,2)));

            % ulvlInds
            ulvlInds=zeros(size(stds,1),1);
            for j = 1:size(stds,1)
                ulvlInds(j)=find(ismember(stdsAll,stds(j,:),'rows'));
            end
            STDX=stdsAll(ulvlInds,:);

            % ULVLS
            ulvls=zeros(size(stds));
            %inds=zeros(size(stds));
            for j = 1:size(stds,2)
                [~,inds,ulvls(:,j)]=unique(stds(:,j));
            end
            %if obj.bFlip
            %    ulvls=fliplr(ulvls);
            %end

            S.lvls=lvls;
            S.ulvlInds=ulvlInds;
            S.STDX=STDX;
            S.ulvls=ulvls;


        end
    end
    function selDataPartition(obj,typ)
        switch typ
            case 'B'
                d=obj.B;
            case 'S'
                d=obj.S;
        end
        obj.partition=typ;
        obj.lvls=d.lvls;
        obj.ulvlInds=d.ulvlInds;
        obj.STDX=d.STDX;
        obj.ulvls=d.ulvls;
    end
%- FITTING
    function get_DVcorr(obj,i,Opts)
        obj.ulvl=obj.ulvls(i,:);
        obj.lvl=obj.lvls(i);

        n=obj.ulvl(1);
        m=obj.ulvl(2);

        [stdX,cmpX,bRCmp]=obj.get_pass_data();
        if isempty(stdX)
            return
        end


        %[stdX,cmpX,bRCmp]=obj.DATA{'stdX',obj.ulvl,'stdX','cmpX','bRCmp'};
        dvcorr=DVcorr(stdX, cmpX, bRCmp, Opts.opts);
        dvcorr.run(false);
        obj.bBoot(n,m)=dvcorr.nBoot > 1;
        if obj.bBoot(n,m)
            s=substruct('()',[':' num2cell([n m]) ]);
            obj.CI=subsasgn(obj.CI,s,dvCorr.CI);
        end
        obj.CMPX(n,m,:)=unique(cmpX);

        obj.RHO(n,m,:)=dvcorr.RHO;
        obj.MU1(n,m,:)=dvcorr.MU1;
        obj.MU2(n,m,:)=dvcorr.MU2;
        obj.CR1(n,m,:)=dvcorr.CR1;
        obj.CR2(n,m,:)=dvcorr.CR2;

        obj.DVCORR{n,m}=dvcorr;
    end
    function [stdX,cmpX,bRCmp]=get_pass_data(obj)
        N=length(obj.passes);
        stdX=[];
        cmpX=[];
        bRCmp=[];

        lvlInd=obj.lvl;
        for j = 1:N
            pass=obj.passes(j);

            [t,b,p,s,c,r]=obj.DATA{'pass',pass,'lvlInd',lvlInd,'trial','blk','pass','stdX','cmpX','bRCmp'};
            Rows=[t,b,c];
            if j>1 && size(stdX,1) ~= size(s,1)
                ind1=~ismember(lRows,Rows,'rows');
                ind2=~ismember(Rows,lRows,'rows');
                if any(ind1)
                    stdX(ind1,:)=[];
                    cmpX(ind1,:)=[];
                    bRCmp(ind1,:)=[];
                end
                if any(ind2)
                    s(ind2,:)=[];
                    c(ind2,:)=[];
                    r(ind2,:)=[];
                end
            end
            if isempty(s)
                %lvlInd
                return
            end
            stdX=[stdX s(:,1)];
            cmpX=[cmpX c(:,1)];
            bRCmp=[bRCmp r];
            lRows=Rows;
        end
    end
    function get_curve(obj,i,Opts)
        obj.ulvl=obj.ulvls(i,:);
        obj.lvl=obj.lvls(i);

        %[n,m]=ind2sub(obj.sz,obj.ulvlInds(i));
        n=obj.ulvl(1);
        m=obj.ulvl(2);
        %[n,m]=ind2sub(obj.sz,i);

        if Opts.bBlk
            o=obj.B.blks(i);
            sels={'lvlInd',obj.lvl,'blk',o};
        else
            sels={'lvlInd',obj.lvl};
            o=1;
        end
        gets={'stdX','cmpX','bRCmp'};
        in=[sels gets];
        [stdX,cmpX,bRCmp]=obj.DATA{in{:}};

        curve=psyCurve(stdX, cmpX, bRCmp, Opts.opts);

        obj.X(n,m)=stdX(1);
        obj.bBoot(n,m,o)=curve.bBoot;
        if obj.bBoot(n,m,o)
            s=substruct('()',[':' num2cell([n m]) ]);
            obj.TCI=subsasgn(obj.TCI,s,curve.tCI);
            obj.T(n,m,o)=curve.tMU;
        else
            obj.T(n,m,o)=curve.tFit;
        end
        obj.CURVE{n,m,o}=curve;
    end
%%
    function adopt_CI(obj,EP2,EP2Inds,EPInds)
        if nargin < 3 || isempty(EP2Inds)
            EP2Inds=obj.ulvls;
        end
        flds={'TCI','bBoot'};
        if nargin < 4
            EPInds=[];
        end
        obj.adopt_fun(EPInds,EP2,EP2Inds,flds);
    end
    function get_thresh(obj)
        obj.selDataPartition('S');
        for i = 1:length(obj.ulvlInds)
            [n,m]=ind2sub(obj.sz,obj.ulvlInds(i));
            obj.T(n,m)=obj.CURVE{n,m}.tFit;
        end
    end
    function adopt_curves(obj,EP2,EP2Inds,EPInds)
        flds={'CURVE','T','X'};
        if nargin < 4
            EPInds=[];
        end
        obj.adopt_fun(EPInds,EP2,EP2Inds,flds);
    end
    function change_CIsz(obj,val)
        for i = 1:size(obj.ulvls)
            N=obj.ulvls(i,1);
            M=obj.ulvls(i,2);
            if isempty(obj.T(N,M));
                continue
            end
            obj.CURVE{N,M}.CIsz=val;
            obj.CURVE{N,M}.pack_boot;


            [n,m]=ind2sub(obj.sz,obj.ulvlInds(i));
            %if ~obj.bFlip
            %    nn=n;
            %    n=m;
            %    m=nn;
            %end
            s=substruct('()',[':' num2cell([n m]) ]);
            obj.TCI=subsasgn(obj.TCI,s,obj.CURVE{N,M}.tCI);
        end
    end
    function adopt_fun(obj,EPInds,EP2,EP2Inds,flds)
        if isempty(EPInds)
            EPInds=EP2Inds;
        end
        for f = 1:length(flds)
            fld=flds{f};
            for i = 1:size(EP2Inds)
                F=EP2.(fld);
                if isempty(F)
                    continue
                end
                inds1=EPInds(i,:);
                inds2=EP2Inds(i,:);
                if strcmp(fld,'TCI')
                    obj.(fld)(:,inds1(1),inds1(2))=Obj.copy(F(:,inds2(1),inds2(2)));
                elseif iscell(F)
                    obj.(fld){inds1(1),inds1(2)}=Obj.copy(F{inds2(1),inds2(2)});
                else
                    obj.(fld)(inds1(1),inds1(2))=Obj.copy(F(inds2(1),inds2(2)));
                end
            end
        end
    end
end
methods(Static,Hidden)
end
end

