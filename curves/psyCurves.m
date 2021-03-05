classdef psyCurves < handle
% TODO
% tight plots
properties
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
    minFuncType

    nTrials
    nLvls
    nTrialsPerLvl
    DT

    markersize
    LineWidth
    Xname
    Xunits
    bPlot
    bPlotCI

    psyXtickFormat
    psyYtickFormat
    psyXtickAngle
    psyYlim
    psyXlim

    threshXtickFormat
    threshYtickFormat
    threshXtickAngle
    threshYlim
    threshXlim

    nX
    lcolors
    bExtendCurve % whether to extend curve
    colormap
    colors
    markerfaces
    shapes
    CIalpha

    color
    markerface
    shape

    curs

end
properties(Hidden=true)
    Opts % for individual fits

    X % n x m x 3 ([std cmp chs ])
    %n std subj x 3
    % diff subplot plot per subject
    M % 1 x m x 4 ([subj exp prj pass])

    n
    m
    nSubj
    nStd

    subjs
    EXPS
    titl

    % selection
    ind

    cmpX
    stdX
    RcmpChs

    cur

    tCI
    tMU
    tSE
    tX

    bTest
    pr

    continueflag
end
methods
    function obj=psyCurves(X,subjs,titl,Opts)
        % XXX after dataTable
        if ~exist('X','var')
            bTest=1;
        elseif isa(X,'dataTable')
            obj.DT=X;
            obj.proc_dataTable();
            bTest=0;
        elseif isstruct(X)
            obj.DT=X;
            obj.proc_struct();
            bTest=0;
        else
            obj.X=X;
            bTest=0;
        end
        if exist('subjs','var') && isa(subjs,'exp_detail')
            obj.EXPS=subjs;
        elseif exist('subjs','var')
            obj.subjs=subjs;
        end
        if exist('titl','var')
            obj.titl=titl;
        end

        if ~exist('Opts','var')
            Opts=struct;
        end

        %%
        obj.parse_opts(Opts);



        if bTest
            obj.bTest=1;
        end
        if obj.bTest
            obj.gen_testX();
        end

        %%

        if ndims(obj.X) == 3
            sz=size(obj.X);
            obj.X=reshape(obj.X,sz(1),sz(2),1,sz(3));
        end

        obj.n=size(obj.X,1);
        obj.nStd=size(obj.X,2);
        obj.nSubj=size(obj.X,3);
        obj.m=obj.nStd*obj.nSubj;

        obj.generate_colors();
    end
    function obj=parse_opts(obj,Opts)
        bIgnore=0;
        bQuiet=0;
        if ~exist('Opts','var')
            Opts=obj;
            bIgnore=1;
            bQuiet=1;
        end
        p=inputParser();
        p.addParameter('mFix',[]);
        p.addParameter('sFix',[]);
        p.addParameter('bFix',1);
        p.addParameter('DPcrt',1); % 1.36? XXX
        p.addParameter('nIntrvl',2);
        p.addParameter('nBoot',1000);
        p.addParameter('CIsz',95);
        p.addParameter('flipCmpChs',0)

        p.addParameter('psyXtickFormat','%.02f');
        p.addParameter('psyYtickFormat','%.02f');
        p.addParameter('psyXtickAngle',30);
        p.addParameter('psyYlim',[0 1]);
        p.addParameter('psyXlim','15');

        p.addParameter('threshXtickFormat','%.02f');
        p.addParameter('threshYtickFormat','%.02f');
        p.addParameter('threshXtickAngle',30);
        p.addParameter('threshYlim',[0 1]);
        p.addParameter('threshXlim','15');

        p.addParameter('prcntUse',100);
        p.addParameter('minFuncType','fmincon');
        p.addParameter('bTest',0);
        p.addParameter('bPlot',0);
        p.addParameter('bPlotCI',0);
        p.addParameter('LineWidth',2);
        p.addParameter('markersize',12);
        p.addParameter('Xname','X');
        p.addParameter('Xunits','');
        p.addParameter('bExtendCurve',1);
        p.addParameter('nX',200);

        p.addParameter('markerface','w');
        p.addParameter('color','k');
        p.addParameter('shape','o');
        p.addParameter('CIalpha',.4);
        p.addParameter('colormap','cmapGray');

        % XXX
        p.addParameter('lcolors','k');
        p.addParameter('markerfaces','w');
        p.addParameter('colors','k');
        p.addParameter('shapes','o');

        p=parseStruct(Opts,p,bIgnore,bQuiet);
        flds=fieldnames(p.Results);

        for i = 1:length(flds)
            fld=flds{i};
            obj.(fld)=p.Results.(fld);
        end
        obj.Opts=p.Results();
        obj.Opts=rmfield(obj.Opts,'psyXtickFormat');
        obj.Opts=rmfield(obj.Opts,'psyYtickFormat');
        obj.Opts=rmfield(obj.Opts,'psyXtickAngle');
        obj.Opts=rmfield(obj.Opts,'psyYlim');
        obj.Opts=rmfield(obj.Opts,'psyXlim');

        obj.Opts=rmfield(obj.Opts,'threshXtickFormat');
        obj.Opts=rmfield(obj.Opts,'threshYtickFormat');
        obj.Opts=rmfield(obj.Opts,'threshXtickAngle');
        obj.Opts=rmfield(obj.Opts,'threshYlim');
        obj.Opts=rmfield(obj.Opts,'threshXlim');

        obj.Opts=rmfield(obj.Opts,'colormap');
        obj.Opts=rmfield(obj.Opts,'lcolors');
        obj.Opts=rmfield(obj.Opts,'colors');
        obj.Opts=rmfield(obj.Opts,'shapes');
        obj.Opts=rmfield(obj.Opts,'markerfaces');
        obj.Opts=rmfield(obj.Opts,'bPlot');
        obj.Opts.bPlot=0;
    end


    function obj=generate_colors(obj)
        if ischar(obj.colormap)
            cmp=eval(obj.colormap);
        elseif isnumeric(obj.colormap)
            cmp=obj.colormap;
        end
        inds=round(linspace(1,size(cmp,1),obj.nStd));
        colors=cmp(inds,:);
        obj.lcolors=colors;
        obj.colors=colors;
    end
%% fit
    function obj=proc_dataTable(obj)
        iflds=fieldnames(obj.DT.index);
        flds=obj.DT.labels{1};
        if ismember('stdX',iflds) && ~ismember('stdX',flds)
            obj.DT.index_to_fld('stdX');
        end
        n=obj.DT.keyInd;
        m=prod(obj.DT.nLabelVals(2:end));
        obj.X=zeros(n,m,3);

        obj.DT.select('cmpX');
        obj.X(:,:,2)=obj.DT.dt;
        obj.X(:,:,3)==obj.DT.dt;

        % NOTE IF DOESNT WORK, STANDARDIZE USING dataTable
        obj.DT.select('RcmpChs');
        obj.DT.select('stdX');
        obj.X(:,:,1)=obj.DT.dt;
    end
    function obj=proc_struct(obj)
        stdsz=size(obj.DT.stdX,1);
        sz=size(obj.DT.cmpX,1);
        if isequal(stdsz,sz)
            s=obj.DT.stdX;
        elseif ~isequal(stdsz,sz) && mod(sz,stdsz)==0
            n=sz/stdsz;
            s=repmat(obj.DT.stdX,n,1,1);
        elseif mod(sz,stdsz)~=0
            error('std size is incompatibile')
        end
        if isfield(obj.DT,'RcmpChosen')
            fld='RcmpChosen';
        elseif isfield(obj.DT,'RcmpChs')
            fld='RcmpChs';
        elseif isfield(obj.DT,'Rchs')
            fld='Rchs';
        elseif isfield(obj.DT,'responses')
            fld='RcmpChs';
            obj.DT.RcmpChs=obj.DT.cmpIntrvl~=obj.DT.responses;
        end
        c=obj.DT.cmpX;
        r=obj.DT.(fld);
        obj.X=cat(4,s,c,r);
    end
    function obj=fit_all(obj)
        obj.tCI=zeros(obj.m,2);
        obj.pr=prog(obj.m,'fit PsyCurves',1);
        obj.curs=cell(obj.m,1);
        for ind = 1:obj.m
            obj.ind=ind;

            obj.select_col();
            if obj.continueflag
                continue
            end
            % obj.get_col_opts() XXX
            obj.fit_col();
            obj.package_cur();
            obj.get_T_ind();
            obj.pr=obj.pr.update();
        end
    end
    function obj=reshape_out(obj)
        if obj.nSubj > 1
            obj.tMU=reshape_fun(obj.tMU);
            obj.tX=reshape_fun(obj.tX);
            obj.tSE=reshape_fun(obj.tSE);
            obj.tCI=reshape(transpose(obj.tCI),2,obj.nStd,obj.nSubj);
            obj.tCI=permute(obj.tCI,[2,1,3]);
        end
        function out=reshape_fun(in)
            out=reshape(in,obj.nStd,obj.nSubj);
        end
    end
    function obj=fit_col(obj,ind)
        if exist('ind','var') && ~isempty(ind)
            obj.ind=ind;
        end
        obj.cur=psyCurve(obj.stdX,obj.cmpX,obj.RcmpChs,obj.Opts);
    end
%% select
    function obj=select_col(obj,ind)
        if exist('ind','var') && ~isempty(ind)
            obj.ind=ind;
        end
        [n,m]=ind2sub([obj.nStd obj.nSubj],obj.ind);
        X=squeeze(obj.X(:,n,m,:));
        X(any(isnan(X), 2), :) = [];

        obj.stdX=X(:,1);
        obj.cmpX=X(:,2);
        obj.RcmpChs=X(:,3);
        if all(isnan(obj.RcmpChs)) || isempty(obj.RcmpChs)
            obj.continueflag=1;
        else
            obj.continueflag=0;
        end
    end
    function obj=select_cur(obj,ind)
        if exist('ind','var') && ~isempty(ind)
            obj.ind=ind;
        end
        obj.cur=obj.curs{obj.ind};
    end
%% T
    function obj = get_T_all(obj)
        for ind = 1:obj.m
            obj.ind=ind;
            obj.select_cur();
            obj.get_T_ind();
        end
    end
    function obj = get_T_ind(obj,ind)
        if exist('ind','var') && ~isempty(ind)
            obj.ind=ind;
            obj.select.cur();
        end
        obj.tCI(obj.ind,:)=obj.cur.tCI;
        obj.tMU(obj.ind)=obj.cur.tMU;
        obj.tSE(obj.ind)=obj.cur.tSE;
        obj.tX(obj.ind)=obj.cur.stdX(1);
    end
%% package
    function obj=package_cur(obj,ind)
        if exist('ind','var') && ~isempty(ind)
            obj.ind=ind;
        end
        obj.curs{obj.ind}=obj.cur;
    end
%% opts
    function Opts=get_psy_plot_opts(obj)
        Opts=struct();
        Opts.xtickFormat=obj.psyXtickFormat;
        Opts.ytickFormat=obj.psyYtickFormat;
        Opts.xtickAngle =obj.psyXtickAngle;
        Opts.ylim       =obj.psyYlim;
        Opts.xlim       =[];
        Opts.ax         ='auto';
    end
    function obj=get_col_opts(obj,ind)
        if exist('ind','var') && ~isempty(ind)
            obj.ind=ind;
        end
        [ind,~]=ind2sub([obj.nStd,obj.nSubj],obj.ind);

        if numel(obj.lcolors) <=1
            obj.Opts.lcolor=obj.lcolors;
        else
            obj.Opts.lcolor=obj.lcolors(ind,:);
        end
        if numel(obj.markerfaces) <=1
            obj.Opts.markerface=obj.markerfaces;
        else
            obj.Opts.markerface=obj.markerfaces{ind,:};
        end
        if numel(obj.colors) <=1
            obj.Opts.color=obj.colors;
        else
            obj.Opts.color=obj.colors(ind,:);
        end
        if numel(obj.shapes) <=1
            obj.Opts.shape=obj.shapes;
        else
            obj.Opts.shape=obj.shapes(:,ind);
        end
        obj.Opts.bTitle=0;
    end
%% PLOT PUBLIC
    function [] = plot_curves(obj,SUBJS,opts)
        if ~exist('opts','var')
            opts=struct();
        end
        if ~exist('SUBJS','var')
            SUBJS=[];
        end
        figure(nFn);
        obj.plot_curve_all_p(SUBJS,opts);
    end
    function [] = plot_thresh(obj)
        figure(nFn)
        obj.plot_thresh_all_p();
    end

end
%% Plot Private
methods(Hidden=true)
%% plot parts
    function obj = plot_curve_all_p(obj,SUBJS,opts)
        if ~exist('opts','var') || isempty(opts)
            opts=struct();
        end
        opts

        obj.generate_colors();
        obj.parse_opts();

        if isfield(opts,'xtitl')
            xtitl=opts.xtitl;
        else
            xtitl=obj.EXPS.get_xtitl();
        end
        ytitl='Proportion Cmp. Chosen';

        if isfield(opts,'titl')
            titl=opts.titl;
        else
            titl=obj.EXPS.get_titl;
        end
        rtitl=obj.EXPS.subjs;
        rtitl
        ctitl=[];


        Opts=obj.get_psy_plot_opts();
        Opts=structCombinePrefer(opts,Opts);

        Opts.xticks=unique(round(obj.tX,2));
        Opts.ax='square';

        bSUBJS= exist('SUBJS','var') && ~isempty(SUBJS);
        if ~bSUBJS
            r=obj.nSubj;
        elseif iscell(SUBJS) || isnumeric(SUBJS)
            r=length(SUBJS);
            rtitl=rtitl(SUBJS);
        end
        r
        if ~isfield(Opts,'syms')
            syms=repmat('o',r);
        else
            syms=Opts.syms;
            Opts=rmfield(Opts,'syms')
        end
        Opts.position=[3 3 1060 350*r+200];

        sp=subPlots([r 1],xtitl,ytitl,titl,rtitl,ctitl,Opts);

        Ymin=[];
        Ymax=[];
        Xmin=[];
        Xmax=[];
        for ind = 1:obj.m
            [std,subj]=ind2sub([obj.nStd obj.nSubj],ind);
            if bSUBJS && ~ismember(subj,SUBJS)
                continue
            elseif bSUBJS
                subj=find(SUBJS==subj);
            end
            obj.ind=ind;
            obj.select_cur();
            if isempty(obj.cur)
                continue
            end
            obj.get_col_opts();
            sp.select(subj,1);
            hold on
            obj.plot_curve_p([],syms(subj));
            hold off
            Ymin=min([Ymin obj.cur.Ymin]);
            Ymax=max([Ymax obj.cur.Ymax]);
            Xmin=min([Ymin obj.cur.Xmin]);
            Xmax=max([Ymax obj.cur.Xmax]);
        end
        %sp.xlim=[Xmin Xmax];
        %sp.ylim=[Ymin Ymax];
        sp.finalize();

    end
    function obj = plot_curve_p(obj,ind,sym)
        if exist('ind','var') && ~isempty(ind)
            obj.ind=ind;
            obj.select.cur();
        end
        obj.cur.lcolor=obj.Opts.color;
        obj.cur.plot_boot_curve(sym);
        hold off
    end
%% THRESH
    function obj=plot_thresh_all_p(obj,bAverage,bSame,opts)
        if ~exist('opts','var') || isempty(opts)
            opts=struct();
        end
        if ~exist('bSame','var') || isempty(bSame)
            bSame=0;
        end
        if ~exist('bAverage','var') || isempty(bAverage)
            bAverage=0;
        end
        obj.generate_colors();
        obj.parse_opts();

        if ~isempty(obj.EXPS)
            ytitl=['Threshold (' obj.EXPS.Xunits ')'];
            xtitl=obj.EXPS.get_xtitl();

            titl=obj.EXPS.get_titl;
            if bAverage
                rtitl=[];
            else
                rtitl=obj.EXPS.subjs;
            end
            if ~bAverage & bSame
                rtitl=join(rtitl,' ');
                rtitl=rtitl{1};
            end
        else
            xtitl=[];
            rtitl=[];
            ytitl=[];
            titl=[];
        end
        ctitl=[];
        Opts=struct();
        Opts.xticks=unique(round(obj.tX,2));
        if bSame
            Opts.position=[3 3 1060 1060];
        else
            Opts.position=[3 3 1060 350*obj.nSubj+200];
        end
        Opts=structCombinePrefer(opts,Opts);
        if isfield(Opts,'xtitl')
            xtitl=Opts.xtitl;
        end
        if isfield(Opts,'titl')
            titl=Opts.titl;
        end
        if isfield(Opts,'syms')
            syms=Opts.syms;
            Opts=rmfield(Opts,'syms');
        else
            syms='o';
        end
        if isfield(opts,'colors')
            colors=Opts.colors;
            Opts=rmfield(Opts,'colors');
        else
            colors='k';
            if ~bAverage
                repmat(colors,obj.nSubj,1);
            end
        end

        if bSame || bAverage
            RC=[1 1];
        else
            RC=[obj.nSubj 1];
        end
        sp=subPlots(RC,xtitl,ytitl,titl,rtitl,ctitl,Opts);

        if bAverage
            sp.select(1,1);
            obj.plot_thresh_subj_avg_stdev_p(colors(1));
            hold on
            obj.plot_thresh_subj_avg_p(syms);
            sp.finalize();
            return
        end

        for j = 1:obj.nSubj
            if j==1 || ~bSame
                sp.select(j,1);
            end
            if size(colors,1) == obj.nSubj
                color=colors(j,:);
            else
                color=colors(1,:)
            end
            if numel(syms) == obj.nSubj
                sym=syms(j);
            else
                sym=syms(1);
            end
            obj.plot_thresh_CI_p(j,color);
            hold on;
            obj.plot_thresh_p(j,sym);

        end
        sp.finalize();
    end
    function obj = plot_thresh_subj_avg_p(obj,sym)
        if ~exist('sym','var') || isempty(sym)
            sym=obj.shape;
        end

        tX=reshape(obj.tX,obj.nStd,obj.nSubj);
        tMU=reshape(obj.tMU,obj.nStd,obj.nSubj);
        tMu=mean(tMU,2);
        tX=mean(tX,2);

        plot(tX,tMu, [obj.color sym],'markerface', obj.markerface, 'markersize',obj.markersize, 'LineWidth',obj.LineWidth);
        hold off
    end
    function obj = plot_thresh_subj_avg_stdev_p(obj,color)
        if ~exist('color','var') || isempty(color)
            color=[.5 .5 5]
        end
        tX=reshape(obj.tX,obj.nStd,obj.nSubj);
        tMU=reshape(obj.tMU,obj.nStd,obj.nSubj);
        tX=mean(tX,2);
        tMu=mean(tMU,2);
        tStd=std(tMU,[],2);

        U=tMu+tStd;
        L=tMu-tStd;
        interv(tX,U,L,color ,'FaceAlpha',obj.CIalpha,'EdgeColor','none');
    end
    function obj = plot_thresh_p(obj,subj,sym)
        if ~exist('sym','var') || isempty(sym)
            sym=obj.shape;
        end

        tX=reshape(obj.tX,obj.nStd,obj.nSubj);
        tMU=reshape(obj.tMU,obj.nStd,obj.nSubj);

        plot(tX(:,subj),tMU(:,subj), [obj.color sym],'markerface', obj.markerface, 'markersize',obj.markersize, 'LineWidth',obj.LineWidth);
        hold off
    end
    function obj = plot_thresh_err(obj,subj,CI,color)
        if ~exist('color','var') || isempty(color)
            color='k';
        end

        tCI=reshape(    CI,obj.nStd,obj.nSubj,2);
        tX =reshape(obj.tX ,obj.nStd,obj.nSubj);
        
        ymin=transpose(tCI(:,subj,1));
        ymax=transpose(tCI(:,subj,2));
        x=transpose(tX(:,subj));

        x=[x    fliplr(x)];
        y=[ymin fliplr(ymax)];

        patch(x,y, color,'FaceAlpha',obj.CIalpha,'EdgeColor','none');
        hold off
    end
    function obj = plot_thresh_CI_p(obj,subj,color)
        if ~exist('color','var') || isempty(color)
            color='k';
        end
        obj.plot_thresh_err(subj,obj.tCI,color);
    end
    function obj = plot_thresh_SE_p(obj)
        obj.plot_thresh_err(subj,obj.tSE);
    end
%% TESTS
    function obj=gen_testX(obj)
        obj.gen_testX1();
        obj.gen_testX2();
        obj.gen_testX3();
        obj.gen_testX4();
        obj.gen_testX5();
        obj.gen_testX6();
    end
    function obj=gen_testX1(obj)
        %basic
        nTrlPerLvl = 100;
        lvls=-2:.5:2;
        stdX = 0;

        obj.gen_test_helper(nTrlPerLvl,lvls,stdX);
    end
    function obj=gen_testX2(obj)
        % more between levels
        nTrlPerLvl = 100;
        lvls=-2:.25:2;
        stdX = 0;

        obj.gen_test_helper(nTrlPerLvl,lvls,stdX);
    end
    function obj=gen_testX3(obj)
        % more trials
        nTrlPerLvl = 200;
        lvls=-2:.5:2;
        stdX = 0;

        obj.gen_test_helper(nTrlPerLvl,lvls,stdX);
    end
    function obj=gen_testX4(obj)
        % more trials & between levels
        nTrlPerLvl = 200;
        lvls=-2:.25:2;
        stdX = 0;

        obj.gen_test_helper(nTrlPerLvl,lvls,stdX);
    end
    function obj=gen_testX5(obj)
        % Different standard 1
        nTrlPerLvl = 200;
        lvls=4:.5:8;
        stdX = 6;

        obj.gen_test_helper(nTrlPerLvl,lvls,stdX);
    end
    function obj=gen_testX6(obj)
        % Different standard and larger sep
        nTrlPerLvl = 200;
        lvls=11:1:19;
        stdX = 15;

        obj.gen_test_helper(nTrlPerLvl,lvls,stdX);
    end

    function nX=gen_test_helper(obj,nTrlPerLvl,lvls,stdX)
        cmpX = repmat(lvls,nTrlPerLvl,1);
        cmpX = cmpX(:);
        RcmpChs = binornd(1,repmat(normcdf(transpose(unique(cmpX))),nTrlPerLvl,1));
        RcmpChs = RcmpChs(:);
        if numel(stdX)==1
            stdX=repmat(stdX,size(cmpX,1),1);
        end
        nX=[stdX cmpX RcmpChs];
        nX=reshape(nX,size(nX,1),1,3);

        obj.X=appendColIsz(obj.X,nX);
    end
end
end
