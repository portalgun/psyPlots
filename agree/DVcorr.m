classdef DVcorr < handle & DVcorr_plots & DVcorr_plots_p
% PER STANDARD.
% By comparison or for all comparisons
properties
    % input
    modelType
    rhoFix
    mu1Fix
    mu2Fix
    cr1Fix
    cr2Fix
    bRhoFixCmp
    bMu1FixCmp
    bMu2FixCmp
    bCr1FixCmp
    bCr2FixCmp

    nParabSim
    nTrlParabSim
    Magr
    nPssInd2use

    % stat params
    nBoot
    CI
    fitOpts

    % plot params
    bPlot
    colorEllipse
    colormap
    CIcolor
    LineWidth
    Xname
    Xunits
    bCombineCmp

    % fitt dv corr
    RHO
    MU1
    MU2
    CR1
    CR2

    PA
    PD

    RSPAGR
    NTRL
    NEGLL

    % FIT FROM SIMULATIONS
    PAsim
    PCsim
    PAsimL
    PAsimU

    % binomial null model
    PAbino
    PCbino
    PAbinoL
    PAbinoU

    %polynomial fitting
    PCfit
    PCfitU
    PCfitL

    % empirical
    PAe
    PCe
    PAboot
    PCboot

    % FROM CURVES
    % COMBINED PASSES
    T
    TCI
    TSE

    TvarI
    TvarE
    TvarT

    varI
    varE
    varT

    %INVIDIUAL PASSES
    Tind
    TindCI
    TindSE

    boot
end
properties(Hidden = true)
    zer
    colors
    nStd
    nCmp

% input

    stdX
    cmpX
    stdXunq
    cmpXunq
    RcmpChs
    RCMPCHS

% select
    ind

    std
    cmp
    DATA

    rho
    mu1
    mu2
    cr1
    cr2

    rspAgr
    Ntrl %Nall
    negLL
    NtrlAll
    negLLall

    pA
    pD

% fitting
    LB
    UB
    param
    param0

    Npp
    Nnn
    Npn
    Nnp

    Pnn
    Pnp
    Ppn
    Ppp

    nR
    nM
    nC

    ri
    m1i
    m2i
    c1i
    c2i


end
methods
    function obj= DVcorr(stdX2,cmpX2,RcmpChs2,Opts)

        % XXX after dataTable
        bTest=0;
        if ~exist('stdX2','var')
            bTest=1;
        else
            obj.stdX=stdX2;
        end
        if ~exist('cmpX2','var')
            bTest=1;
        else
            obj.cmpX=cmpX2;
        end
        if ~exist('RcmpChs2','var')
            bTest=1;
        else
            obj.RcmpChs=RcmpChs2;
        end

        if ~exist('Opts','var')
            Opts=struct;
        end

        p=inputParser();
        p.addParameter('modelType','RMM');
        p.addParameter('rhoFix'   ,[]);
        p.addParameter('mu1Fix'   ,[]);
        p.addParameter('mu2Fix'   ,[]);
        p.addParameter('cr1Fix'   ,[]);
        p.addParameter('cr2Fix'   ,[]);
        p.addParameter('bRhoFixCmp',0);
        p.addParameter('bMu1FixCmp',0);
        p.addParameter('bMu2FixCmp',0);
        p.addParameter('bCr1FixCmp',0);
        p.addParameter('bCr2FixCmp',0);

        p.addParameter('bPlot'    ,0);
        p.addParameter('nBoot'    ,1000);
        p.addParameter('colorEllipse',[0 0 0]);
        p.addParameter('colormap','cmapGray2');
        p.addParameter('CI',95);
        p.addParameter('CIcolor',0.9.*[1 1 1]);
        p.addParameter('LineWidth',2);

        p.addParameter('minFuncType','fmincon');
        p.addParameter('bParallel',0);
        p.addParameter('nPssInd2use',[1 2]);
        p.addParameter('Magr',[]);
        p.addParameter('nParabSim',1000);
        p.addParameter('nTrlParabSim',[]);
        p.addParameter('Xname',[]);
        p.addParameter('Xunits',[]);


        p=parseStruct(Opts,p);
        flds=fieldnames(p.Results);
        for i = 1:length(flds)
            fld=flds{i};
            if strcmp(fld,'minFuncType') || strcmp(fld,'bParallel')
                continue
            end
            obj.(fld)=p.Results.(fld);
        end
        if isempty(obj.Magr)
            obj.Magr=size(obj.RcmpChs,2);
        end


        % fit opts
        fitOpts = optimset(p.Results.minFuncType);
        fitOpts.Algorithm      = 'active-set';
        fitOpts.LargeScale     = 'off';
        if p.Results.bParallel
           fitOpts.UseParallel = 'always';
        else
            fitOpts.UseParallel = 'never';
        end
        fitOpts.Display        = 'off';
        fitOpts.MaxIter        = 250;
        fitOpts.TolFun         = 1e-8; % 1e-6
        fitOpts.TolX           = 1e-8;
        fitOpts.TolCon         = 1e-8;
        obj.fitOpts=fitOpts;

        if bTest
            obj.DATA = mvnrnd([0 0],[1 .5; .5 1],1000);
        end

        if isempty(obj.nTrlParabSim)
            obj.nTrlParabSim=size(obj.cmpX,1);
        end
        try
            obj.rm_nans();
        end

        obj.cmpXunq=unique(obj.cmpX);
        obj.stdXunq=unique(obj.stdX);
        obj.nCmp = numel(obj.cmpXunq);
        obj.nStd = numel(obj.stdXunq);
        if obj.nStd > 1
            error('more than one standard')
            % TODO write code that combines these?
        end
    end
    function obj=rm_nans(obj)
        ii=any(isnan(obj.RcmpChs),2) | any(isnan(obj.cmpX),2) | any(isnan(obj.stdX),2);
        obj.cmpX(ii,:)=[];
        obj.stdX(ii,:)=[];
        obj.RcmpChs(ii,:)=[];
    end
    function obj=run(obj)
        obj.get_DV_corr();
        obj.get_magr();
        obj.get_binom_parab();
        if numel(unique(obj.RHO)) > 1
            obj.get_fit_parab();
        else
            obj.get_sim_parab();
        end
    end
    function obj = get_param_bounds(obj)

        if obj.bCombineCmp
            tol=1e-2;
        else
            tol=0.05;
        end
        
        switch(obj.modelType)
        case 'R'
            % FIT RHO
            obj.LB      = [-1]+tol;
            obj.UB      = [+1]-tol;
            obj.ri=1;
        case 'RM'
            % FIT RHO,MU1=MU2
            obj.LB      = [-1 -4 ]+tol;
            obj.UB      = [ 1  4 ]-tol;
            obj.ri=1;
            obj.m1i=2;
            obj.m2i=2;
        case 'RC'
            % FIT RHO,Cr1=CR2
            obj.LB      = [-1 -3 ]+tol;
            obj.UB      = [ 1  3 ]-tol;
            obj.ri=1;
            obj.c1i=2;
            obj.c2i=2;
        case 'RCC'
            % FIT RHO,Cr1,CR2
            obj.LB      = [-1 -3 -3]+tol;
            obj.UB      = [ 1  3  3]-tol;
            obj.ri=1;
            obj.c1i=2;
            obj.c2i=3;
        case 'RMC'
            % FIT RHO,MU1=MU2,CR1=CR2
            obj.LB      = [-1 -4 -3]+tol;
            obj.UB      = [ 1  4  3]-tol;
            obj.ri=1;
            obj.m1i=2;
            obj.m2i=2;
            obj.c1i=3;
            obj.c2i=3;
        case 'RMM'
            % FIT RHO,MU1,MU2
            obj.LB      = [-1 -4 -4 ]+tol;
            obj.UB      = [ 1  4  4 ]-tol;
            obj.ri =1;
            obj.m1i=2;
            obj.m2i=3;
        case 'RMCC'
            obj.LB      = [-1 -4 -3 -3]+0.05;
            obj.UB      = [ 1  4  3  3]-0.05;
            obj.ri =1;
            obj.m1i=2;
            obj.m2i=2;
            obj.c1i=3;
            obj.c2i=4;
        case 'RMMC'
            obj.LB      = [-1 -4 -4 -3]+0.05;
            obj.UB      = [ 1  4  4  3]-0.05;
            obj.ri =1;
            obj.m1i=2;
            obj.m2i=3;
            obj.c1i=4;
            obj.c2i=4;
        case 'RMMCC'
            % FIT RHO,MU1,MU2,CR1,CR2
            obj.LB      = [-1 -4 -4 -3 -3]+0.05;
            obj.UB      = [ 1  4  4  3  3]-0.05;
            obj.ri =1;
            obj.m1i=2;
            obj.m2i=3;
            obj.c1i=4;
            obj.c2i=5;
        otherwise
            error(['psyFitDecisionVariableCorrParamBounds: WARNING! unhandled modelType=' num2str(obj.modelType)]);
        end

        if obj.bCombineCmp
            obj.LB=repmat(obj.LB,obj.nCmp,1);
            obj.UB=repmat(obj.UB,obj.nCmp,1);
        end

    end
    function obj=boot_package_init(obj,flds)
        obj.boot=struct();
        if obj.bCombineCmp
            C=1;
        else
            C=obj.nCmp;
        end
        for i = 1:length(flds)
            if ismember(flds{i},{'NTRL'})
                continue
            end
            obj.boot.(flds{i})=zeros(C,obj.nBoot);
        end

    end
    function obj=boot_package(obj,flds,b)
        for i =1:length(flds)
            if isstruct(obj.(flds{i}))
                continue
            end
            obj.boot.(flds{i})(:,b)=obj.(flds{i});
        end
    end
    function obj=boot_complete(obj,flds)
        CIlohi = 0.5*(1-obj.CI/100) + [0 obj.CI/100];
        for i = 1:length(flds)
            if ismember(flds{i},{'NTRL'})
                continue
            end
            fld=flds{i};
            cifld=[fld '_CI'];
            sefld=[fld '_SE'];
            val=obj.boot.(fld)( ~isnan(any(obj.boot.(fld),2)),:);
            obj.boot.(cifld)=quantile(val, CIlohi,2);
            obj.boot.(sefld)=std(val,[],2);
            if endsWith(flds{i},'_CI')
                obj.(fld)=mean(obj.boot.(fld),3);
            else
                obj.boot.(fld)=mean(obj.boot.(fld),2);
            end
        end
    end
    function obj=init_package(obj,C)
        obj.NTRL    = zeros(C,1);

        obj.NEGLL   = zeros(C,1);
        obj.RSPAGR  = cell(C,1);
        obj.PA  = zeros(C,1);
        obj.PD  = zeros(C,1);

        obj.RHO  = zeros(C,1);
        obj.MU1  = zeros(C,1);
        obj.MU2  = zeros(C,1);
        obj.CR1  = zeros(C,1);
        obj.CR2  = zeros(C,1);
    end

    function obj=package(obj,c)
        obj.RSPAGR{c,1}=obj.rspAgr;


        if obj.bCombineCmp
            obj.PA=obj.pA;
            obj.PD=obj.pD;
            obj.NEGLL = obj.negLL;
            obj.NTRL = obj.Ntrl;
            obj.RHO=obj.rho;
            obj.MU1=obj.mu1;
            obj.MU2=obj.mu2;
            obj.CR1=obj.cr1;
            obj.CR2=obj.cr2;
        else

            obj.PA(c,1)    =obj.pA;
            obj.PD(c,1)    =obj.pD;
            obj.NEGLL(c,1) =obj.negLL;
            obj.NTRL(c,1)  =obj.Ntrl;

            obj.RHO(c,1)   =obj.rho;
            obj.MU1(c,1)   =obj.mu1;
            obj.MU2(c,1)   =obj.mu2;
            obj.CR1(c,1)   =obj.cr1;
            obj.CR2(c,1)   =obj.cr2;
        end

    end
    function obj=select_data_only(obj,c)
        IND=obj.cmpX==obj.cmpXunq(c);
        obj.DATA=obj.RcmpChs(IND,:);
    end
    function obj=select(obj,c)
        IND=obj.cmpX==obj.cmpXunq(c);

        obj.DATA=obj.RcmpChs(IND,:);
        obj.std=obj.stdX(IND);
        obj.cmp=obj.cmpX(IND);

        %obj.Ntrl=  obj.NTRL(c,1);
        %obj.negLL= obj.NEGLL(c,1);
        %obj.pA=    obj.PA(c,1);
        %obj.pD=    obj.PD(c,1);
        obj.rho=   obj.RHO(c,1);
        obj.mu1=   obj.MU1(c,1);
        obj.mu2=   obj.MU2(c,1);
        obj.cr1=   obj.CR1(c,1);
        obj.cr2=   obj.CR2(c,1);

        flds=fieldnames(obj.RSPAGR);
        for i = 1:length(flds)
            fld=flds{i};
            obj.rspAgr.(fld)=obj.RSPAGR.(fld)(c);
        end
    end
    function obj= proc_model_type(obj)
        obj.nR=sum(ismember(obj.modelType,'R'));
        obj.nM=sum(ismember(obj.modelType,'M'));
        obj.nC=sum(ismember(obj.modelType,'C'));
    end
    function obj = init_params(obj)
        if obj.bCombineCmp
            val=1/8;
            C=obj.nCmp;
        else
            val=1/4;
            C=1;
        end

        obj.param0=zeros(C,size(obj.LB,2));
        for c = 1:C
            obj.param0(c,:)  = rand_interval_fun([obj.LB(c,:); obj.UB(c,:)]).*val;
        end

        function r=rand_interval_fun(range)
            m=size(range,2);
            n=1;
            p=1;
            if size(range,1)==1 && size(range,2)==2, range = transpose(range);  end
            r = bsxfun(@plus,range(1,:) + 0.5.*(1-p).*(range(2,:)-range(1,:)),    p.*bsxfun(@times,(range(2,:)-range(1,:)),rand(n,m)) );
        end
    end
    function obj=init_fix(obj)
        flds={'rhoFix', 'mu1Fix', 'mu2Fix', 'cr1Fix', 'cr2Fix'};

        for i = 1:length(flds)
            fld=flds{i};
            if ~isempty(obj.(fld)) && size(obj.(fld),1) == 1
                obj.(fld)=repmat(obj.(fld),obj.nCmp,1);
            end
        end
    end
    function obj=get_DV_corr(obj)
        obj.get_bCombineCmp;
        obj.get_param_bounds;

        if obj.bCombineCmp
            obj.init_fix();
            % XXX
            %obj.normalize_model_means();
            C=1;
        else
            C=obj.nCmp;
        end
        obj.init_package(C);

        for c = 1:C

            if obj.bCombineCmp
                obj.DATA=obj.RcmpChs;
            else
                IND=obj.cmpX==obj.cmpXunq(c);
                obj.DATA=obj.RcmpChs(IND,:);
            end

            obj.fit();
            obj.get_rsp_agr();
            [obj.negLL,obj]=obj.get_negLL(obj.param);
            obj.get_final_params();
            obj.package(c);
        end
        obj.RSPAGR=structMerge(obj.RSPAGR{:});
    end
    function obj=get_bCombineCmp(obj)
        flds= {'bRhoFixCmp', 'bMu1FixCmp', 'bMu2FixCmp', 'bCr1FixCmp', 'bCr2FixCmp'};

        obj.bCombineCmp=0;
        for i = 1:length(flds)
            fld=flds{i};
            if obj.(fld)
                obj.bCombineCmp=1;
                break
            end
        end
    end
    function obj=normalize_model_means(obj)
        Rideal = obj.mu1Fix;
        Opts=struct('sFix',1);
        curve=psyCurve(obj.stdX,obj.cmpX,Rideal,Opts);
        obj.mu1Fix = curve.DPfit.*sqrt(.43)./sqrt(2);
    end
    function obj = fit(obj)
        obj.proc_model_type();
        obj.zer=zeros(obj.nCmp,1);

        obj.init_params();

        %param = [obj.mu1 obj.mu2 obj.rho obj.cr1 obj.cr2];
        f=@(p) obj.get_negLL(p);
        [obj.param,~] = fmincon(f,obj.param0,[],[],[],[],obj.LB,obj.UB,[],obj.fitOpts);
    end
    function obj = fit_boot(obj)
        flds={'NTRL','NEGLL','RSPAGR','PA', 'PD', 'RHO', 'MU1', 'MU2', 'CR1', 'CR2'};

        obj.boot_package_init(flds);
        obj.RCMPCHS=obj.RcmpChs;
        for b = 1:obj.nBoot

            % RESAMPLE DATA W. REPLACEMENT
            [~,ii] = datasample(obj.RCMPCHS,size(obj.RCMPCHS,1));
            obj.RcmpChs=obj.RCMPCHS(ii,:);

            obj.get_DV_corr();
            obj.boot_package(flds,b);
        end
        obj.boot_complete(flds);
    end
 %%%%%%%%%%%%%%%%
%% PRIVATE PLOT
   function [negLLall,obj] = get_negLL(obj,param)

        % XXX bSorted = psyDataTrialSortCheck(stdX,cmpX,nPssInd);
        %if bSorted == 0
        %    disp(['psyFitDecisionVariableCorrAllNegLL: WARNING! data not sorted!!!']);
        %end

        % UNPACK PARAMETERS
        % PROBABILITY OF AGREEMENTS AND DISAGREEMENTS

        % DETERMINE WHETHER DATA STORES RESPONSES OR DECISION VARIABLE VALUES
        if obj.bCombineCmp
            C=obj.nCmp;
        else
            C=1;
        end
        for c = 1:C

            if obj.bCombineCmp
                obj.select_data_only(c);
            end
            obj.get_params(param,c);

            % for rho(c,s)
            obj.get_joint_p(c);

            DATAminmax = [min(obj.DATA(:)) max(obj.DATA(:))];
            if DATAminmax(1) == 0 && DATAminmax(2) == 1
                % NUMBER OF EACH TYPE OF RESPONSE
                obj.Npp(c) = sum(obj.DATA(:,1) == 1 & obj.DATA(:,2) == 1);
                obj.Nnn(c) = sum(obj.DATA(:,1) == 0 & obj.DATA(:,2) == 0);
                obj.Npn(c) = sum(obj.DATA(:,1) == 1 & obj.DATA(:,2) == 0);
                obj.Nnp(c) = sum(obj.DATA(:,1) == 0 & obj.DATA(:,2) == 1);
            else
                % COMPARE DECISION VARIALBE TO CRITERION
                obj.Npp(c) = sum(obj.DATA(:,1) >= obj.cr1 & obj.DATA(:,2) >= obj.cr2);
                obj.Nnn(c) = sum(obj.DATA(:,1) <  obj.cr1 & obj.DATA(:,2) <  obj.cr2);
                obj.Npn(c) = sum(obj.DATA(:,1) >= obj.cr1 & obj.DATA(:,2) <  obj.cr2);
                obj.Nnp(c) = sum(obj.DATA(:,1) <  obj.cr1 & obj.DATA(:,2) >= obj.cr2);
            end
        end

        % TOTAL NUMBER OF RESPONSES
        obj.Ntrl = (obj.Npp + obj.Nnn + obj.Npn + obj.Nnp);
        obj.NtrlAll = sum(obj.Ntrl(:));

        % TOTAL NUMBER OF TRIALS

        % NEGATIVE LOG-LIKELIHOOD PER DATA POINT
        obj.negLL = -( obj.Npp.*log(obj.Ppp) + obj.Nnn.*log(obj.Pnn) + obj.Npn.*log(obj.Ppn) + obj.Nnp.*log(obj.Pnp) )./obj.Ntrl; % XXX no frac?
        obj.negLLall = sum(obj.negLL(~isnan(obj.negLL(:))));

        negLLall=obj.negLLall;
    end
    function obj=get_joint_p(obj,c)
        obj.Pnn(c) = mvncdf([obj.cr1 obj.cr2],[obj.mu1 obj.mu2],[1 obj.rho; obj.rho 1]);        % NEG/NEG
        obj.Pnp(c) = mvncdf([obj.cr1 Inf]    ,[obj.mu1 obj.mu2],[1 obj.rho; obj.rho 1]) - obj.Pnn(c);  % NEG/POS
        obj.Ppn(c) = mvncdf([Inf     obj.cr2],[obj.mu1 obj.mu2],[1 obj.rho; obj.rho 1]) - obj.Pnn(c);  % POS/NEG
        obj.Ppp(c) = 1 - obj.Pnn(c) - obj.Ppn(c) - obj.Pnp(c);                               % POS/POS

        obj.pA(c) = obj.Ppp(c) + obj.Pnn(c);
        obj.pD(c) = obj.Ppn(c) + obj.Pnp(c);
    end
    function obj=get_params(obj,param,c)
        if obj.bRhoFixCmp
            r=1;
        else
            r=c;
        end
        if obj.bMu1FixCmp
            m1=1;
        else
            m1=c;
        end
        if obj.bMu2FixCmp
            m2=1;
        else
            m2=c;
        end
        if obj.bCr1FixCmp
            c1=1;
        else
            c1=c;
        end
        if obj.bCr2FixCmp
            c2=1;
        else
            c2=c;
        end

        %RHO
        if ~isempty(obj.rhoFix)
            obj.rho = obj.rhoFix(r);
        elseif obj.nR > 0
            obj.rho=param(r,obj.ri);
        else
            obj.rho=obj.zer(r);
        end

        %Mu1
        if ~isempty(obj.mu1Fix)
            obj.mu1=obj.mu1Fix(m1);
        elseif  obj.nM > 0
            obj.mu1=param(m1,obj.m1i);
        else
            obj.mu1=obj.zer(m1);
        end

        %Mu2
        if ~isempty(obj.mu2Fix)
            obj.mu2=obj.mu2Fix(m2);
        elseif obj.nM > 0
            obj.mu2=param(m2,obj.m2i);
        else
            obj.mu2=obj.zer(m2);
        end

        %Cr1
        if ~isempty(obj.cr1Fix)
            obj.cr1=obj.cr1Fix(c1);
        elseif obj.nC > 0
            obj.cr1=param(c1,obj.c1i);
        else
            obj.cr1=obj.zer(c1);
        end

        %Cr2
        if ~isempty(obj.cr2Fix)
            obj.cr2=obj.cr2Fix(c2);
        elseif obj.nC > 0
            obj.cr2=param(c2,obj.c2i);
        else
            obj.cr2=obj.zer(c2);
        end

    end
    function obj=get_final_params(obj)
        c=1:size(obj.param,1);
        obj.get_params(obj.param,c);
        if obj.bRhoFixCmp
            obj.rho=repmat(obj.rho,obj.nCmp,1);
        end
        if obj.bMu1FixCmp
            obj.mu1=repmat(obj.mu1,obj.nCmp,1);
        end
        if obj.bMu2FixCmp
            obj.mu2=repmat(obj.mu2,obj.nCmp,1);
        end
        if obj.bCr1FixCmp
            obj.cr1=repmat(obj.cr1,obj.nCmp,1);
        end
        if obj.bCr2FixCmp
            obj.cr2=repmat(obj.cr2,obj.nCmp,1);
        end
    end
    function obj = get_rsp_agr(obj)
        obj.rspAgr = struct('Npp',obj.Npp,'Nnn',obj.Nnn,'Npn',obj.Npn,'Nnp',obj.Nnp,...
                            'Ppp',obj.Ppp,'Pnn',obj.Pnn,'Ppn',obj.Ppn,'Pnp',obj.Pnp);
    end
%% MAGR
    function obj=get_magr(obj)

        C=obj.nCmp;
        %n=size(obj.cmpX,1);

        PCboot=zeros(C,obj.nBoot);
        PAboot=zeros(C,obj.nBoot);
        for c = 1:C
            obj.select_data_only(c);
            for b = 1:obj.nBoot
                N=size(obj.DATA,1);
                ii=datasample(1:N,N,'Replace',true);
                R=obj.DATA(ii,:);

                PCboot(c,b)=sum(R(:))./numel(R);
                PAboot(c,b)=sum(sum(diff(R,[],2),2)==0,1)/size(R,1);
            end
        end
        obj.PAboot=PAboot;
        obj.PCboot=PCboot;
        obj.PAe=mean(PAboot,2);
        obj.PCe=mean(PCboot,2);
    end
    function obj=get_fit_parab(obj)

        PC=obj.PCboot;
        PA=obj.PAboot;
        order=10;

        CI = [0.5.*(1-obj.CI/100) 1-0.5.*(1-obj.CI/100)];
        xfix=[-.5 .5];
        yfix=[1 1];


        nBoot=size(PC,2);
        p=cell(nBoot,1);
        o=order;
        while true
            for i = 1:nBoot
                p{i}=polyfix(PC(:,i)-.5,PA(:,i),2,xfix,yfix);
            end
            sz=cellfun(@numel,p);
            if all(diff(sz)==0)
                o=sz(1);
                p=vertcat(p{:});
                break
            else
                b=unique(sz);
                [count]=hist(sz,b);
                o=b(b==max(count));
            end
        end
        Uyi=quantile(p(:,end),CI(1));
        Lyi=quantile(p(:,end),CI(2));
        p=zeros(nBoot,o);
        for i = 1:nBoot
            pU(i,:)=polyfix(PC(:,i)-.5,PA(:,i),o-1,[xfix 0],[yfix Uyi]);
            pL(i,:)=polyfix(PC(:,i)-.5,PA(:,i),o-1,[xfix 0],[yfix Lyi]);
        end

        U=mean(pU,1);
        L=mean(pL,1);
        m=mean([U;L],1);

        obj.PCfit=m;
        obj.PCfitU=U;
        obj.PCfitL=L;
    end
    function obj=get_sim_parab(obj)
        if numel(unique(obj.RHO))==1
            rh=obj.RHO(1);
        else
            rh=mean(obj.RHO);
        end
        cov=[1 rh; rh 1];
        mu=repmat(transpose((-2:0.2:2)),[1 2]);
        [obj.PAsim,obj.PCsim,obj.PAsimL,obj.PAsimU]=obj.get_parab(mu,cov);
    end
    function obj=get_binom_parab(obj)
        cov=[1 0; 0 1];
        mu=repmat(transpose((-2:0.2:2)),[1 2]);
        [obj.PAbino,obj.PCbino,obj.PAbinoL,obj.PAbinoU]=obj.get_parab(mu,cov);
    end
    function [PA,PC,PAl,PAu]=get_parab(obj,mu,COV)
        N=numel(unique(obj.nPssInd2use));
        CI = [0.5.*(1-obj.CI/100) 1-0.5.*(1-obj.CI/100)];

        flag=0;
        if isequal(size(COV),[N N])
            cov=COV;
        elseif size(COV,1)==size(mu,1) && size(COV,2)==1
            flag=1;
        end

        p=pr(length(mu),length(mu)/100,'Getting Parabola');
        for i = 1:length(mu) % FOR EACH MEAN (cmp)
            if flag
                cov=COV(i)*eye(N);
                cov(cov==0)=1;
                cov=fliplr(cov);
            end
            p.u();

            % COMPUTE AGREEMENT AND CMP CHOSEN MONTE CARLO SIMULATIONS
            [PAtmp, PCtmp] = obj.get_parab_mc(N,obj.Magr,mu(i,:),cov,obj.nParabSim,obj.nTrlParabSim);
            % STORE ACTUAL P(A) AND P(C)
            PC(i,:) = mean(PCtmp);
            PA(i,:) = mean(PAtmp);
            % STORE LOWER BOUNDS ON P(A)
            PAl(i,:) = quantile(PAtmp,CI(1));
            PAu(i,:) = quantile(PAtmp,CI(2));
        end
        p.c();
    end
    function [PA,PC]=get_parab_mc(obj,N,M,mu,cov,numSim,numTrlPerSim)
    % simulates 2AFC trials for which the observer's decision variable is
    % Gaussian distributed, then computes proportion agreement and proportion
    % cmp chosen
        if sum(cov(boolean(1-eye(N))))==0 % IF THERE IS NO CORRELATION BETWEEN DECISION VARIABLES, CAN COMPUTE FAST
            % GENERATE DECISION VARIABLES
            Z = normrnd(mu(1),cov(1,1),[numTrlPerSim N numSim]);
            % COMPUTE PROPORTION CMP CHOSEN
            PC = transpose(sum(reshape(Z,[size(Z,1)*size(Z,2) size(Z,3)])>0))./(numTrlPerSim*N);
            % COMPUTE AGREEMENT
            PA = squeeze(sum(sum(Z>=0,2)==M | sum(Z<0,2)==M)./numTrlPerSim);
        else % IF THERE IS CORRELATION BETWEEN DECISION VARIABLES, LOOP (SLOWER)
            for i = 1:numSim % FOR EACH SIMULATION
                % GENERATE DECISION VARIABLES
                Z = mvnrnd(mu,cov,numTrlPerSim);
                % COMPUTE PC
                PC(i,:) = sum(Z(:)>0)./numel(Z);
                % COMPUTE M-AGREEMENT
                bZmoreThan0 = Z>=0;
                bZlessThan0 = Z<0;
                PA(i,:) = sum(sum(bZmoreThan0,2)==M | sum(bZlessThan0,2)==M)./numTrlPerSim;
            end
        end
    end
    function obj=get_noise_sources(obj,PsyCurve)
        %C=obj.nCmp;
        %obj.varE=zeros(C,1);
        %obj.varI=zeros(C,1);
        %obj.varT=zeros(C,1);


        %for c = 1:C
        %    obj.select(c);
        %    dp=PsyCurve.DPfit(c);
        %    if obj.cmp==obj.std
        %        continue
        %    end
        %    [obj.varE(c),obj.varI(c),obj.varT(c)]=DVcorr.getNoiseSourcesInd(obj.cmp,obj.std,obj.rho,dp);
        %end


        dp=PsyCurve.DPcrt;
        T=PsyCurve.tFit;
        rho=obj.rho;
        [obj.TvarE,obj.TvarI,obj.TvarT]=DVcorr.getNoiseSources(T,rho,dp);
        obj.T=T;

        t=T-obj.tSE;
        [obj.TvarE, obj.TvarI, obj.TvarT]=DVcorr.getNoiseSources(t,rho,dp);

        t=T+obj.tSE;
        [obj.TvarE, obj.TvarI, obj.TvarT]=DVcorr.getNoiseSources(t,rho,dp);

        t=obj.tCI(1);
        [obj.TvarE, obj.TvarI, obj.TvarT]=DVcorr.getNoiseSources(t,rho,dp);

        t=obj.tCI(2);
        [obj.TvarE, obj.TvarI, obj.TvarT]=DVcorr.getNoiseSources(t,rho,dp);
    end
end
methods(Static=true)
    function [varE,varI,varT]=getNoiseSourcesInd(cmp,std,rho,dp)
        if numel(unique(cmp))==1 && numel(unique(std))==1
            cmp=cmp(1);
            std=std(1);
        end
        varT=((cmp-std)/dp).^2;
        varE=rho*varT;
        varI=varT-varE;
    end
    function [TvarE,TvarI,TvarT]=getNoiseSources(T,rho,dp)
        TvarT=T^2;
        TvarE=rho*TvarT;
        TvarI=TvarT-TvarE;

    end
end
end
