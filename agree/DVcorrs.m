classdef DVcorrs < handle & DVcorrs_plots_p & DVcorrs_plots
properties
    cell     % individual
    cellSubj % combined subjects
    EXPS=struct();

    stdXall    % ntrl x ncmp x nsubj x npass
    cmpXall    % ntrl x ncmp x nsubj x npass
    RcmpChsAll % ntrl x ncmp x nsubj x npass

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
    Magr % XXX
    %nobj.nSubjPssInd2use XXX?
    nPssAtOnce

    nBoot
    CI
    fitOpts

    bPlot
    colorEllipse
    colormap
    CIcolor
    LineWidth

    bFlipRC


    bSaveFig
    figSaveTypes={'.eps','.fig','.png'};
end
properties(Hidden=true)
    nStd
    nSubj
    nPss
    colors

    stdX
    cmpX
    RcmpChs

    PAbino
    PCbino
    PAbinoL
    PAbinoU
end
methods
    function obj=DVcorrs(S_or_stdXall,EXPS_or_cmpXall,RcmpChsAll,Opts,EXPS)
    % DVcorrs(S_or_stdXall,EXPS_or_cmpXall,RcmpChsAll,Opts,EXPS)
    % S_or_stdXall = dt OR merged dt or stdXall matrix
        if exist('EXPS','var')
           obj.EXPS=EXPS;
        end
        if ~exist('Opts','var')
            Opts=struct;
        end
        if iscellstruct(S_or_stdXall)
            S=structMergeNanFill(S_or_stdXall{:});
            obj.proc_struct(S);
        elseif isstruct(S_or_stdXall)
            obj.proc_struct(S_or_stdXall);
        else
            obj.stdXall=S_or_stdXall;
            obj.RcmpChsAll=RcmpChsAll;
        end
        if isa(EXPS_or_cmpXall,'exp_details')
            obj.EXPS=EXPS_or_cmpXall;
        else
            obj.cmpXall=EXPS_or_cmpXall;
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
        p.addParameter('LineWidth',2);
        p.addParameter('CI',95);
        p.addParameter('minFuncType','fmincon');
        p.addParameter('bParallel',0);
        %p.addParameter('nPssInd2use',[1 2]);
        p.addParameter('CIcolor',0.9.*[1 1 1]);
        p.addParameter('nPssAtOnce',2);
        p.addParameter('bFlipRC',1);

        p=parseStruct(Opts,p);
        flds=fieldnames(p.Results);
        for i = 1:length(flds)
            fld=flds{i};
            if strcmp(fld,'minFuncType') || strcmp(fld,'bParallel')
                continue
            end
            obj.(fld)=p.Results.(fld);
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
        fitOpts.TolFun         = 1e-8;
        fitOpts.TolX           = 1e-8;
        fitOpts.TolCon         = 1e-8;
        obj.fitOpts=fitOpts;

        obj.nStd=size(obj.RcmpChsAll,2);
        obj.nSubj=size(obj.RcmpChsAll,3);
        obj.nPss=size(obj.RcmpChsAll,4);
    end
    function [IND]=get_ind_cmp(obj)
        PASSES=obj.get_passes();
        nCmp=obj.get_nCmp();

        IND=Set.distribute(1:obj.nStd,1:obj.nSubj,1:size(PASSES,1),1:nCmp);
        ind=repelem(transpose(1:numel(obj.cell)),size(IND,1)/numel(obj.cell),1);
        IND=[IND ind];
    end
    function nCmp=get_nCmp(obj)
        n=zeros(numel(obj.cell),1);
        for i = 1:numel(obj.cell)
            n(i)=obj.cell{i}.nCmp;
        end
        nCmp=max(n);
    end
    function cmps =get_cmps_all(obj)
        nCmp=obj.get_nCmp();
        cmps=nan(numel(obj.cell),nCmp);
        for i = 1:numel(obj.cell)
            vals=transpose(obj.cell{i}.cmpXunq);
            cmps(i,1:length(vals))=vals;
        end
    end
    function stds=get_stds_all(obj)
        stds=nan(numel(obj.cell),1);
        for i = 1:numel(obj.cell)
            stds(i,1)=transpose(obj.cell{i}.stdunq);
        end
    end
    function [stds,cmps]=get_stds_cmps_unq(obj)
        cmps=obj.get_cmps_all;
        stds=obj.get_stds_all;
        both=unique([stds, cmps],'rows');
        stds=both(:,1);
        cmps=both(:,2:end);
    end
    function PASSES=get_passes(obj)
        PASSES=nchoosek([1:obj.nPss],obj.nPssAtOnce);
    end
    function [IND]=get_ind_ind(obj)
        PASSES=obj.get_passes();
        IND=Set.distribute(1:obj.nStd,1:obj.nSubj,1:size(PASSES,1));
    end
    function [IND]=get_ind_subj(obj)
        PASSES=obj.get_passes();
        IND=Set.distribute(1:obj.nStd,1:size(PASSES,1));
    end
    function [val,valid]=get_val_all(obj,name,cellName,varargin)
        % varargin = passses
        % RHO [std x cmp]
        N=size(obj.(cellName),2);

        if isempty(varargin)
            passes=1:obj.nPss;
        else
            passes=[varargin{:}];
        end

        PASSES=obj.get_passes();
        switch cellName
            case 'cellSubj'
                [IND]=obj.get_ind_subj();
            case 'cell'
                [IND]=obj.get_ind_ind();
        end

        for i = 1:size(IND,1)
            ind=IND(i,end);
            if any(~ismember(passes,PASSES(ind,:)))
                continue
            end
            valid=i;
            val{i}=obj.(cellName){i}.(name);
        end
        if isstruct(val{1})
            val=transpose({val{:}});
        else
            val=[val{:}];
        end
        % std x cmp
    end
    function [boot,valid]=get_boot_all(obj,cellName,varargin)
        [boot,valid]=obj.get_val_all('boot',cellName,varargin{:});
        boot=structMerge(boot{:});

        flds=fieldnames(boot);
        for i = 1:length(flds)
            fld=flds{i};
            if size(boot.(fld),3) > 1
                boot.(fld)=permute(boot.(fld),[3,1,2]);
            elseif  size(boot.(fld),2) > 1
                boot.(fld)=permute(boot.(fld),[2,1]);
            end
        end
    end
    function [RHO,valid]=get_rho_all(obj,cellName,varargin)
        [RHO,valid]=obj.get_val_all('RHO',cellName,varargin{:});
    end
    function [VARI,VARE,valid,STD]=get_var_src_all(obj,cellName,varargin)
        % varargin = passses
        % RHO [std x cmp]
        N=size(obj.(cellName),2);

        if isempty(varargin)
            passes=1:obj.nPss;
        else
            passes=[varargin{:}];
        end

        PASSES=obj.get_passes();
        switch cellName
            case 'cellSubj'
                [IND]=obj.get_ind_subj();
            case 'cell'
                [IND]=obj.get_ind_ind();
        end

        for i = 1:size(IND,1)
            ind=IND(i,end);
            if any(~ismember(passes,PASSES(ind,:)))
                continue
            end
            valid=i;
            VARI{i}=sqrt(obj.(cellName){i}.TvarI);
            VARE{i}=sqrt(obj.(cellName){i}.TvarE);
            STD(i)=IND(i,1);
        end
        VARI=transpose([VARI{:}]);
        VARE=transpose([VARE{:}]);
        % std x cmp
    end
    function obj=proc_struct(obj,S)
        n=size(S.cmpX,1);
        m=size(S.stdX,1);
        if m < n && mod(n,m)==0
            r=n/m;
            S.stdX=repmat(S.stdX,r,1,1,1);
        end
        % NOTE IF DOESNT WORK, standardize using dataTable
        obj.stdXall=S.stdX;
        obj.cmpXall=S.cmpX;
        obj.RcmpChsAll=S.RcmpChs;
    end
%% ind
    function obj=run(obj)
        obj.run_ind(0);
        obj.run_subj_comb();
    end
    function obj=run_boot(obj)
        obj.run_ind(1);
    end
    function obj=init_ind(obj)
        IND=obj.get_ind_ind();
        PASSES=obj.get_passes();
        obj.cell=cell(size(IND,1),1);
        Opts=obj.get_Opts('cell');
        for i = 1:size(IND,1)
            std=IND(i,1);
            subj=IND(i,2);
            ind=IND(i,3);
            passes=PASSES(ind,:);

            obj.select_ind(std,subj,passes);

            obj.cell{i}=DVcorr(obj.stdX,obj.cmpX,obj.RcmpChs,Opts);
        end
    end
    function Opts=get_Opts(obj,cellName)
        flds=fieldnames(obj);
        Opts=struct();
        dummy=DVcorr(1,1,1);
        badlist={'fitOpts'};
        for i = 1:length(flds)
            fld=flds{i};
            if ~isprop(dummy,fld) || ismember(fld,badlist)
                continue
            end
            Opts.(fld)=obj.(fld);
        end
    end
    function obj=get_noise_sources(obj,P)
        %DT - trl, std, subj, pass

        % std, subj, ind
        IND=obj.get_ind_ind();
        PASSES=obj.get_passes();
        subjs=repelem([1:obj.nSubj],1,obj.nStd);
        for i = 1:size(IND,1)
            std=IND(i,1);
            subj=IND(i,2);
            ind=IND(i,3);
            passes=PASSES(ind,:);

            imgDim=obj.EXPS.imgDim;
            %obj.EXPS.AimgDim

            % match DIMENSION & PASSES
            for j = 1:numel(P.passInds)
                PimgDim=P.passInds{j}.EXPS.imgDim;
                Ppasses=P.passInds{j}.EXPS.pass;

                if strcmp(imgDim, PimgDim) && all(ismember( passes, Ppasses ))

                    break
                end
            end
            obj.select_ind(std,subj,passes);
            stdX=unique(obj.stdX);
            % STD SUBJ
            p=P.passInds{j};
            subjs==subj;
            k= find((p.tX==stdX) & (subjs==subj));


            % COMBINED PASSES
            T=p.tMU(k);
            rho=obj.cell{i}.rho;
            obj.cell{i}.T = T;
            obj.cell{i}.TvarT = T .^ 2;
            obj.cell{i}.TvarE = rho .* T .^ 2;
            obj.cell{i}.TvarI = obj.cell{i}.TvarT - obj.cell{i}.TvarE;

            % INDIVIDUAL PASSES
            obj.cell{i}.Tind=zeros(numel(P.inds),1);
            for j = 1:numel(P.inds)
                p=P.inds{j}.curs{i};

                obj.cell{i}.Tind(j)       = p.tMU;
                obj.cell{i}.TindCI(j,:)   = p.tCI;
                obj.cell{i}.TindSE(j,:)   = p.tSE;
            end

            %Ppasses
            %p=P.inds{j}.curs{i};
        end
    end
    function obj=run_ind(obj,bBoot)
        if ~exist('bBoot','var') || isempty(bBoot)
            bBoot=0;
        end
        obj.init_ind();
        IND=obj.get_ind_ind();
        PASSES=obj.get_passes();
        p=Pr(size(IND,1),'Fitting all agreement',[],1);
        for i = 1:size(IND,1)
            p.u();
            std=IND(i,1);
            subj=IND(i,2);
            ind=IND(i,3);
            passes=PASSES(ind,:);
            if bBoot
                obj.cell{i}.fit_boot();
            else
                obj.cell{i}.run();
            end
        end
        p.c();
    end
    function obj=select_ind(obj,std,subj,passes)
        obj.stdX=var_select(obj.stdXall,std,subj,passes,1);
        obj.cmpX=var_select(obj.cmpXall,std,subj,passes,1);
        obj.RcmpChs=var_select(obj.RcmpChsAll,std,subj,passes,0);

        function out=var_select(var,std,subj,passes,bSmash)
            a=cell(numel(passes),1);
            for i = 1:numel(passes)
                val=squeeze(var(:,std,subj,passes(i)));
                val=val(:);
                a{i}=val;
                if i == 1
                    beq=1;
                else
                    beq=isequaln(last,val) & beq;
                end
                last=val;
            end
            if beq && bSmash
                out=val;
                return
            elseif bSmash
                error('data does not lineup')
            end
            out=cat(2,a{:});
        end
    end
%% COMB SUBJ
    function obj=init_subj_comb(obj)
        PASSES=obj.get_passes();
        IND=obj.get_ind_subj();
        obj.cellSubj=cell(size(IND,1),1);
        Opts=obj.get_Opts('cellSubj');
        for i = 1:size(IND,1)
            std=IND(i,1);
            ind=IND(i,2);
            passes=PASSES(ind,:);

            obj.select_subj_comb(std,passes);

            obj.cellSubj{i}=DVcorr(obj.stdX,obj.cmpX,obj.RcmpChs,Opts);
        end
    end
    function obj=run_subj_comb(obj)
        obj.init_subj_comb();
        PASSES=obj.get_passes();
        IND=obj.get_ind_subj();
        p=Pr(size(IND,1),'Fitting subj comb agreement',[],1);
        for i = 1:size(IND,1)
            p.u();
            std=IND(i,1);
            ind=IND(i,2);
            passes=PASSES(ind,:);


            obj.cellSubj{i}.run();
        end
        p.c();
    end
    function obj=select_subj_comb(obj,std,passes)
        obj.stdX=var_select(obj.stdXall,std,passes,1);
        obj.RcmpChs=var_select(obj.RcmpChsAll,std,passes,0);

        function out=var_select(var,std,passes,bSmash)
            a=cell(numel(passes),1);
            for i = 1:numel(passes)
                val=squeeze(var(:,std,:,passes(i)));
                val=val(:);
                a{i}=val;
                if i == 1
                    beq=1;
                else
                    beq=isequal(last,val) & beq;
                end
                last=val;
            end
            if beq && bSmash
                out=val;
                return
            elseif bSmash
                error('data does not lineup')
            end
            out=cat(2,a{:});
       end
    end
%% MAIN
end
end
