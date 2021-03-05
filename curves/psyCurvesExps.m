classdef psyCurvesExps < handle
properties
    dt
    EXPS
    Opts

    inds
    comb
    passComb %combined subjects and pass
    passInds % combined pass

    bSaveFig
    figSaveTypes={'.svg','.fig','.png'};
end
methods
    function obj=psyCurvesExps(dt,EXPS,Opts)
        obj.dt=dt;
        obj.EXPS=EXPS;
        obj.Opts=Opts;
    end
    function obj=rep_inds(obj)
        for i = 1:length(obj.dt)
            obj.dt{i}=structRepInd(obj.dt{i},1);
        end
    end
%% fit
    function obj=fit_inds(obj)
    % COMBINE NOTHING
        N=length(obj.dt);
        obj.inds=cell(N,1);
        for i = 1:N
            EXP=obj.EXPS.get_EXP(i);
            obj.inds{i}=psyCurves(obj.dt{i},EXP,[],obj.Opts);
            obj.inds{i}.fit_all();
        end
    end
    function obj=fit_pass_inds(obj)
    % COMBINE PASSES
        ind=obj.EXPS.get_unique_exp_ind();
        n=max(ind);
        EXPS=obj.EXPS.combine_passes();
        obj.passInds=cell(n,1);
        for i = 1:n
            d=structMerge(obj.dt{ind==i});
            d=structShrink(d,1,4);
            EXP=EXPS.get_EXP(i);

            obj.passInds{i}=psyCurves(d,EXP,obj.Opts);
            obj.passInds{i}.fit_all();
        end
    end

    function obj=fit_pass_comb(obj)
    % COMBINE SUBJECTS AND PASSES
        ind=obj.EXPS.get_unique_exp_ind();
        n=max(ind);
        EXPS=obj.EXPS.combine_passes();
        obj.passComb=cell(n,1);
        for i = 1:n
            d=structMerge(obj.dt{ind==i});
            d=structShrink(d,1,4);
            d=structShrink(d,1,3);
            EXP=EXPS.get_EXP(i);

            obj.passComb{i}=psyCurves(d,EXP,obj.Opts);
            obj.passComb{i}.fit_all();
            obj.passComb{i}
        end
    end
    function obj=fit_comb(obj)
    % COMBINE SUBJECTS
        N=length(obj.dt);
        obj.comb=cell(N,1);
        for i = 1:N
            d=structShrink(obj.dt{i},1,3);
            EXP=obj.EXPS.get_EXP(i);
            obj.comb{i}=psyCurves(d,EXP,obj.Opts);
            obj.comb{i}.fit_all();
        end
    end

%% opts
    function obj=apply_EXP_all(obj)
        obj.apply_EXP('inds');
        obj.apply_EXP('comb');
        obj.apply_EXP('passComb');
        obj.apply_EXP('passInds');
    end
    function obj=apply_EXP(obj,fld)
        if startsWith(fld,'pass')
            EXPS=obj.EXPS.combine_passes();
        else
            EXPS=obj.EXPS;
        end
        N=size(obj.(fld),1);
        for i = 1:N
            EXP=EXPS.get_EXP(i);
            obj.(fld){i}.EXPS=EXP;
        end
    end
    function obj=parse_opts_inds(obj)
        N=size(obj.inds,1);
        for i = 1:N
            obj.inds{i}.parse_opts();
        end
    end
%% THRESH
    function obj=plot_thresh(obj,fld,bAverage,bSame)
        if ~exist('bAverage','var') || isempty(bAverage)
            bAverage=0;
        end
        if ~exist('bSame','var') || isempty(bSame)
            bSame=0;
        end
        if bSame
            str='Same';
        else
            str='';
        end
        N=size(obj.(fld),1);
        for i = 1:N
            fig=figure(nFn);
            obj.(fld){i}.plot_thresh_all_p(bAverage,bSame);
            if obj.bSaveFig
                obj.save_fig(fig,i,['thesh' str],fld)
            end
        end
    end
%% PSY
    function obj=plot_curves(obj,fld)
        N=size(obj.(fld),1);
        for i=1:N
            fig=figure(nFn);
            obj.(fld){i}.plot_curve_all_p();
            if obj.bSaveFig
                obj.save_fig(fig,i,'psyCurve',fld);
            end
        end
    end
%% SAVE
    function obj=save_fig(obj,fig,ind,figType,dataType)
        name=obj.gen_fname(ind,figType,dataType);
        dir=BLdirs('fig');
        for i = 1:length(obj.figSaveTypes)
            ext=obj.figSaveTypes{i};
            fname=[dir name ext];
            saveas(fig,fname);
        end
        close(fig);
    end
    function name=gen_fname(obj,ind,figType,dataType)
        if startsWith(dataType,'pass')
            EXPS=obj.EXPS.combine_passes();
        else
            EXPS=obj.EXPS;
        end
        EXP=EXPS.get_EXP(ind);
        name=[figType '-' dataType und EXP.name];
    end
    function obj=save_off(obj)
        obj.bSaveFig=0;
    end
    function obj=save_on(obj)
        obj.bSaveFig=1;
    end
end
end
