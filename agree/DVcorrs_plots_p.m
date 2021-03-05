classdef DVcorrs_plots_p < handle & DVcorr_plots_p_opts
methods(Hidden=true)

    function plot_scatter_rho_p(obj,bSame,bComb,Opts,varargin)

        Opts=obj.get_rho_scatter_opts(Opts);
        SP=sp_tools(obj,0,bSame,bComb,Opts,varargin);

        % RHO
        R=obj.get_RHO(0,SP);

        for i = 1:SP.nInd
            SP.get_iter(i,R);

            if SP.continueflag; continue; end

            hold on
            SP.plot_interv();
            if SP.MODE==1;
                obj.cell{i}.plot_ratio_scatter_p(sym,stagger);
            else
                SP.plot();
            end
        end

        SP.finalize();
    end
    function plot_scatter_ratio_p(obj,bRatio,bSame,bComb,Opts,varargin)
        if bRatio
            Opts=obj.get_ratio_scatter_opts(Opts);
        else
            Opts=obj.get_rho_scatter_opts(Opts);
        end
        SP=sp_tools(obj,0,bSame,bComb,Opts,varargin);

        % RHO
        R=obj.get_RHO(bRatio,SP);
        R.RHO

        for i = 1:SP.nInd
            SP.get_iter(i,R);

            if SP.continueflag; continue; end

            hold on
            if isfield(Opts,'bStagger') && Opts.bStagger & bRatio
                obj.cell{i}.plot_ratio_scatter_p(SP.iter.sym,SP.iter.stagger);
            elseif isfield(Opts,'bStagger') && Opts.bStagger
                obj.cell{i}.plot_rho_scatter2_p(SP.iter.sym,SP.iter.stagger);
            else
                %SP.plot_interv();
                plot_fun(SP);

            end

        end

        SP.finalize();
        function plot_fun(SP)
            I=SP.iter;
            h2=interv(I.X,I.U,I.L,I.color); hold on;
            plot(I.X,I.m,[I.sym 'k'],'MarkerFaceColor','w','MarkerSize',10,'LineWidth',2); hold on
        end

    end
    function plot_magr_p(obj,Opts,varargin)
        [subj,sub,passes]=obj.parse_varargin(varargin);
        Opts=obj.get_magr_opts(Opts,sub);
        SP=sp_tools(obj,1,0,0,Opts,subj,passes);

        for i = 1:SP.nInd
            SP.get_iter(i);

            if SP.continueflag==1; continue; end

            hold on
            obj.cell{i}.plot_magr_p();
        end
        SP.finalize();
    end
    function plot_counts_p(obj,subj,Opts,varargin)
        SP=sp_tools(obj.EXPS,0,0,0,Opts,subj)
        Opts=obj.get_counts_opts(Opts);
        SP=sp_tools(obj.EXPS,1,0,0,Opts,varargin)


        for i = 1:SP.nInd
            SP.get_iter(i);

            if SP.continueflag==1; continue; end

            obj.cell{ind}.select(cmp);
            obj.cell{ind}.plot_response_count_p();
            hold on
        end
        SP.finalize();


    end
    function plot_ellipse_p(obj,bInd,bSame,bComb,Opts,varargin)
        Opts=obj.get_ellipse_opts(Opts,varargin{1});
        varargin{1}
        SP=sp_tools(obj,bInd,bSame,bComb,Opts,varargin{:});
        for i = 1:SP.nInd
            SP.get_iter(i);
            SP.continueflag

            if SP.continueflag==1; continue; end

            obj.cell{i}.plot_ellipse_all_same_p();
            hold on
        end
        SP.finalize();
    end
    function plot_scatter_abs_p(obj,bInd,bSame,bComb,bSplit,Opts,varargin)
        SP=sp_tools(obj,bInd,bSame,bComb,Opts,varargin);
        R=obj.get_RHO(1,SP);

        bComb
        cm=cmap('hsv',obj.nStd);
        for i = 1:SP.nInd
            SP.get_iter(i,R);

            if SP.continueflag==1; continue; end

            hold off

            if SP.MODE==3
                plot_fun(SP);
            elseif SP.MODE==2
                plot_fun(SP);
            end
        end
        SP.finalize();
        sgtitle('Stimulus variability and internal noise estimates');
        function plot_fun(SP)
            I=SP.iter;
            h1=interv(I.X,I.UI,I.LI,'b'); hold on;
            h2=interv(I.X,I.UE,I.LE,'r'); hold on;
            plot(I.X,I.mI,[I.sym 'k'],'MarkerFaceColor','w','MarkerSize',10,'LineWidth',2); hold on
            plot(I.X,I.mE,[I.sym 'k'],'MarkerFaceColor','w','MarkerSize',10,'LineWidth',2);
            legend([h1,h2],{'\sigma_{Int}^2','\sigma_{Ext}^2'});

            %plot(I.X, I.mI./I.mE, [I.sym 'k'],'MarkerFaceColor','w','MarkerSize',10,'LineWidth',2);

            % SANITY CHECK
            % should match rel
            %plot(I.X, I.mI./I.mE, [I.sym 'k'],'MarkerFaceColor','w','MarkerSize',10,'LineWidth',2);

            % SANITY CHECK
            % SHOULD MATCH THRESHOLDS
            %plot(I.X, sqrt(I.mI+I.mE), [I.sym 'k'],'MarkerFaceColor','w','MarkerSize',10,'LineWidth',2);
        end

    end
    function plot_bin_rho_p(obj,bRatio,bInd,bSame,bComb,Opts,varargin)
        % I=ismember(IND(:,2),subj);
        % ind  - 0 0 0, R=RHO(I,:);
        % all  - 0 1 1, R=transpose(RHO(:));
        % std  - 0 1 0, R=RHO(I,:); R=R(:);
        % subj - 0 0 1  r=RHO(ind,:); R(i,:)=r(:)';

        [subj,sub,passes]=obj.parse_varargin(varargin);
        if bRatio
            Opts=obj.get_bin_ratio_opts(Opts,sub);
        else
            Opts=obj.get_bin_rho_opts(Opts,sub);
        end

        SP=sp_tools(obj,bInd,bSame,bComb,Opts,subj,passes);
        R=obj.get_RHO(bRatio,SP);

        BE=[];
        if isfield(Opts,'nBin') && ~isempty(Opts.nBin)
            BE=Opts.nBin;
        end

        bEdge=0;
        if isfield(Opts,'edges') && ~isempty(Opts.edges)
            BE=Opts.edges;
            bEdge=1;
        end


        if SP.MODE==3
            R=R.RHO(:,SP.I);
            R=transpose(R(:));
        else
            R=transpose(R.RHO(:));
        end

        %histogram(R,8)
        hst=histo(R,BE,'bLog',Opts.bLog);
        if ~bEdge
            hst=histo(R,hst.edges,'bLog',Opts.bLog,'bSameColor',1);
        end
        hst.edges
        for i = 1:size(R,1)
            i
            SP.get_iter(i);

            if SP.continueflag==1; continue; end

            hold on
            hst.plot_ind(i);
        end
        hst.c();
        SP.finalize();
    end
    function plot_bin_ratio_p(obj,bInd,bSame,bComb,Opts,varargin)
        % I=ismember(IND(:,2),subj);
        % ind  - 0 0 0, R=RHO(I,:);
        % all  - 0 1 1, R=transpose(RHO(:));
        % std  - 0 1 0, R=RHO(I,:); R=R(:);
        % subj - 0 0 1  r=RHO(ind,:); R(i,:)=r(:)';

        [subj,sub,passes]=obj.parse_varargin(varargin);
        Opts=obj.get_bin_rho_opts(Opts,sub);
        SP=sp_tools(obj,bInd,bSame,bComb,Opts,subj,passes);
        R=obj.get_RHO(1,SP);


        if SP.MODE==3
            R=R.RHO(:,SP.I);
            R=transpose(R(:));
        end

        hst=histo(R,[],'bLog',Opts.bLog);
        hst=histo(R,hst.edges,'bLog',Opts.bLog,'bSameColor',1);
        for i = 1:size(R,1)
            SP.get_iter(i);

            if SP.continueflag==1; continue; end

            hold on
            hst.plot_ind(i);
        end
        hst.c();
        SP.finalize();
    end
    function plot_thresh_p(obj,subj)

        [subj,sub,passes]=obj.parse_varargin(varargin);
        Opts=obj.get_bin_rho_opts(Opts,sub);
        SP=sp_tools(obj,bInd,bSame,bComb,Opts,subj,passes);
        R=obj.get_RHO(1,SP);

        T=zeros(obj.nSubj,obj.nStd, obj.nPss);
        TCI=zeros(obj.nSubj,obj.nStd, obj.nPss,2);

        hold off

        Xt=TCI(:,:,1,:);
        Yt=TCI(:,:,2,:);
        X=T(:,:,1);
        Y=T(:,:,2);
        for i = 1:max(IND(:,2))
            hold on

            x=squeeze(X(i,:));
            y=squeeze(Y(i,:));
            dx=squeeze(Xt(i,:,:));
            dy=squeeze(Yt(i,:,:));

            errorbarSane(x,y,dx,dy,'k-','LineWidth',obj.LineWidth)
            %p(i)=plot(T(i,:,1),T(i,:,2),[shapes{i} 'k'],'LineWidth',obj.LineWidth,'MarkerFaceColor','w','MarkerSize',10)
            text{i}=['Observer ' num2str(i)];
        end
        for i = 1:max(IND(:,2))
            x=squeeze(X(i,:));
            y=squeeze(Y(i,:));
            p(i)=plot(x,y,[shapes{i} 'k'],'LineWidth',obj.LineWidth,'MarkerFaceColor','w','MarkerSize',10);
        end
        legend(p,text,'location','northwest')

        %sp.finalize();

    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    function get_thresh(obj,SP)
        T=zeros(obj.nSubj,obj.nStd, obj.nPss);
        TCI=zeros(obj.nSubj,obj.nStd, obj.nPss,2);
        for i = 1:SP.nInd
            T(subjInd,std,passes(1))=obj.cell{i}.Tind(passes(1));
            T(subjInd,std,passes(2))=obj.cell{i}.Tind(passes(2));

            TCI(subjInd,std,passes(1),:)=obj.cell{i}.TindCI(passes(1));
            TCI(subjInd,std,passes(2),:)=obj.cell{i}.TindCI(passes(2));

            TSE(subjInd,std,passes(1))=obj.cell{i}.TindSE(passes(1));
            TSE(subjInd,std,passes(2))=obj.cell{i}.TindSE(passes(2));

            %obj.cell{i}.Tind
            % XXX HERE
        end
    end
    function R=get_RHO(obj,bRatio,SP)
        R=struct;
        tmp=num2cell(SP.passesSel);
        [VARI,VARE,~,STDS]=obj.get_var_src_all('cell',tmp{:});
        RHO=obj.get_rho_all('cell',tmp{:});
        boot=[];
        try
            boot=obj.get_boot_all('cell',tmp{:});
        end

        bRmNeg=0;
        if bRatio
            RHO=sqrt((1-RHO)./(RHO));
            bRmNeg=1;
        end


        if SP.MODE==3
            [m,U,L]=SP.average_stat(transpose(RHO),1,bRmNeg);

            if ~isempty(VARI)
                [mI,UI,LI]=SP.average_stat(VARI,1,bRmNeg);
                [mE,UE,LE]=SP.average_stat(VARE,1,bRmNeg);
            end
        elseif SP.MODE==2
            %RHO
            RHO=transpose(mean(RHO,1));
            [m,U,L]=SP.average_stat(transpose(RHO),0,bRmNeg);
            if ~isempty(VARI)

                [mI,UI,LI]=SP.average_stat(transpose(VARI),0,bRmNeg);
                [mE,UE,LE]=SP.average_stat(transpose(VARE),0,bRmNeg);
            end
       % elseif ~isempty(boot)
       %     L=boot.RHO_CI(:,:,1);
       %     U=boot.RHO_CI(:,:,2);
       %     if bRatio
       %         L=sqrt((1-L)./(L));
       %         U=sqrt((1-U)./(U));
       %     end
       %     m=RHO;
        else
            R.RHO=RHO;
            return
        end
        R.RHO=RHO;
        R.m=m;
        R.U=U;
        R.L=L;
        R.RHO

        if ~isempty(VARI) && exist('mI','var')
            R.mI=mI;
            R.UI=UI;
            R.LI=LI;

            R.mE=mE;
            R.UE=UE;
            R.LE=LE;
        end
    end

%% FIG HANDLING
    function obj=save_off(obj)
        obj.bSaveFig=0;
    end
    function obj=save_on(obj)
        obj.bSaveFig=1;
    end
    function obj=save_fig(obj,fig,ind,figType,dataType)
        name=obj.gen_fname(ind,figType,dataType);
        dir=BLdirs('fig');
        for i = 1:length(obj.figSaveTypes)
            ext=obj.figSaveTypes{i};
            fname=[dir name ext];
            fname
            saveas(fig,fname);
        end
        close(fig);
    end
    function name=gen_fname(obj,ind,figType,dataType)
        if ~isempty(ind)
            subj=obj.EXPS.subjs{ind};
            name=[figType '-' dataType '-' subj und obj.EXPS.name];
        else
            name=[figType '-' dataType und obj.EXPS.name];
        end
    end
    function [subj,sub,passes]=parse_varargin(obj,vin)
        if isempty(vin)
            subj=[];
            sub=[];
            passes=[];
            return
        end

        if length(vin) > 0
            subj=vin{1};
            sub=obj.EXPS.subjs(subj);
        else
            subj=[];
            sub=[];
        end
        if length(vin) > 1
            passes=vin{2};
        else
            passes=[1,2];
        end

    end
end
end
