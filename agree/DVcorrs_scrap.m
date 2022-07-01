%% COMB
    function plot_ellipse_comb_subj(obj)
        fig=Fig.new();
        obj.plot_ellipse_subj_p();
        if obj.bSaveFig
            obj.save_fig(fig,[],'DVellipse','comb');
        end
    end
    function plot_magr_comb_subj(obj)
        fig=Fig.new();
        obj.plot_magr_subj_p();
        if obj.bSaveFig
            obj.save_fig(fig,[],'magr','comb');
        end
    end
% private
    function plot_ellipse_subj_p(obj)
        PASSES=obj.get_passes();
        IND=obj.get_ind_subj();
        sz=[numel(unique(IND(:,1))) numel(unique(IND(:,2)))];

        if ~nflds(obj.EXPS) > 0
            ctitl=obj.EXPS.get_rc_title('imgDim','pass');
            rtitl=obj.EXPS.get_rc_title('std');
            ctitl=join(nchoosek(ctitl,obj.nPssAtOnce));
        else
            ctitl=strcat('Passes ', strsplit(Num.toStr(PASSES),';'));
            rtitl=strcat('Std ',   strsplit(Num.toStr([1:obj.nStd]),','));
        end
        if obj.bFlipRC
            r=rtitl;
            c=ctitl;
            rtitl=c;
            ctitl=r;
            sz=fliplr(sz);
        end
        sp=SubPlots(sz,'1st pass decision variable','2nd pass decision variable','DV correlation - combined subject data',rtitl,ctitl);
        for i = 1:size(IND,1)
            ind=IND(i,end);
            std=IND(i,1);
            ind=IND(i,end);
            passes=PASSES(ind,:);

            if obj.bFlipRC
                sp.select(ind,std);
            else
                sp.select(std,ind);
            end
            hold on
            obj.cellSubj{i}.plot_ellipse_all_same_p();
        end
        sp.finalize();
    end
    function plot_scatter_subj_p(obj)
        PASSES=obj.get_passes();
        IND=obj.get_ind_subj();
        sz=[1 obj.nPss];

        ylabl=obj.cell{1}.ylabel_rho_scatter_p();
        xlabl=obj.cell{1}.xlabel_rho_scatter_p();
        x=unique(obj.stdXall);
        y=[.01 .03 .1 .3  1 3 10];

        ctitl=strcat('Passes ', strsplit(Num.toStr(PASSES),';'));
        rtitl=strcat('Subj ',   strsplit(Num.toStr([1:obj.nSubj]),','));
        sp=SubPlots(sz,xlabl,ylabl,'DV correlation',rtitl,ctitl,[],x,y);
        for i = 1:size(IND,1)
            ind=IND(i,end);
            passes=PASSES(ind,:);

            %subPlot(sz,subj,ind,1,x);
            sp.select(1,ind);
            hold on
            obj.cellSubj{i}.plot_rho_scatter_p();

            obj.cellSubj{i}.format_rho_scatter_p();
            rang=max(x)-min(x);
            xlim([min(x)-rang*.15 max(x)+rang*.15]);
            %Fig.format();

        end
        sp.finalize();
    end
    function plot_magr_subj_p(obj)
        PASSES=obj.get_passes();
        IND=obj.get_ind_subj();
        sz=[numel(unique(IND(:,1))) obj.nPss/obj.nPssAtOnce];

        ylabl='% Comparison chosen';
        xlabl='% Agreement';

        if nflds(obj.EXPS) > 0
            ctitl=obj.EXPS.get_rc_title('imgDim','pass');
            rtitl=obj.EXPS.get_rc_title('std');
            ctitl=join(nchoosek(ctitl,obj.nPssAtOnce));
        else
            ctitl=strcat('Passes ', strsplit(Num.toStr(PASSES),';'));
            rtitl=strcat('Std ',   strsplit(Num.toStr([1:obj.nStd]),','));
        end
        sp=SubPlots(sz,xlabl,ylabl,'Pass Agreement',rtitl,ctitl,[]);
        sp.ylim=[.4 1];
        for i = 1:size(IND,1)
            std=IND(i,1);
            ind=IND(i,end);
            passes=PASSES(ind,:);


            sp.select(std,ind);
            hold on
            obj.cellSubj{i}.plot_magr_p();
        end
        sp.finalize();
    end
