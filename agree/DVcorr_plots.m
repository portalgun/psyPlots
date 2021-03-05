classdef DVcorr_plots < handle
%%%%%%%%
%% PLOT
    % response
    % ellipse
    % magr
    % rho bin
    % rho scatter
    % ratio bin
    % ratio scatter

methods
%% RESPONSE
    function plot_response_count_all(obj)

        cmp=round(unique(obj.cmpX),2);
        for i = 1:obj.nCmp
            c=cmp(i);
            c=num2str(c,'%.2f');
            ctitl{i}=['Cmp ' c ];
        end


        Opts=struct();
        Opts.xticks=[1 2 3 4];
        Opts.xtickLabels={'--','-+','+-','++'};
        Opts.yticks=[.0 .2 .4 .6 .8];
        Opt.ylimSpace=.15;

        C=obj.nCmp;
        sp=subPlots([1,C],'Joint Response Type','Proportion Chosen','','',ctitl,Opts);
        for c = 1:C
            sp.select(1,c);
            obj.select(c);
            obj.plot_response_count_p();
        end
        sp.finalize();
    end
    function plot_response_count(obj,ind)
        if exist('ind','var') && isempty(ind)
            obj.select(ind);
        end
        figure(nFn)
        obj.plot_response_count_p();
        obj.label_response_count_p();
    end
%% ELIPSE
    function plot_ellipse_all_same(obj)
        obj.plot_ellipse_all_same_p();
        %obj.format_ellipse_all_same_p();
    end
    function plot_ellipse_all_same_p(obj)
        C=obj.nCmp;
        obj.get_colors();

        flag=0;
        if numel(unique(obj.RHO))==1
            flag=1;
            txt=['\rho =' num2str(obj.RHO(1),'%0.2f')];
        end


        for c = 1:C
            obj.get_color(c);
            obj.select(c);
            obj.plot_ellipse_p(); hold on
            %obj.label_ellipse_p();
            if flag
                text(1,-3,txt,'fontSize',16);
            end
        end
        hold off
    end
    function function_ellipse_all_same(obj)
        str=sprintf('%.2f',unique(round(obj.std,2)));
        title(['Std ' str]);
    end

    function plot_ellipse_all(obj)
        C=obj.nCmp;
        for c = 1:C
            obj.select(c);
            subPlot([1,C],1,c);
            obj.plot_ellipse_p();
            obj.label_ellipse_p();
            if c~=1
                ylabel('');
                yticklabels('');
            end
            obj.label_cmp();
        end
    end
    function plot_ellipse(obj,ind)
        if exist('ind','var') && isempty(ind)
            obj.select(ind);
        end
        figure(nFn)
        obj.plot_ellipse_p();
        obj.label_ellipse_p();
    end
%% MAGR
    function plot_magr(obj)
        obj.plot_magr_p;
        obj.format_magr;
    end
%% RHO
    function plot_rho_scatter(obj)
        obj.plot_rho_scatter_p();
        obj.format_rho_scatter_p();
    end
    function plot_rho_bin(obj)
        obj.plot_rho_bin_p();
        obj.format_rho_bin_format_p();
    end
%% RATIO
    function plot_ratio_scatter(obj)
        obj.plot_ratio_scatter_p();
        obj.format_ratio_scatter_format_p();
    end
    function plot_ratio_bin(obj)
        obj.plot_ratio_bin_p();
        obj.format_ratio_bin_format_p();
    end
%% ABS
    function plot_abs_scatter(obj)
        obj.plot_abs_scatter_p();
        obj.format_abs_scatter_format_p();
    end
end
end
