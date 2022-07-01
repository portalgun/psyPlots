classdef DVcorrs_plots < handle
methods
    function plot_all(obj,bSave)
        if exist('bSave','var') && ~isempty(bSave)
            obj.bSaveFig=bSave;
        end
        %obj.plot_counts_ind_all();
        obj.plot_ellipse_ind_all();
        obj.plot_magr_ind_all();

        %% IND
        obj.plot_bin_rho_ind_all();
        obj.plot_bin_ratio_ind_all();
        obj.plot_scatter_rho_ind_all();
        obj.plot_scatter_ratio_ind_all();

        %% ALL
        obj.plot_bin_rho_comb_all();
        obj.plot_bin_ratio_comb_all();
        obj.plot_scatter_rho_comb_all();
        obj.plot_scatter_ratio_comb_all();

        % STD
        obj.plot_bin_rho_comb_std_all();
        obj.plot_bin_ratio_comb_std_all();

        % CMP
        obj.plot_scatter_rho_comb_cmp();
        obj.plot_scatter_ratio_comb_cmp();

        % SUBJ
        obj.plot_bin_rho_comb_subj();
        obj.plot_bin_ratio_comb_subj();
        obj.plot_scatter_rho_comb_subj();
        obj.plot_scatter_ratio_comb_subj();
    end

%% IND ALL
    function plot_counts_ind_all(obj)
        for i = 1:obj.nSubj
            fig=Fig.new();
            obj.plot_counts_p(i);
            if obj.bSaveFig
                obj.save_fig(fig,i,'counts','inds');
            end
        end
    end
    function plot_ellipse_ind_all(obj)
        for i = 1:obj.nSubj
            fig=Fig.new();
            obj.plot_ellipse_p(i);
            if obj.bSaveFig
                obj.save_fig(fig,i,'DVellipse','inds');
            end
        end
    end
    function plot_magr_ind_all(obj,Opts)
        % HERE
        for i = 1:obj.nSubj
            fig=Fig.new();
            obj.plot_magr_p(Opts,i);
            if obj.bSaveFig
                obj.save_fig(fig,i,'magr','inds');
            end
        end
    end
    function plot_scatter_abs_thresh_ind_all(obj,Opts,varargin)
        % HERE
        if ~exist('Opts','var'); Opts=struct(); end
        fig=Fig.new();
        obj.plot_scatter_abs_p(0,0,0,0,Opts,varargin{:});
        if obj.bSaveFig
            obj.save_fig(fig,[],'scatterAbs','inds');
        end
    end
    function plot_scatter_abs_ind_all(obj,Opts,varargin)
        % HERE
        if ~exist('Opts','var'); Opts=struct(); end
        fig=Fig.new();
        obj.plot_scatter_abs_p(1,0,0,0,Opts,varargin{:});
        if obj.bSaveFig
            obj.save_fig(fig,[],'scatterAbs','inds');
        end
    end
    function plot_scatter_abs_split_ind_all(obj,Opts,varargin)
        % HERE
        if ~exist('Opts','var'); Opts=struct(); end
        fig=Fig.new();
        obj.plot_scatter_abs_p(1,0,0,1,Opts,varargin{:});
        if obj.bSaveFig
            obj.save_fig(fig,[],'scatterAbsSplit','inds');
        end
    end
%% STD AL[],L
    function plot_bin_ratio_comb_std_all(obj,varargin)
        for i = 1:obj.nSubj
            fig=Fig.new();
            obj.plot_bin_rho_comb_std_ind(i,varargin{:})
            if obj.bSaveFig
                obj.save_fig(fig,[],'binRatio','std');
            end
        end
    end
%% IND
    function plot_counts_ind(obj,ind)
        Fig.new()
        obj.plot_counts_p(ind);
        if obj.bSaveFig
            obj.save_fig(fig,ind,'counts','inds');
        end
    end
    function plot_ellipse_ind(obj,Opts,subj)
        Fig.new()
        obj.plot_ellipse_p(1,0,0,Opts,subj);
        if obj.bSaveFig
            obj.save_fig(fig,subj,'DVellipse','inds');
        end
    end
    function plot_magr_ind(obj,Opts,ind)
        Fig.new()
        obj.plot_magr_p(Opts,ind);
        if obj.bSaveFig
            obj.save_fig(fig,ind,'magr','inds');
        end
    end
    % no scatter
% COMB STD IND
    function plot_bin_ratio_comb_std_ind(obj,subj,varargin)
        obj.plot_bin_rho_ratio_p('ratio','std',subj,varargin{:});
    end
%%%%%%%%%%%%%%%
%%  COMB ALL
    %
    function plot_scatter_rho_comb_all(obj,Opts)
    % HERE
        fig=Fig.new();
        if ~exist('Opts','var'); Opts=struct(); end
        obj.plot_scatter_rho_p(1,1,Opts);
        if obj.bSaveFig
            obj.save_fig(fig,[],'scatterRho','all');
        end
    end
%%%
    function plot_scatter_abs_thresh_comb_all(obj,varargin)
        fig=Fig.new();
        obj.plot_scatter_abs_p(0,1,1,0,varargin{:});
        if obj.bSaveFig
            obj.save_fig(fig,[],'scatterAbs','inds');
        end
    end
    function plot_scatter_abs_thresh_comb_subj(obj,varargin)
        fig=Fig.new();
        obj.plot_scatter_abs_p(0,0,1,1,varargin{:});
        if obj.bSaveFig
            obj.save_fig(fig,[],'scatterAbs','inds');
        end
    end
    % no scatter
% COMB CMP
    % no bin
    function plot_scatter_rho_comb_cmp(obj)
        fig=Fig.new();
        obj.plot_scatter_rho_ratio_p('rho',0,1);
        if obj.bSaveFig
            obj.save_fig(fig,[],'scatterRho','cmp');
        end
    end
% COMB subj new
    function plot_scatter_rho_comb_subj(obj,Opts)
        if ~exist('Opts','var'); Opts=struct(); end
        fig=Fig.new();
        obj.plot_scatter_rho_p(0,1,0,Opts);
        if obj.bSaveFig
            obj.save_fig(fig,[],'scatterRatio','subj');
        end
    end

    function plot_bin_rho_comb_std_all(obj,varargin)
        for i = 1:obj.nSubj
            fig=Fig.new();
            obj.plot_bin_rho_comb_std_ind(i,varargin{:})
            if obj.bSaveFig
                obj.save_fig(fig,[],'binRho','std');
            end
        end
    end
%% BIN
    function plot_bin_rho_same_all(obj,Opts,varargin)
        fig=Fig.new();
        obj.plot_bin_rho_p(0,0,1,0,Opts,varargin{:});
        if obj.bSaveFig
            obj.save_fig(fig,[],'binRho','all');
        end
    end
    function plot_bin_rho_comb_all(obj,Opts,varargin)
        fig=Fig.new();
        obj.plot_bin_rho_p(0,0,1,1,Opts,varargin{:});
        if obj.bSaveFig
            obj.save_fig(fig,[],'binRho','all');
        end
    end
    function plot_bin_ratio_comb_all(obj,Opts,varargin)
        fig=Fig.new();
        obj.plot_bin_rho_p(1,0,1,1,Opts,varargin{:});
        if obj.bSaveFig
            obj.save_fig(fig,[],'binRho','all');
        end
    end
    function plot_bin_ratio_same_all(obj,Opts,varargin)
        fig=Fig.new();
        obj.plot_bin_rho_p(1,0,1,0,Opts,varargin{:});
        if obj.bSaveFig
            obj.save_fig(fig,[],'binRho','all');
        end
    end



    function plot_scatter_rho_ind_comb(obj,Opts)
    % HERE
        if ~exist('Opts','var'); Opts=struct(); end
        fig=Fig.new();
        obj.plot_scatter_ratio_p(0,1,0,Opts);
        if obj.bSaveFig
            obj.save_fig(fig,[],'scatterRatio','subj');
        end
    end

    function plot_scatter_ratio_ind_all(obj,Opts)
    % HERE
        if ~exist('Opts','var'); Opts=struct(); end
        fig=Fig.new();
        obj.plot_scatter_ratio_p(1,0,0,Opts);
        if obj.bSaveFig
            obj.save_fig(fig,[],'scatterRatio','inds');
        end
    end
    function plot_scatter_ratio_ind_comb(obj,Opts)
    % HERE
        if ~exist('Opts','var'); Opts=struct(); end
        fig=Fig.new();
        obj.plot_scatter_ratio_p(1,1,0,Opts);
        if obj.bSaveFig
            obj.save_fig(fig,[],'scatterRatio','subj');
        end
    end
    function plot_scatter_ratio_comb_subj(obj,Opts)
    % HERE
        if ~exist('Opts','var'); Opts=struct(); end
        fig=Fig.new();
        obj.plot_scatter_ratio_p(1,0,1,Opts);
        if obj.bSaveFig
            obj.save_fig(fig,[],'scatterRatio','subj');
        end
    end
    function plot_scatter_ratio_comb_all(obj,Opts)
    % HERE
        fig=Fig.new();
        if ~exist('Opts','var'); Opts=struct(); end
        obj.plot_scatter_ratio_p(1,1,1,Opts);
        if obj.bSaveFig
            obj.save_fig(fig,[],'scatterRatio','all');
        end
    end
end
end
