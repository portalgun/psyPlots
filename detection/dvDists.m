classdef dvDists < handle
properties
    mDists

    A
    B
    xB
    xA
    fA
    fB
    fAA
    fAB
    fBA
    fBB

    LLRa
    LLRb
    LLRc

    LLRaCtrs
    LLRbCtrs
    LLRcCtrs

    ctrs
    x
    y
    fPrs
    tPrs
    AUCs

    IN
    INname
    OUT

    method='spline';
    names
end
properties(Hidden = true)
    n
    fignum=999
end
methods
    function obj=dvDists(varargin)
        obj.method='spline';
        obj.fignum=999;
        obj.mDists=varargin;
        obj.n=length(obj.mDists);
        if obj.n==0
            return
        end

        obj.get_names();
        obj.get_ctrs();
        obj.get_y();

        obj.get_conditionals();
        obj.get_LLRs();
        obj.get_ROCs();
        obj.get_AUCs();
    end
    function obj=init_OUT(obj)
        obj.OUT.LLRb=struct();
        obj.OUT.LLR_IN=struct();
        obj.OUT.PC=struct();
    end
    function obj=get_IN(obj,measDist)
        if ~isa(measDist,'measDist')
            error('IN requires a measDist');
        end
        obj.IN=measDist.y;
        obj.INname=measDist.name;
        obj.get_LLRs_IN();
    end
    function obj = get_names(obj)
        obj.names=cell(obj.n,1);
        for i = 1:obj.n
            obj.names{i}=obj.mDists{i}.name;
        end
    end
    function obj = get_ctrs(obj)
        obj.ctrs=cell(obj.n,1);
        for i = 1:obj.n
            obj.x{i}=obj.mDists{i}.x;
            obj.ctrs{i}=obj.mDists{i}.ctrs;
        end
    end
    function obj = get_y(obj,n)
        %for i = 1:obj.n
        %    obj.y{i}=obj.mDists{i}.y;
    %    ynum(i)=numel(obj.y{i});
        %end

        if ~exist('n','var') || isempty(n)
            n=10000;
        end
        for i = 1:obj.n
            obj.y{i}=obj.mDists{i}.dst.rnd([n,1]);
        end
    end
    function obj = get_conditionals(obj)
        obj.fAA=cell(obj.n,1);
        for i = 1:obj.n
            obj.get_fAA(i);

            obj.get_fAB(i);

            obj.get_fBA(i);

            obj.get_fBB(i);

        end
    end
    function obj=get_A(obj,ind)
        obj.A =obj.y{ind};
    end
    function obj=get_B(obj,ind)
        obj.B=[];
        for i = 1:obj.n
            if i == ind
                continue
            end
            obj.B=[obj.B; obj.y{i}(:)];
        end
    end
    function get_fA(obj,ind)
        obj.xA=obj.x{ind};
        obj.fA=obj.mDists{ind}.pdf(:);
        obj.fA=obj.fA./sum(obj.fA(:));
    end
    function obj=get_fB(obj,ind)
        fB=[];
        xB=[];
        for i = 1:obj.n
            if i == ind
                continue
            elseif isempty(fB)
                fB=obj.mDists{i}.pdf(:);
                obj.xB=obj.x{i};
            else
                fB=[fB + obj.mDists{i}.pdf(:)];
            end
        end
        obj.fB=fB;
        obj.fB=fB./sum(fB(:)); % normalize B
    end
    function out = get_fAA(obj,ind)
    %likelihood of dataA given distribution A
    % USE ALL DATA
        obj.get_fA(ind);
        obj.get_A(ind);

        out=obj.interp_fun(obj.xA,obj.fA,obj.A,obj.method);
        %out=out/sum(out(:));
        obj.fAA{ind}=out;
    end
    function obj = get_fBA(obj,ind)

        obj.get_fB(ind);
        obj.get_A(ind);

        out=obj.interp_fun(obj.xB,obj.fB,obj.A,obj.method);
        obj.fBA{ind}=out;
    end
    function obj = get_fAB(obj,ind)
        obj.get_fA(ind);
        obj.get_B(ind);

        out=obj.interp_fun(obj.xA,obj.fA,obj.B,obj.method);
        %out=out/sum(out(:));
        obj.fAB{ind}=out;
    end
    function obj = get_fBB(obj,ind)
        obj.get_fB(ind);
        obj.get_B(ind);

        out=obj.interp_fun(obj.xB,obj.fB,obj.B,obj.method);
        %out=out/sum(out(:));
        obj.fBB{ind}=out;
    end
    function obj = get_LLRs(obj)
        obj.LLRa=cell(obj.n,1);
        obj.LLRb=cell(obj.n,1);
        for i = 1:obj.n
            obj.get_LLR(i);
        end
    end
    function obj = get_LLR(obj,ind)
        obj.OUT.LLRa{ind}=log(obj.fAA{ind}./obj.fBA{ind});
        obj.OUT.LLRb{ind}=log(obj.fAB{ind}./obj.fBB{ind});
        ctrs=Hist.binWidthsFD([obj.OUT.LLRa{ind}; obj.OUT.LLRb{ind}]);
        [LLRa,obj.LLRaCtrs{ind}]=hist(obj.OUT.LLRa{ind},ctrs);
        [LLRb,obj.LLRbCtrs{ind}]=hist(obj.OUT.LLRb{ind},ctrs);
        LLRa=LLRa./sum(LLRa(:));
        LLRb=LLRb./sum(LLRb(:));
        obj.LLRa{ind}=LLRa;
        obj.LLRb{ind}=LLRb;
    end
    function get_ROCs(obj)
        obj.fPrs=cell(obj.n,1);
        obj.tPrs=cell(obj.n,1);
        for i = 1:obj.n
            obj.get_ROC(i);
        end
    end
    function [fPr,tPr]=get_ROC(obj,ind)
        A=obj.LLRa{ind};
        B=obj.LLRb{ind};

        x=union(obj.LLRaCtrs{ind},obj.LLRbCtrs{ind});
        fPr=zeros(length(x),1);
        tPr=zeros(length(x),1);

        for i=1:length(x)
            c=x(i);
            tP=sum(A(x>=c));
            fP=sum(B(x>=c));

            fN=sum(A(x<c));
            tN=sum(B(x<c));

            tPr(i)=tP./(tP+fN);
            fPr(i)=fP./(fP+tN);
        end
        if tPr(end)~=0 || fPr(end)~=0
            tPr=[tPr; 0];
            fPr=[fPr; 0];
        end
        obj.fPrs{ind}=fPr;
        obj.tPrs{ind}=tPr;
    end
    function obj=get_AUCs(obj)
        obj.AUCs=zeros(obj.n,1);
        for i = 1:obj.n
            obj.get_AUC(i);
        end
    end
    function out=get_AUC(obj,ind)
        [fPr,i]=sort(obj.fPrs{ind});
        tPr=obj.tPrs{ind};
        tPr=tPr(i);
        out=trapz(fPr,tPr);
        obj.AUCs(ind)=out;
    end
%%%%%%
    function obj=get_LLRs_IN(obj)
        obj.LLRc=cell(obj.n,1);
        obj.OUT.LLR_IN=cell(obj.n,1);
        obj.OUT.PC=zeros(obj.n,1);
        for i = 1:obj.n
            obj.get_LLR_IN(i);
            obj.get_PC(i);
        end
    end
    function obj = get_PC(obj,ind)
        c=obj.LLRcCtrs{ind}>=0;
        obj.OUT.PC(ind)=sum(obj.LLRc{ind}(c))./sum(obj.LLRc{ind}(:));
    end
    function obj=get_LLR_IN(obj,ind)
        obj.get_fA(ind);
        obj.get_fB(ind);
        A=obj.interp_fun(obj.xA,obj.fA,obj.IN,'pchip');
        B=obj.interp_fun(obj.xB,obj.fB,obj.IN,'pchip');
        LLR=log(A./B);
        obj.OUT.LLR_IN{ind}=LLR;

        A
        ctrs=Hist.binWidthsFD(LLR);
        [LLR]=hist(LLR,ctrs);
        LLR=LLR./sum(LLR(:));

        obj.LLRc{ind}=LLR;
        obj.LLRcCtrs{ind}=ctrs;
    end
%% PLOT HELPERS

    function name=get_name(obj,ind)
        if isempty(obj.names{ind})
            name=Num.alphaU(ind);
        else
            name=obj.names{ind};
        end
    end

    function leg=get_legend_conditional(obj,ind)
        name=obj.get_name(ind);
        leg{1}=['p(data ' name ' | ' name ')'];
        leg{2}=['p(not data ' name ' | ' name ')'];
        leg{3}=['p(data ' name ' | not ' name ')'];
        leg{4}=['p(not data ' name ' | not ' name ')'];
        legend(leg{:});
    end

    function leg=get_legend_LLR(obj,ind)
        name=obj.get_name(ind);
        leg{1}=['p( LLR | ' name ' )'];
        leg{2}=['p( LLR | ~' name ' )'];
        legend(leg{:});
    end

%% PLOT

    function plot(obj)
        figure(obj.fignum)
        hold off;
        obj.plot_measurements_same();

        figure(obj.fignum+1)
        hold off;
        obj.plot_conditionals(); hold off

        figure(obj.fignum+2)
        hold off;
        obj.plot_LLRs();

        figure(obj.fignum+3)
        hold off;
        obj.plot_ROCs();
    end

    function plot_measurements_same(obj,bNoFormat)
        if ~exist('bNoFormat','var')
            bNoFormat=0;
        end
        m=[];
        for i = 1:length(obj.n)
            m=[m obj.mDists{i}.pdf];
        end
        y=max(m(:))*.05;
        for i = 1:obj.n
            obj.plot_measurement(i,i,y,bNoFormat); hold on
            if i == 1 && ~bNoFormat
                Fig.format('','Probability','Measurement Distributions');
            elseif ~bNoFormat
                Fig.format('','Probability','');
            end
        end
        if ~bNoFormat
            Fig.format('x','Probability','Measurement Distributions');
        end
    end

    function plot_measurements(obj,bNoFormat)
        if ~exist('bNoFormat','var')
            bNoFormat=0;
        end
        m=[];
        for i = 1:length(obj.n)
            m=[m obj.mDists{i}.pdf];
        end
        y=max(m(:))*.05;
        for i = 1:obj.n
            subplot(obj.n,1,i); hold on
            obj.plot_measurement(i,i,y,bNoFormat);
            if i == 1 && ~bNoFormat
                Fig.format('','Probability','Measurement Distributions');
            elseif ~bNoFormat
                Fig.format('','Probability','');
            end
        end
        Fig.format('x','Probability','Measurement Distributions');
    end

    function plot_measurement(obj,ind,i,y,bNoFormat)
        if ~exist('bNoFormat','var')
            bNoFormat=0;
        end

        if ~exist('i','var') || isempty(i)
            i=1;
        end
        if ~exist('y','var') || isempty(y)
            y=max(obj.mDists{ind}.pdf(:).*.05);
        end
        obj.mDists{ind}.plot_all(i,y);
        if ~bNoFormat
            Fig.format('x','Probability','Measurement Distributions');
        end
    end

    function plot_conditionals(obj)
        for i = 1:obj.n
            subplot(obj.n,1,i); hold off
            obj.plot_conditional(i);
            if i == 1
                Fig.format('','Probability','Conditional Distributions');
            else
                Fig.format('','Probability');
            end
        end
        Fig.format('y','Probability','Conditional Distributions');
    end
    function plot_conditional(obj,ind,bNoFormat)
        if ~exist('bNoFormat','var')
            bNoFormat=0;
        end
        plot(obj.fAA{ind}); hold on
        plot(obj.fBA{ind});
        plot(obj.fAB{ind});
        plot(obj.fBB{ind});
        obj.get_legend_conditional(ind);
        if ~bNoFormat
            Fig.format('y','Probability','Conditional Distributions');
        end
    end
    function plot_LLRs(obj)
        for i = 1:obj.n
            subplot(obj.n,1,i); hold off
            obj.plot_LLR(i);
            Fig.format('','LLR');
            if i ==1
                Fig.format('','Ratio','');
            else
                Fig.format('','LLR');
            end
        end
        Fig.format('LLR','Probability');
    end
    function plot_LLR(obj,ind,bNoFormat)
        if ~exist('bNoFormat','var')
            bNoFormat=0;
        end
        plot(obj.LLRaCtrs{ind},obj.LLRa{ind}); hold on
        plot(obj.LLRbCtrs{ind},obj.LLRb{ind});
        obj.get_legend_LLR(ind);
        if ~bNoFormat
            Fig.format('LLR','Probability','');
        end
    end
    function plot_ROCs(obj,bNoFormat)
        if ~exist('bNoFormat','var')
            bNoFormat=0;
        end
        for i = 1:obj.n
            subplot(obj.n,1,i); hold off
            obj.plot_ROC(i);
            Fig.format('','Hit Rate (TP)',obj.names{i});
        end
        if ~bNoFormat
            Fig.format('False Alarm Rate (FP)','Hit Rate (TP)');
        end
    end
    function plot_ROC(obj,ind,bNoFormat)
        if ~exist('bNoFormat','var')
            bNoFormat=0;
        end
        plot(obj.fPrs{ind},obj.tPrs{ind},'LineWidth',2);hold on
        plot([0,1],[0,1],'k--');
        xlim([0,1]);
        ylim([0,1]);
        axis square;
        if ~bNoFormat
            Fig.format('False Alarm Rate (FP)','Hit Rate (TP)');
        end
        t=['AUC = ' num2str(sprintf('%.3f',round(obj.AUCs(ind),3))) ];
        text(.25,.1,t,'FontSize',18);
    end
    function plot_IN(obj,ind,bNoFormat)
        if ~exist('bNoFormat','var')
            bNoFormat=0;
        end
        plot(obj.LLRcCtrs{ind},obj.LLRc{ind}); hold on
        m=max(obj.LLRc{ind})*1.05;
        plot([0 0],[0 m],'k:');
        ylim([0 m]);
        if ~bNoFormat
            Fig.format('LLR','Probability',[num2str(obj.OUT.PC{ind}.*100) '% chosen']);
        end
        obj.get_legend_LLR_IN(ind);
        hold off;
    end
%%%%%%
end
methods(Hidden = true)
    function ha=init_panel(obj)
        ha=panel();
        ha.pack('h',{1/obj.n []});
        ha(1).pack(obj.n);
        ha.fontsize=18;
        ha.margin=[25 20 2 2];
        ha.de.margin=2;
        ha.select('all');
    end
    function out=interp_fun(obj,x,y,q,method)
        try
            out=interp1(x,y,q,method);
        catch ME
            rethrow(ME);
        end
        i=(out<=0 | isnan(out));
        if sum(i)>0
            val=min(out(~i));
            out(i)=val;
            %out(i)=interp1(x,y,q(i),'nearest',val);
        end

        i=isnan(out);
        if sum(i)>0
            val=max(out(~i));
            out(i)=val;
            %out(i)=interp1(x,y,q(i),'nearest',val);
        end
    end
end
end
