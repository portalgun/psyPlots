classdef measDist < handle
properties
    A
    B
    x
    y
    h
    pdf
    ctrs
    method
    dst
    name
    nanind
end
methods
    function obj=measDist(name,x,y,binCtrs,method,bRmOutliers)
    %function obj=measDist(name,x,y,binCtrs,method,bRmOutliers)

    % IF WORKING WITH MULTIPLE MEASDIST, CHOOSE X MANUALLY
        obj.name=name;
        obj.nanind=isnan(y);
        obj.y=y;
        if ~isvar('bRmOutliers')
            bRmOutliers=0;
        end
        if bRmOutliers
            y=rmoutliers(y(:),'percentile',[2.5 97.5]);
        else
            y=y(:);
        end
        if ~isvar('x')
            obj.x=1:size(y,1);
        else
            obj.x=x;
        end
        if ~isvar('binCtrs')
            obj.ctrs=bin_widths_FD(obj.y);
        else
            obj.ctrs=binCtrs;
        end
        if ~isvar('method')
            obj.method='normal';
        else
            obj.method=method;
        end
        obj.get_hist();
        obj.get_pdf();

    end
    function obj=get_hist(obj)
        obj.h=hist(obj.y,obj.ctrs);
        %obj.h=h./sum(h(:));
    end
    function obj=get_pdf(obj)
        switch obj.method
        case {'spline','linear','nearest','pchip'}
            obj.get_pdf_interp();
        otherwise
            y=obj.y;
            y(obj.nanind)=[];
            dst=dstC_selector(obj.method,y);
            dst.get_fit();
            dst.get_pdf(obj.x);
            obj.pdf=dst.pdf;
            if nansum(obj.pdf) == 0
                warning('Bad fit')
            end
            obj.dst=dst;
        end
    end
    function obj=get_pdf_interp(obj)
        out=interp1(obj.ctrs,obj.h,obj.x,obj.method,0);
        ind=out(out<0);
        if sum(ind)>=0
            interp1(obj.ctrs,obj.h,obj.x(ind),'nearest',0);
        end
        out=out./sum(out(:));
        obj.pdf=out;
    end
    function out=get_pdf_at(obj,X)
        out=interp1(obj.x,obj.pdf,X,obj.method,0);
        %out(out<0)=0;
        %out=out./sum(out(:));
    end
    function plot_h(obj)
        figure(832)
        plot(obj.ctrs,obj.h,'k')
        formatFigure('x','count')
    end
    function plot(obj)
        figure(831)
        hold off
        plot(obj.x,obj.pdf,'k'); hold on
        %plot(obj.y,zeros(size(obj.y)),'.')
        formatFigure('x','p')
    end
    function plot_all(obj,i,y)
        if ~exist('i','var') || isempty(i)
            i=1;
        end
        if ~exist('y','var') || isempty(y)
            y=max(obj.pdf(:).*.05);
        end
        h=obj.h./sum(obj.h(:))./mean(diff(obj.ctrs(:))).*mean(diff(obj.x(:)));
        plot(obj.ctrs,h,':'); hold on
        plot(obj.x,obj.pdf,'Color',lastColor); hold on
        %plot(obj.x,obj.pdf); hold on
        plot(obj.y,i*y*ones(size(obj.y)),'.','Color',lastColor);
    end
end
end
