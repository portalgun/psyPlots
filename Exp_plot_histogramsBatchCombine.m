function fig=Exp_plot_histogramsBatchCombine(obj,opts,figNum)
% SINGLE HISTOGRAM, DIFFERENT EXPERIMENTS WITH DIFFERENT LINE

fig=figure(figNum)
bReset=1;

subjAll=[];
subjSeen={''};
expSeen={''};
opts.titl='';

xlims=zeros(size(obj.subsAllLimit,1),2);
for i = 1:size(obj.subsAllLimit,1)
    sub=obj.subsAllLimit(i,:);

    %GET DATA
    SS=obj.select(sub);
    data=SS.Data.(opts.field);
    if isempty(data)
        continue
    end

    %GET COLORS ETC
    [~, opts.color, titl]=def_Exp_plotParamsByExp(obj.shapeFunc,sub);

    %CREATE TITLE
    if ~contains(expSeen,titl)
        opts.titl=[ opts.titl ', ' titl ];
        expSeen{end+1}=titl;
    end
    %GET SUBJ DEFINED LINES
    if ~contains(subjSeen,SS.Subjs)
        subjAll=[ subjAll ', ' SS.Subjs ];
        subjSeen{end+1}=SS.Subjs;
    end
    opts.lines=def_Exp_subj(SS.Data.subjName);

    %PLOT
    [~,~,xlims(i,:)]=plot_histograms(data,[],[],opts,0,figNum,bReset);
    bReset=0;
end

%FORMAT
if ~isfield(opts,'xlims')
    opts.xlims(1)=minall(xlims);
    opts.xlims(2)=maxall(xlims);
end

opts.titl(1)=[];
opts.titl(end)=[];
subjAll(1)=[];
opts.xtitl=opts.field;
opts.figNum=figNum;
opts.titl=[opts.titl ';' subjAll];

gAxis(opts);
gFigure(opts);

save_plot(fig,opts);
