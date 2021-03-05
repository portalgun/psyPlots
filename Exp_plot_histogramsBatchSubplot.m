function []=Exp_plot_histogramsBatchSubplot(obj,opts,figNum)
figure(figNum)

if ~isfield(opts,'subPlotXSub')
    opts.subPlotXsub=1;
end

flds=fieldnames(opts.F2);
for i = 1:length(flds)
    fld=flds{i};
    opts.(fld)=opts.F2.(fld);
end
rmfield(opts,'F2');

counts=zeros(size(obj.subsAllLimit,1),opts.nbins);
ctrs=zeros(size(obj.subsAllLimit,1),opts.nbins);
xl=zeros(size(obj.subsAllLimit,1),2);
for i = 1:size(obj.subsAllLimit,1)
    sub=obj.subsAllLimit(i,:);

    %GET DATA
    SS=obj.select(sub);
    data=SS.(opts.field);

    %GET CTRS AND COUNTS, BUT DON'T PLOT
    [counts(i,:),ctrs(i,:),xl(i,:)]=plot_histograms(data,[],[],opts,1,figNum,1);
end

% FORMAT
opts.xlims=[minall(xl) maxall(xl)];
opts.ylims=[0 maxall(counts)];
opts.figNum=figNum;
opts.bSubPlot=1;
opts.xtitl=opts.field;

bReset=1;
for i = 1:size(obj.subsAllLimit,1)
    subplotSubsExps(figNum,sub,obj.subsAllLimit,opts.subPlotXSub);
    sub=obj.subsAllLimit(i,:);
    [~, opts.color, opts.titl]=def_Exp_plotParamsByExp(obj.shapeFunc,sub);
    plot_histograms([],counts(i,:),ctrs(i,:),opts,0,figNum,bReset);
    opts.titl=[opts.titl ' ' SS.Subjs];
    gAxis(opts);
    bReset=0;
end
