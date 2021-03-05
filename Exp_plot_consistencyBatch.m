function [figsAll,inds] = Exp_plot_consistencyBatch(obj,opts,passORsubjORexp,figNum)
bReset=1;
if strcmp(passORsubjORexp,'subj')
    n=obj.noSubjLimit;
    v=obj.limits{1};
elseif strcmp(passORsubjORexp,'pass')
    n=obj.noPassLimit;
    v=obj.limits{end};
end
figsAll=cell(size(n,1),1);
titlsAll=cell(size(n,1),1);
inds=cell(size(n,1),1);
for i = 1:size(n,1)
    if ~opts.hold || bReset
        figNum=figNum+1;
    end
    loop=n(i,:);

    [figsAll{i},titlsAll{i},inds{i},]=plot_consistency(obj,opts,v,loop,passORsubjORexp,figNum);
    if all(cellfun(@(x) isempty(x),titlsAll{i}))
        continue
    end
    opts.figNum=figNum;
    gFigure(opts);

    for j = 1:length(titlsAll{i})
        opts.titl=titlsAll{i}{j};
        save_plot(figsAll{i}{j},opts);
    end

    obj=obj.save_inds(passORsubjORexp,loop,inds{i});

    bReset=0;
end
