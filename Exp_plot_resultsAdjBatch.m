function plots = Exp_plot_resultsAdjBatch(obj,opts,figNum)
plots=cell(size(obj.subsAllLimit,1),1);
bReset=1;
lastsub=[];
if ~isempty(opts.holdType)
    subsAllLimit=sort_fun(obj.subsAllLimit,opts.holdType);
else
    subsAllLimit=obj.subsAllLimit;
end
for i = 1:size(obj.subsAllLimit,1)
    %SUB
    if i > 1
        lastsub=sub;
    end
    sub=subsAllLimit(i,:);

    holdCond=holdCond_fun(lastsub,sub,opts.holdType);

    %FIGNUM
    opts.hold;
    if ~opts.hold
        figNum=figNum+1;
    elseif opts.hold && i > 1 && holdCond
        figNum=figNum+1;
    end

    %GET SS
    opts.figNum=figNum;
    SS=obj.select(sub);
    if ~strcmp(opts.sortFld,'')
        SS=sortDataByField(SS.Data,opts.sortFld);
    end
    if isempty(SS.X) || (~isequal(length(SS.X),length(SS.Xhat)))
        continue
    end

    %SUBJ/EXP SPECIFICS
    [opts.shape,opts.color,titl]=def_Exp_plotParamsByExp(obj.shapeFunc,sub);
    if opts.hold && strcmp(opts.holdType,'Pass') && sub(end)==2
        opts.color='b';
    end

    % TITLE
    if opts.hold && strcmp(opts.holdType,'Subj')
        sb=join(obj.Subjs(obj.prjInd{1}),'-');
        sb=sb{1};
        opts.titl=[sb ' ' titl ' pass' num2str(sub(4))];
    else
        if size(SS.subjName,1)>1
            SS.subjName=cellstr(SS.subjName);
            SS.subjName=join(SS.subjName,'-');
            SS.subjName=SS.subjName{1};
        end
        opts.titl=[SS.subjName ' ' titl ' pass' num2str(sub(4))];
    end

    %PLOT
    fig=plot_resultsAdj(SS,opts);
    gAxis(opts);
    plots{i}=fig;
    bReset=0;

    %save_plot(fig,opts);
end
end

function subsAllLimit=sort_fun(subsAllLimit,holdType)
    switch holdType
    case 'Subj'
        order=[3 2 4 1];
    case 'Exp'
        order=[1 4 3 2];
    case 'Pass'
        order=[3 2 1 4];
    end
    subsAllLimit=sortrows(subsAllLimit,order);
end

function holdCond=holdCond_fun(lastsub,sub,holdType)
    switch holdType
        case 'Subj'
            %NEW FIGURE ON EXPS OR PASS
            holdCond=~isequal(lastsub(2:end),sub(2:end))
        case 'Exp'
            %NEW FIGURE ON SUBJ OR PASS
            holdCond=~isequal([lastsub(1) lastsub(end)],[sub(1) sub(end)]);
        case 'Pass'
            % NEW FIGURE ON SUBJ OR EXP
            holdCond=~isequal([lastsub(1:end-1)],[sub(1:end-1)]);
    end
end
