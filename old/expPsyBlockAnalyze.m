function [] = expPsyBlockAnalyze(prjCodesAll,expIDsAll,SUBJ)
%loop accross
% train
% prelim
% experiment
% pass
fsLbl=12;
fsTck=12;
xtckSD=0;
ytckSD=1;
FILES=cell(0,1);
%color=[0,0,1];
%shape=markerDefs(subj);
% -------------------------------------------------------------------------------
%init field if doesn't exist

%   E=loadExpStruct(prjCode,expID);
%   nStdX=length(E.stdXunqAll);
%
%   disp(newline);
%   for i = 1:nStdX
%       disp([ num2str(i) ' ' num2str(E.stdXunqAll(i)) ]);
%   end
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
if ~iscell(prjCodesAll)
    prjCodesAll={prjCodesAll};
end
if ~iscell(expIDsAll)
    expIDsAll={expIDsAll};
end
if ~iscell(SUBJ)
    SUBJ={SUBJ};
end

bStart=1;

%PRE PRJCODES+EXPIDS IN ORDER TO LIMIT REDUNDANCY
passListAll=cell(length(prjCodesAll),1);
for W=1:length(prjCodesAll)
    prjCode=prjCodesAll{W};
    expID=expIDsAll{W};

    E=loadExpStruct(prjCode,expID);
    %Put all passes in list, with desired pass first
    if isfield(E,'link') && isfield(E.link,'pass') && ~isempty(E.link.pass)
        passListAll{W}=cell(length(E.link.pass)+1,1);
        passListAll{W}{1:length(E.link.pass),1}=E.link.pass{:};
        passListAll{W}{end,1}=[prjCode ' ' expID];
        passListAll{W}=sort(passListAll{W});
    else
        passListAll{W}=cell(1,1);
        passListAll{W}{1}=[prjCode ' ' expID];
    end

end

AllFirst=flip(cellfun(@(x) x(1),passListAll));
AllLast=cellfun(@(x) x(end),passListAll);


%unique to each set
Alltmp=setxor(AllFirst,AllLast);
bidx=ismember(Alltmp,AllLast);

idx=[];
for i = 1:length(Alltmp)
    idx(end+1)=find(contains(AllLast,Alltmp{i}),1);
end

AllLast=unique(flip(Alltmp(bidx)));
AllFirst=unique(AllFirst(idx));

prjCodesAll=cell(length(AllFirst),1);
expIDsAll=cell(length(AllFirst),1);
prjCodesTit=cell(length(AllFirst),1);
expIDsTit=cell(length(AllFirst),1);
for i = 1:length(AllFirst)
    strs=strsplit(AllLast{i});
    prjCodesAll{i}=strs{1};
    expIDsAll{i}=strs{2};

    %FOR TITLES
    strs=strsplit(AllFirst{i});
    prjCodesTit{i}=strs{1};
    expIDsTit{i}=strs{2};
end
% -------------------------------------------------------------------------------

for W=1:length(prjCodesAll)
    prjCode=prjCodesAll{W};
    expID=expIDsAll{W};

    superTitle=[ prjCodesTit{W} '_' expIDsTit{W} ];
    superTitle=strrep(superTitle,'_','-');
    [b,e]=regexp(superTitle,'_pass[0-9]+');
    if isempty(b)
        superTitlePRE=superTitle;
    elseif length(superTitle) > e
        superTitlePRE=[superTitle(1:b-1) superTitle(e+1:end)];
    else
        superTitlePRE=superTitle(1:b-1);
    end

    E=loadExpStruct(prjCode,expID);
    if isfield(E,'link') && isfield(E.link,'pass') && ~isempty(E.link.pass)
        passList=E.link.pass;
        passList{end+1}=[prjCode ' ' expID];
        passList=sortrows(passList,2);
    else
        passList={[prjCode ' ' expID]};
    end
    [prjCodes,expIDs]=cellfun(@(s) strtok(s,' '),passList(:), 'UniformOutput',false);
    for i = 1:length(expIDs)
        expIDs{i}=strrep(expIDs{i},' ','');
    end


    fnum=nFn;
    for Z = 1:length(SUBJ)
        subj=SUBJ{Z};
        shape=markerDefs(subj);

        if length(SUBJ)==1
            superTitle=[ superTitlePRE ' ' subj ];
        else
            superTitle=superTitlePRE;
        end
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        errorflag=0;
        while true
            if errorflag==1 || errorflag==2
                response=cellstr(response');
                disp('Invalid response(s): ');
                disp(['    ' response{idxInv}]);
            else
                E=loadExpStruct(prjCode,expID);
                nStdX=length(E.stdXunqAll);

                if bStart==1
                    disp(newline);
                    for i = 1:nStdX
                        disp([ num2str(i) ' ' num2str(E.stdXunqAll(i)) ]);
                    end
                end
            end
            if bStart==1
                response=input('Enter space seperated corresponding numbers to standards you want to check, or ''a'' for all: ','s');
                if strcmp(response,'a') || strcmp(response,'A') || strcmp(response,'all') || strcmp(response,'ALL') || strcmp(response,'All')
                    response=num2str(1:nStdX);
                end
                if strcmp(response,'')
                    return
                end

                idxInv=valChars(response,[' ';numV]);
                if ~isempty(idxInv)
                    errorflag=1;
                    continue
                end
                idxInv=valRange(response,1:nStdX);
                if ~isempty(idxInv)
                    errorflag=2;
                    continue
                end

                response=strsplit(response)';
                response(strcmp('',response))=[];
                break
            else
                break
            end
        end
        if bStart==1
            bStart=0;
            idxStd=str2double(response);
            stdX=E.stdXunqAll(idxStd);
            blkPerExp=E.trlPerExp/E.trlPerBlk;
        end

        % -------------------------------------------------------------------------------


        %GET BLOCK INDICES & FILES
        Dtypes={'train','prelim','experiment'};

        count=0;
        blkidx=zeros(length(prjCodes),3,length(stdX),blkPerExp);
        SS=cell(length(prjCodes),3,length(stdX),blkPerExp,10);
        EE=cell(length(prjCodes),1);
        ESubj=cell(length(prjCodes),3);
        CurReRun=zeros(length(prjCodes),length(stdX),blkPerExp);
        for h = 1:length(prjCodes)
            prjCode=prjCodes{h};
            expID=expIDs{h};
            E=loadExpStruct(prjCode,expID);
            EE{h}=E;
            for d = 1:length(Dtypes)

                Dtype=Dtypes{d};
                if strcmp(Dtype,'experiment')
                    stdfldsAll=fieldnames(E.SubjData.(subj).stdX);
                    stdflds=stdfldsAll(idxStd);
                    ESubj{h,d}=E.SubjData.(subj);
                    fdir=strcat(E.serverDataDir,subj,filesep);
                elseif strcmp(Dtype,'prelim')
                    if ~isfield(E.SubjData.(subj),'prelim')
                        continue
                    end
                    stdfldsAll=fieldnames(E.SubjData.(subj).prelim.stdX);
                    stdflds=stdfldsAll(idxStd);
                    ESubj{h,d}=E.SubjData.(subj).prelim;
                    fdir=strcat(E.serverDataDir,subj,filesep,'prelim',filesep);
                elseif strcmp(Dtype,'train')
                    if ~isfield(E.SubjData.(subj),'prelim')
                        continue
                    end
                    stdfldsAll=fieldnames(E.SubjData.(subj).train.stdX);
                    stdflds=stdfldsAll(idxStd);
                    ESubj{h,d}=E.SubjData.(subj).train;
                    fdir=strcat(E.serverDataDir,subj,filesep,'train',filesep);
                end

                %INIT FLAG FIELD IF DOESNT EXIST
                for i = 1:length(E.stdXunqAll)
                    if ~isfield(ESubj{h,d}.stdX.(stdfldsAll{i}),'flag')
                        ESubj{h,d}.stdX.(stdfldsAll{i}).flag=zeros(blkPerExp,1);
                    end
                end

                for i = 1:length(stdX)
                    files=strcat(fdir,ESubj{h,d}.stdX.(stdflds{i}).data);
                    filesold=strcat(fdir,'old',filesep,ESubj{h,d}.stdX.(stdflds{i}).data);
                    for j = 1:blkPerExp
                        Split=strsplit(files{j},'-');
                        baseidx=regexp(files{j},'-[0-9]-[0-9]\.mat');
                        ReExpidx=regexp(files{j},'-[0-9]\.mat');
                        filebase=filesold{j}(1:baseidx);
                        CurTmp=Split{end-1};
                        CurReRun(h,i,j)=str2double(CurTmp)+1; %starts with 1 instead of zero
                        CurReExp=files{j}(ReExpidx+1:end-4);
                        %if CurRun
                        if ~isempty(ESubj{h,d}.stdX.(stdflds{i}).data{j}) && ESubj{h,d}.stdX.(stdflds{i}).blk(j)==1 && CurReRun(h,i,j)<=10
                            count=count+1;
                            blkidx(h,d,i,j)=count;
                            if ~ismember(files{j},FILES)
                                SS{h,d,i,j,CurReRun(h,i,j)}=load(files{j});
                                FILES{end+1}=files{j};%ADD TO LOADED LIST
                            end
                            for c=1:10
                                try
                                    if ~strcmp(num2str(c),num2str(CurReRun(h,i,j)))
                                        possFile=strcat(filebase, sprintf('%03d',j) ,'-', num2str(c-1), '-', CurReExp, '.mat'); %block - trial rerun (10 max) - experiment rerun; c-1 because files start at zero
                                        if ~ismember(possFile,FILES)
                                            try
                                                SS{h,d,i,j,c}=load(possFile);
                                            catch
                                                continue
                                            end
                                            FILES{end+1}=possFile;%ADD TO LOADED LIST
                                        end
                                    end
                                catch
                                    break
                                end
                            end
                        end
                    end
                end

            end
        end

        % -------------------------------------------------------------------------------
        figure(fnum);
        hold on
        if length(SUBJ) == 1
            SPC=ceil(sqrt(length(stdX)));
            SPR=floor(sqrt(length(stdX)));
        else
            SPC=length(stdX);
            SPR=length(SUBJ);
        end

        fprintf('Please wait. \nBootstrapping ')

        tCI=cell(length(prjCodes),3,length(stdX),blkPerExp,10);
        tFit=zeros(length(prjCodes),3,length(stdX),blkPerExp,10);
        pcData=cell(length(prjCodes),3,length(stdX),blkPerExp,10);
        Count=zeros(length(prjCodes),length(stdX));
        for i = 1:length(stdX)
        for j = 1:blkPerExp
            for h=1:length(prjCodes)
            for d=1:3
            for c=1:10
                if blkidx(h,d,i,j)~=0
                    S=SS{h,d,i,j,c};
                    if ~isempty(S) && isfield(S,'RcmpChosen')
                        [~,~,~,tCI{h,d,i,j,c},~,~,~,~]  = psyfitgengaussBootstrap(S.stdX, S.cmpX, S.RcmpChosen,[],[],1,1.36,2,100,68,[],0,1);
                        %fits
                        [~,~,~,tFit(h,d,i,j,c),pcData{h,d,i,j,c},~,~] = psyfitgengauss(S.stdX, S.cmpX, S.RcmpChosen, [],[],1,1.36,2,0);
                        Count(h,i)=Count(h,i)+1;
                    else
                        %tCI{h,d,i,j,c}=[-10,10];
                        %tFit(h,d,i,j,c)=0;
                        %pcData{h,d,i,j,c}=zeros(length(E.ncmp),1);
                        %ESubj{h,d}.stdX.(stdflds{i}).flag(j)=1;
                        %saveExpStruct(E);

                        tCI{h,d,i,j,c}=[];
                        pcData{h,d,i,j,c}=[];
                        tFit(h,d,i,j,c)=NaN;
                        pcData{h,d,i,j,c}=[];
                    end
                    fprintf('.')
                else
                    continue
                end
            end
            end
            end

        end
        end

        %Bootstrap mean
        figure(fnum);
        siz=[length(prjCodes),3,length(stdX),blkPerExp,10];
        for i = 1:length(stdX)


            flagidx=find(ESubj{h,d}.stdX.(stdflds{i}).flag==1);
            subplot(SPR,SPC,i+(Z-1)*length(stdX))
            hold on
            COUNT=1;
            for j = 1:blkPerExp
                gInd=cell(0);
                COUNTb=COUNT;

                for h=1:length(prjCodes)
                for d=1:3
                for c=1:10
                    if isempty(tCI{h,d,i,j}) || isnan(tFit(h,d,i,j,c))
                        continue
                    end
                    %NOTE PLOT MEAN POINTS
                    %FLAG = RED
                    if ismember(j,flagidx) && c==CurReRun(h,i,j) && d==3
                        plot(COUNT,tFit(h,d,i,j,c),'rs','linewidth',2,'MarkerFaceColor','w');
                    %TRAIN = BLUE
                    elseif d==1
                        plot(COUNT,tFit(h,d,i,j,c),'bs','linewidth',2,'MarkerFaceColor','w');
                    %PRELIM = CYAN
                    elseif d==2
                        plot(COUNT,tFit(h,d,i,j,c),'cs','linewidth',2,'MarkerFaceColor','w');
                    %RERUN = MAGENTA
                    elseif c~=CurReRun(h,i,j)
                        plot(COUNT,tFit(h,d,i,j,c),'ms','linewidth',2,'MarkerFaceColor','w');
                    %NORMAL = BLACK
                    else
                        plot(COUNT,tFit(h,d,i,j,c),'ks','linewidth',2,'MarkerFaceColor','w');
                    end
                    text(COUNT,tFit(h,d,i,j,c),num2str(h),'VerticalAlignment','top','HorizontalAlignment','center','FontSize',10);
                    gInd{end+1}=[h,d,i,j,c];
                    if COUNT==1
                        first=sub2ind(siz,h,d,i,j,c);
                    end
                    last=sub2ind(siz,h,d,i,j,c);
                    COUNT=COUNT+1;
                end
                end
                end

                ind=cell(length(gInd),1);
                for g = 1:length(gInd)
                    ind{g}=sub2ind(siz,gInd{g}(1),gInd{g}(2),gInd{g}(3),gInd{g}(4),gInd{g}(5));
                end
                ind=[ind{:}];

                tCIup=zeros(length(ind),1);
                tCIdown=zeros(length(ind),1);
                for t = 1:length(ind)
                    tCIup(t)  =tCI{ind(t)}(1);
                    tCIdown(t)=tCI{ind(t)}(2);
                end

                hold on
                if length(COUNTb:COUNT-1)==1
                    %NOTE PLOT ERROR BARS
                    errorbar(COUNTb,(tCIup+tCIdown)/2,tCIdown',tCIup');
                else
                    %NOTE PLOT FILLED ERROR
                    plotfillederror(COUNTb:COUNT-1,tCIdown',tCIup',[.2,.2,.2],0);
                end

                if mod(length(COUNTb:COUNT-1),2)==0 && exist('last','var')==1 && exist('first','var')==1
                    text((COUNTb+COUNT-1)/2,(tFit(last)+tFit(first))/2,num2str(j),'VerticalAlignment','bottom','HorizontalAlignment','center','FontSize',14,'FontWeight','Bold');
                elseif ~isempty(tCIup) && ~isempty(tCIdown)
                    ind2=ceil(length(COUNTb:COUNT-1)/2);
                    y=(tCIup(ind2)+tCIdown(ind2))/2*1.1;
                    text((COUNTb+COUNT-1)/2,y,num2str(j),'VerticalAlignment','bottom','HorizontalAlignment','center','FontSize',14,'FontWeight','Bold');
                end

            end

            sX=stdX(i);
            axis normal
            if i == 1 && length(SUBJ)>1
                ylbl=({subj;'Threshold (arcmin)'});
            else
                ylbl=('Threshold (arcmin)');
            end
            if i == ceil(length(stdX)/2) && length(SUBJ)>1
                xlbl=('Block Number');
            else
                xlbl=('Block Number');
            end
            ttl= [ num2str(sX) ' arcmins' ];
            xtck=1000;
            formatFigure(xlbl,ylbl,ttl,0,0,fsLbl,fsTck,xtck,[],xtckSD,ytckSD);
            xlim([0,COUNT]);
            hold off
        end

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    end
    suptitle(superTitle);
end
