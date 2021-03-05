function []= expPsyCurveAll(prjCodes,expIDs,SUBJ,Dtype,bSubPlot)
%ALL BLOCKS

fsLbl=12;
fsTck=12;
xtckSD=2;
ytckSD=1;
colors=[[18,133,158];[49,97,107];[45,209,148];[215,103,100];[158,18,73];[18,133,158]]; %Compound
colors=colors ./ 255;
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
if ~iscell(prjCodes)
    prjCodes={prjCodes};
end
if ~iscell(expIDs)
    expIDs={expIDs};
end
if ~iscell(SUBJ)
    SUBJ={SUBJ};
end

for W=1:length(prjCodes)
    prjCode=prjCodes{W};
    expID=expIDs{W};

    superTitle=[ prjCode '_' expID ];
    superTitle=strrep(superTitle,'_','-');
    [b,e]=regexp(superTitle,'_pass[0-9]+');
    if isempty(b)
        superTitlePRE=superTitle;
    elseif length(superTitle) > e
        superTitlePRE=[superTitle(1:b-1) superTitle(e+1:end)];
    else
        superTitlePRE=superTitle(1:b-1);
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
        E=loadExpStruct(prjCode,expID);
        stdflds=fieldnames(E.SubjData.(subj).stdX);
        blkPerExp=E.trlPerExp/E.trlPerBlk;
        if strcmp(Dtype,'experiment')
            ESubj=E.SubjData.(subj);
            fdir=strcat(E.serverDataDir,filesep,subj,filesep);
        elseif strcmp(Dtype,'prelim')
            ESubj=E.SubjData.(subj).prelim;
            fdir=strcat(E.serverDataDir,filesep,subj,filesep,'prelim',filesep);
        elseif strcmp(Dtype,'train')
            ESubj=E.SubjData.(subj).train;
            fdir=strcat(E.serverDataDir,filesep,subj,filesep,'train',filesep);
        end
        % -------------------------------------------------------------------------------
        % SELECT STDS
        if Z==1 && W==1
            disp(newline);
            for i = 1:length(E.stdXunqAll)
                disp([ num2str(i) ' ' num2str(E.stdXunqAll(i)) ]);
            end

            errorflag=0;
            while true
                if errorflag==1 || errorflag==2
                    response=cellstr(response');
                    disp('Invalid response(s): ');
                    disp(['    ' response{idxInv}]);
                end
                response=input('Enter space seperated corresponding numbers to standards you want to check, or ''a'' for all: ','s');
                if strcmp(response,'a') || strcmp(response,'A') || strcmp(response,'all') || strcmp(response,'ALL') || strcmp(response,'All')
                    response=num2str(1:length(E.stdXunqAll));
                end
                if strcmp(response,'')
                    return
                end

                idxInv=valChars(response,[' ';numV]);
                if ~isempty(idxInv)
                    errorflag=1;
                    continue
                end
                idxInv=valRange(response,1:length(E.stdXunqAll));
                if ~isempty(idxInv)
                    errorflag=2;
                    continue
                end

                response=strsplit(response)';
                response(strcmp('',response))=[];
                break
            end
            idxStd=str2double(response);
            stdX=E.stdXunqAll(idxStd);
            nStds=length(stdX); %number of stds/fields
        end

        % -------------------------------------------------------------------------------
        %GET BLOCK INDICES & FILES
        count=0;
        count2=0;
        %files=cell(length(stdflds)*blkPerExp,1);
        files=cell(1);
        blkidx=zeros(length(E.stdXunqAll),blkPerExp);
        for i = 1:length(E.stdXunqAll)
            for j = 1:blkPerExp
                if ~isempty(ESubj.stdX.(stdflds{i}).data{j}) && ESubj.stdX.(stdflds{i}).blk(j)==1
                    count=count+1;
                    blkidx(i,j)=count;
                    strnum=fldstr2num(stdflds{i});
                    if ismember(strnum,stdX)
                        count2=count2+1;
                        files{count2,1}=strcat(fdir,ESubj.stdX.(stdflds{i}).data{j});
                    end
                end
            end
        end
        % -------------------------------------------------------------------------------

        XS=struct();
        for i = 1:length(files)
            S=load(files{i});
            str=strrep(num2str(S.stdXunq),'-','n'); %make stds into valid field names
            str=strrep(str,'.','_');               %...
            r=ceil(i/blkPerExp);
            c=mod(i,blkPerExp);
            if c==0
                c=blkPerExp;
            end
            %XXX prevent different cmps for same standard
            if isfield(XS,str) %if field has been created
                try
                    XS.(num2str(str)).RcmpChosen = [XS.(num2str(str)).RcmpChosen; S.RcmpChosen];
                    XS.(num2str(str)).stdX       = [XS.(num2str(str)).stdX; S.stdX];
                    XS.(num2str(str)).cmpX       = [XS.(num2str(str)).cmpX; S.cmpX];
                    XS.(num2str(str)).Rcorrect   = [XS.(num2str(str)).Rcorrect; S.Rcorrect];
                catch
                    disp(files{i});
                    if ESubj.stdX.(stdflds{r}).blk(c)~=0
                        ESubj.stdX.(stdflds{r}).flag(c)=1;
                        saveExpStruct(E);
                    end
                end
            else                 %if field has not been created
                try
                    XS.(num2str(str)).RcmpChosen = S.RcmpChosen;
                    XS.(num2str(str)).stdX       = S.stdX;
                    XS.(num2str(str)).cmpX       = S.cmpX;
                    XS.(num2str(str)).cmpXunq    = S.cmpXunq';
                    XS.(num2str(str)).Rcorrect   = S.Rcorrect;
                catch
                    disp(files{i});
                    if ESubj.stdX.(stdflds{r}).blk(c)~=0
                        ESubj.stdX.(stdflds{r}).flag(c)=1;
                        saveExpStruct(E);
                    end
                end
            end
        end

        %Get all fieldnames and sort
        fldNames=fieldnames(XS);

        % -------------------------------------------------------------------------------
        track=sum(blkidx,2);
        nnzeros=nnz(track);
        mCI    = cell(nnzeros,1);
        sCI    = cell(nnzeros,1);
        tCI    = cell(nnzeros,1);
        mFit   = cell(nnzeros,1);
        sFit   = cell(nnzeros,1);
        bFit   = cell(nnzeros,1);
        tFit   = zeros(nnzeros,1);
        pcData = cell(nnzeros,1);
        pcFit  = cell(nnzeros,1);

        fprintf('Please wait. \nBootstrapping ')
        count=0;
        for i = 1:nStds
            if track(i)~=0
                count=count+1;
                %confidence intervals
                [mCI{count},sCI{count},~,tCI{count},~,~,~,~]  = psyfitgengaussBootstrap(XS.(fldNames{count}).stdX,XS.(fldNames{count}).cmpX,XS.(fldNames{count}).RcmpChosen,[],[],1,1.36,2,100,68,[],0,1);
                %fits
                [mFit{count},sFit{count},bFit{count},tFit(count),pcData{count},pcFit{count},~] = psyfitgengauss(XS.(fldNames{count}).stdX,XS.(fldNames{count}).cmpX,XS.(fldNames{count}).RcmpChosen,[],[],1,1.36,2,0);
            end
        %                 psyfitgengaussBootstrap(S.stdX,S.cmpX,S.RcmpChs,[],[],[1],1.36,1000,68,[],1);
            fprintf('.')
        end
        % -------------------------------------------------------------------------------
        disp([ newline newline 'Ratio Correct' ]);
        str=fieldnames(XS);
        for i = 1:length(str)
            disp(['    ' str{i} '    ' num2str(mean( XS.(str{i}).Rcorrect))]);
        end

        figure(fnum); %NOTE NOT TESTED INDIVIDUALLY
        %legend(gca,'show')
        subplot(length(SUBJ),1,Z)

        if length(SUBJ)>1
            ylbl=({subj; '% Comparison Chosen'});
            if Z == length(SUBJ)
                xlbl=('Disparity (arcmin)');
            else
                xlbl=('');
            end
        else
            ylbl=('% Comparison Chosen');
            xlbl=('Disparity (arcmin)');
        end

        if Z == 1
            ttl={[prjCode '-' expID]; 'All Blocks Combined'};
            ttl=strrep(ttl,'_','-');
        else
            ttl='';
        end

        expPsyCurvePlotAll(fsLbl,fsTck,xtckSD,ytckSD,shape,colors,XS,fldNames,pcData,mFit,sFit,tFit,bFit,bSubPlot,xlbl,ylbl,ttl)
        if length(SUBJ)==1 && length(prjCodes)==1
            input([newline 'Press Return: '],'s')
        end
        %legend({subj},'Location','northwest')

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    end
    %if bPlot==1
    %    suptitle(superTitle);
    %    %SUBJJ=repelem(SUBJ,2);
    %    %lgg.String=SUBJJ;
    %end
end

if exist('files','var')~=1
    disp(newline)
    disp('No data available to calculate thresholds.');
    input('Press Return.');
end
