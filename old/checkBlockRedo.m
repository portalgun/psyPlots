function [E] = checkBlockRedo(prjCodes,expIDs,SUBJ,Dtype)

fsLbl=12;
fsTck=12;
xtckSD=0;
ytckSD=1;
errorbars=2;
color=[0,0,1];
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
if ~iscell(prjCodes)
    prjCodes={prjCodes};
end
if ~iscell(expIDs)
    expIDs={expIDs};
end
if ~iscell(SUBJ)
    SUBJ={SUBJ};
end

bStart=1;
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
            str=num2fldstr(stdX);
        end
        if strcmp(Dtype,'experiment')
            stdflds=fieldnames(E.SubjData.(subj).stdX);
            ESubj=E.SubjData.(subj);
        elseif strcmp(Dtype,'prelim')
            stdflds=fieldnames(E.SubjData.(subj).prelim.stdX);
            ESubj=E.SubjData.(subj).prelim;
        elseif strcmp(Dtype,'train')
            stdflds=fieldnames(E.SubjData.(subj).train.stdX);
            ESubj=E.SubjData.(subj).train;
        end

        blkPerExp=E.trlPerExp/E.trlPerBlk;

        %INIT FLAG FIELD IF DOESNT EXIST
        for i = 1:length(E.stdXunqAll)
            if ~isfield(ESubj.stdX.(stdflds{i}),'flag')
                ESubj.stdX.(stdflds{i}).flag=zeros(blkPerExp,1);
            end
        end

        %GET BLOCK INDICES
        count=0;
        blkidx=zeros(length(E.stdXunqAll),blkPerExp);
        for i = 1:length(E.stdXunqAll)
            for j = 1:blkPerExp
                if ~isempty(ESubj.stdX.(stdflds{i}).data{j}) && ESubj.stdX.(stdflds{i}).blk(j)==1
                    count=count+1;
                    blkidx(i,j)=count;
                end
            end
        end

        % -------------------------------------------------------------------------------
        if strcmp(Dtype,'experiment')
            fdir=strcat(E.serverDataDir,filesep,subj,filesep);
        elseif strcmp(Dtype,'prelim')
            fdir=strcat(E.serverDataDir,filesep,subj,filesep,'prelim',filesep);
        elseif strcmp(Dtype,'train')
            fdir=strcat(E.serverDataDir,filesep,subj,filesep,'train',filesep);
        end
        figure(fnum);
        if length(SUBJ) == 1
            SPC=ceil(sqrt(length(stdX)));
            SPR=floor(sqrt(length(stdX)));
        else
            SPC=length(stdX);
            SPR=length(SUBJ);
        end
        count=0;
        AtCI=cell(length(stdX));
        AtFit=cell(length(stdX));
        pcData=cell(length(stdX));
        Arange=cell(length(stdX));
        fprintf('Please wait. \nBootstrapping ')
        yl=zeros(length(E.stdXunqAll),2);
        for i = 1:length(E.stdXunqAll)
            if ismember(E.stdXunqAll(i),stdX)
                files=strcat(fdir,ESubj.stdX.(stdflds{i}).data);
            else
                continue
            end

            %'Learning' Threshold plots
            tCI=cell(length(files),1);
            tFit=zeros(length(files),1);
            ApcData=cell(length(files),1);
            %count=length(E.stdXunqAll)*(idxStd(i)-1);
            %range=count+1:count+length(files);
            if sum(blkidx(i,:))>0
                range=min(blkidx(i,:)):max(blkidx(i,:));
            end

            for j = 1:blkPerExp
                if blkidx(i,j)~=0
                %count=count+1;
                    S=load(files{j});
                    if isfield(S,'RcmpChosen')
                        [~,~,~,tCI{j},~,~,~,~]  = psyfitgengaussBootstrap(S.stdX, S.cmpX, S.RcmpChosen,[],[],1,1.36,2,100,68,[],0,1);
                        %fits
                        [~,~,~,tFit(j),pcData{j},~,~] = psyfitgengauss(S.stdX, S.cmpX, S.RcmpChosen, [],[],1,1.36,2,0);
                    else
                        tCI{j}=[-10,10];
                        tFit(j)=0;
                        pcData{j}=zeros(length(E.ncmp),1);
                        ESubj.stdX.(stdflds{i}).flag(j)=1;
                        %saveExpStruct(E);
                    end
                    fprintf('.')
                else
                    continue
                end
            end
            count=count+1;
            sX=E.stdXunqAll(i);
            subplot(SPR,SPC,count+(Z-1)*length(stdX));
            flagidx=find(ESubj.stdX.(stdflds{i}).flag==1);
            if i == 1 && length(SUBJ)>1 && i==ceil(length(E.stdXunqAll)/2) && length(SUBJ)>1
                [yl(i,:)]=expPsyThreshBlk(flagidx,range,shape,errorbars,fsLbl,fsTck,xtckSD,ytckSD,sX,S,pcData,tFit,tCI,subj,1,1);
            elseif i == 1 && length(SUBJ)>1
                [yl(i,:)]=expPsyThreshBlk(flagidx,range,shape,errorbars,fsLbl,fsTck,xtckSD,ytckSD,sX,S,pcData,tFit,tCI,subj,1,0);
            elseif i==ceil(length(E.stdXunqAll)/2) && length(SUBJ)>1
                [yl(i,:)]=expPsyThreshBlk(flagidx,range,shape,errorbars,fsLbl,fsTck,xtckSD,ytckSD,sX,S,pcData,tFit,tCI,[],0,1);
            elseif length(SUBJ)>1
                [yl(i,:)]=expPsyThreshBlk(flagidx,range,shape,errorbars,fsLbl,fsTck,xtckSD,ytckSD,sX,S,pcData,tFit,tCI,[],0,0);
            else
                [yl(i,:)]=expPsyThreshBlk(flagidx,range,shape,errorbars,fsLbl,fsTck,xtckSD,ytckSD,sX,S,pcData,tFit,tCI);
            end
            Arange{count}=range;
            AtCI{count}=tCI;
            AtFit{count}=tFit;
            ApcData{count}=pcData;
        end

        [~,Y,CIn,CIp]=expPsyThreshSep(prjCode,expID,subj,Dtype,0);

        count=0;
        for i = 1:length(E.stdXunqAll)
            if ismember(E.stdXunqAll(i),stdX)
                files=strcat(fdir,ESubj.stdX.(stdflds{i}).data);
            else
                continue
            end
            count=count+1;

            X=[1:blkPerExp]+((count-1)*blkPerExp);
            Yline=ones(1,blkPerExp).*Y(count);
            CInLine=ones(1,blkPerExp).*CIn(count);
            CIpLine=ones(1,blkPerExp).*CIp(count);

            subplot(SPR,SPC,count+(Z-1)*length(stdX)); hold on
            plot(X,Yline,'b','linewidth',2); hold on
            plotfillederror(X,CInLine,CIpLine,color,0);
            ylim([min(yl(:,1)),max(yl(:,2))]);
        end

        hold off
        % -------------------------------------------------------------------------------
        if ~length(prjCodes)>1 || ~length(SUBJ)>1 %% XXX AT SOME POINT MAKE IT SO MULTIPLE PLOTS CAN BE HANDLED IN HERE
        while true
            %TOGGLE FLAGS
            while true
                if errorflag==1 || errorflag==2
                    response=cellstr(response');
                    disp('Invalid response(s): ');
                    disp(['    ' response{idxInv}]);
                end
                disp([newline 'See plot. The blue line indicates threshold for all blocks. Red indicates marked for rerun. Flags are also set by checkTrlCompletion.']);
                response=input('Enter space seperated corresponding run numbers to toggle flag for rerun, or press Return: ','s');
                if strcmp(response,'')
                    return
                end

                idxInv=valChars(response,[' ';numV]);
                if ~isempty(idxInv)
                    errorflag=1;
                    continue
                end
                idxInv=valRange(response,1:nStdX*length(files));
                if ~isempty(idxInv)
                    errorflag=2;
                    continue
                end

                response=strsplit(response)';
                response(strcmp('',response))=[];
                break
            end
            response=str2double(response);
            for i=1:length(response)
                [r,c]=find(blkidx==response(i));
                %toggle
                if ESubj.stdX.(stdflds{r}).flag(c)==1
                    ESubj.stdX.(stdflds{r}).flag(c)=0;
                    E.SubjData.(subj).stdX.(stdflds{r}).flag(c)=0;
                elseif ESubj.stdX.(stdflds{r}).flag(c)==0
                    ESubj.stdX.(stdflds{r}).flag(c)=1;
                    E.SubjData.(subj).stdX.(stdflds{r}).flag(c)=1;
                end
            end
            saveExpStruct(E,prjCode,expID);

            %REPLOT
            figure(fnum);
            count=0;
            for i = 1:length(E.stdXunqAll)
                if ismember(E.stdXunqAll(i),stdX)
                    files=strcat(fdir,ESubj.stdX.(stdflds{i}).data);
                    count=count+1;
                else
                    continue
                end
                sX=E.stdXunqAll(i);
                flagidx=find(ESubj.stdX.(stdflds{i}).flag==1);
                subplot(SPR,SPC,count+(Z-1)*length(stdX));
                expPsyThreshBlk(flagidx,Arange{count},shape,errorbars,fsLbl,fsTck,xtckSD,ytckSD,sX,S,ApcData{count},AtFit{count},AtCI{count});
            end
            count=0;
            for i = 1:length(E.stdXunqAll)

                if ismember(E.stdXunqAll(i),stdX)
                    files=strcat(fdir,ESubj.stdX.(stdflds{i}).data);
                else
                    continue
                end
                count=count+1;

                X=[1:blkPerExp]+((count-1)*blkPerExp);
                Yline=ones(1,blkPerExp).*Y(count);
                CInLine=ones(1,blkPerExp).*CIn(count);
                CIpLine=ones(1,blkPerExp).*CIp(count);

                subplot(SPR,SPC,count+(Z-1)*length(stdX)); hold on
                plot(X,Yline,'b','linewidth',2); hold on
                plotfillederror(X,CInLine,CIpLine,color,0);
                ylim([min(yl(:,1)),max(yl(:,2))]);
            end
            hold off
        end
        end

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    end
    suptitle(superTitle);
end
