function []= expPsyCurveInd(prjCodes,expIDs,SUBJ,Dtype)
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
                if errorflag==1
                    response=cellstr(response');
                    disp('Invalid response(s): ');
                    disp(['    ' response{idxInv}]);
                end
                response=input('Enter corresponding number for the standard you want to check: ','s');
                if strcmp(response,'')
                    return
                end
                %if strcmp(response,'a')
                %    response=num2str(1:length(E.stdXunqAll));
                %end
                tmp=strsplit(response);
                if size(tmp(~cellfun('isempty',tmp)))~=1
                    errorflag=1;
                    continue
                end

                idxInv=valChars(response,[' ';numV]);
                if ~isempty(idxInv)
                    errorflag=1;
                    continue
                end
                idxInv=valRange(response,1:length(E.stdXunqAll));
                if ~isempty(idxInv)
                    errorflag=1;
                    continue
                end

                response=strsplit(response)';
                response(strcmp('',response))=[];
                break
            end
            idxStd=str2double(response);
            stdX=E.stdXunqAll(idxStd);
            fldnames=fieldnames(ESubj.stdX);
            fld=fldnames{idxStd};
        end

        % -------------------------------------------------------------------------------
        %GET BLOCK INDICES & FILES
        count=0;
        files=cell(1);
        blkidx=zeros(blkPerExp,1);
        for j = 1:blkPerExp
            if ~isempty(ESubj.stdX.(fld).data{j}) && ESubj.stdX.(fld).blk(j)==1
                count=count+1;
                blkidx(j)=count;
                files{count,1}=strcat(fdir,ESubj.stdX.(fld).data{j});
            end
        end
        % -------------------------------------------------------------------------------
        %LOAD FILES AND PLACE IN STRUCT

        XS=struct();
        l=length(files);
        XS.RcmpChosen = cell(l,1);
        XS.stdX       = cell(l,1);
        XS.cmpX       = cell(l,1);
        XS.cmpXunq    = cell(l,1);

        for i = 1:l
            S=load(files{i});
            XS.RcmpChosen{i} = S.RcmpChosen;
            XS.stdX{i}       = S.stdX;
            XS.cmpX{i}       = S.cmpX;
            XS.cmpXunq{i}    = S.cmpXunq';
        end

        % -------------------------------------------------------------------------------
        %GET CONFIDENCE INTERVALS AND FITS FROM BOOTSTRAPS
        mCI    = cell(l,1);
        sCI    = cell(l,1);
        tCI    = cell(l,1);
        mFit   = cell(l,1);
        sFit   = cell(l,1);
        bFit   = cell(l,1);
        tFit   = zeros(l,1);
        pcData = cell(l,1);
        pcFit  = cell(l,1);

        fprintf('Please wait. \nBootstrapping ')
        for j = 1:l
            %confidence intervals
            [mCI{j},sCI{j},~,tCI{j},~,~,~,~]  = psyfitgengaussBootstrap(XS.stdX{j}, XS.cmpX{j}, XS.RcmpChosen{j},[],[],1,1.36,2,100,68,[],0,1);
            %fits
            [mFit{j},sFit{j},bFit{j},tFit(j),pcData{j},pcFit{j},~] = psyfitgengauss(XS.stdX{j}, XS.cmpX{j}, XS.RcmpChosen{j}, [],[],1,1.36,2,0);
            %                 psyfitgengaussBootstrap(S.stdX,S.cmpX,S.RcmpChs,[],[],[1],1.36,1000,68,[],1);
            fprintf('.')
        end

        figure(fnum);
        %legend(gca,'show')
        %length(XS)
        if length(SUBJ)>1
            ylbl={subj;'% Comparison Chosen'};
            I=(Z-1)*length(mFit);
            expPsyCurvePlotInd(fsLbl,fsTck,xtckSD,ytckSD,shape,colors,XS,pcData,mFit,sFit,tFit,bFit,length(SUBJ),length(mFit),I,Z,ylbl);
        else
            expPsyCurvePlotInd(fsLbl,fsTck,xtckSD,ytckSD,shape,colors,XS,pcData,mFit,sFit,tFit,bFit);
        end

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    end
    suptitle({superTitle; ['Standard at ' num2str(stdX)]});
end
