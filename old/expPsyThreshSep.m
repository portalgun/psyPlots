function [X,Y,CIn,CIp] = expPsyThreshSep(prjCodes,expIDs,SUBJ,Dtype,bPlot)

if ~exist('bPlot','var')==1 || isempty(bPlot)
    bPlot=1;
end
%THRESHOLD TREND (error colors)
fsLbl=12;
fsTck=12;
xtckSD=2;
ytckSD=1;
errorbars=0;
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

bAvgSub=0;
if length(SUBJ)>1
    while true
        response=input('Would you like to average across subjects?: ','s');
        switch strYN(response)
            case 0
                break
            case 1
                bAvgSub=1;
                break
            case 2
                return
        end
    end
end

bAvgPrj=0;
if length(prjCodes)>1
    while true
        response=input('Would you like to average across experiments?: ','s');
        switch strYN(response)
            case 0
                break
            case 1
                bAvgPrj=1;
                break
            case 2
                return
        end
    end
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
        stdflds=fieldnames(ESubj.stdX);
        blkPerExp=E.trlPerExp/E.trlPerBlk;
        count=0;
        for i = 1:length(E.stdXunqAll)
            for j = 1:blkPerExp
                if ~isempty(ESubj.stdX.(stdflds{i}).data{j})
                    count=count+1;
                    files{count}=strcat(fdir,ESubj.stdX.(stdflds{i}).data{j});
                end
            end
        end
        % -------------------------------------------------------------------------------
        if exist('files','var')~=1
            continue
        end

        if bPlot == 1
            figure(fnum);
            lgg=legend(gca,'show');
        end
        % -------------------------------------------------------------------------------

        %Diclude trials
        stdXunqAll=zeros(length(files),1);
        %sort data based on standard into XS
        XS=struct();
        for i = 1:length(files)
            S=load(files{i});
            stdXunqAll(i)=S.stdXunq;
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
                    files{i}
                    if ESubj.stdX.(stdflds{r}).blk(c)~=0
                        ESubj.stdX.(stdflds{r}).flag(c)=1;
                        saveExpStruct(E);
                    end
                end
            end
        end

        %Get all fieldnames and sort
        stdXunqAll=unique(stdXunqAll,'stable');
        [stdXunqAll,idx]=sort(stdXunqAll);
        fldNames=fieldnames(XS);
        fldNames=fldNames(idx);

        %Get all standards in numerical type
        stdXAll=cell(length(fldNames),1);
        for i=1:length(fldNames)
            stdXAll{i}=strrep(fldNames{i},'_','.');
            stdXAll{i}=str2double(strrep(stdXAll{i},'n','-'));
        end
        stdXAll=cell2mat(stdXAll);
        nStds=length(stdXAll); %number of stds/fields

        mCI=cell(nStds,1);
        sCI=cell(nStds,1);
        tCI=cell(nStds,1);

        mFit=cell(nStds,1);
        sFit=cell(nStds,1);
        bFit=cell(nStds,1);

        tFit=zeros(nStds,1);
        pcData=cell(nStds,1);
        pcFit=cell(nStds,1);

        if bAvgSub==1
            if (Z==1 && bAvgPrj==0) || (Z==1 && W==1)
                stdXSub       = cell(nStds,1);
                cmpXSub       = cell(nStds,1);
                RcmpChosenSub = cell(nStds,1);
            end
            for i = 1:nStds
                stdXSub{i}       = [stdXSub{i};       XS.(fldNames{i}).stdX];
                cmpXSub{i}       = [cmpXSub{i};       XS.(fldNames{i}).cmpX];
                RcmpChosenSub{i} = [RcmpChosenSub{i}; XS.(fldNames{i}).RcmpChosen];
            end
            if bPlot==1 && ( (Z==length(SUBJ) && bAvgPrj==0) || (Z==length(SUBJ) && W==length(prjCodes)) )
                for i = 1:nStds
                    [mCI{i},sCI{i},~,tCI{i},~,~,~,~]  = psyfitgengaussBootstrap(stdXSub{i},cmpXSub{i},RcmpChosenSub{i},[],[],1,1.36,2,100,68,[],0,1);
                    [mFit{i},sFit{i},bFit{i},tFit(i),pcData{i},pcFit{i},~] = psyfitgengauss(stdXSub{i},cmpXSub{i},RcmpChosenSub{i},[],[],1,1.36,2,0);
                end
				XSS=XS;
				StCI=tCI;
				StFit=tFit;
				SpcData=pcData;
				if bPlot==1
					expPsyThreshSepPlot([],[],errorbars,fsLbl,fsTck,xtckSD,ytckSD,SUBJ,nStds,XSS,SpcData,StFit,StCI,fldNames,E.stdXunqAll);
				end

				X=zeros(length(StFit),1);
				for i=1:nStds
					X(i)=XSS.(fldNames{i}).stdX(1);
					CIn(i)=StCI{i}(1);
					CIp(i)=StCI{i}(2);
				end
				[X,idx]=sort(X);
				Y=StFit;
				Y=Y(idx);
				CIn=CIn(idx);
				CIp=CIp(idx);
            end
		elseif bAvgSub==1 && bAvgPrj==0
            if Z==1 && W==1
                stdXSub       = cell(nStds,length(SUBJ));
                cmpXSub       = cell(nStds,length(SUBJ));
                RcmpChosenSub = cell(nStds,length(SUBJ));
            end
            for i = 1:nStds
                stdXSub{i,Z}       = [stdXSub{i,Z};       XS.(fldNames{i}).stdX];
                cmpXSub{i,Z}       = [cmpXSub{i,Z};       XS.(fldNames{i}).cmpX];
                RcmpChosenSub{i} = [RcmpChosenSub{i}; XS.(fldNames{i}).RcmpChosen];
            end
            if bPlot==1 && W==length(prjCodes) && Z==length(prjCodes)
				for ZZ = 1:length(prjCodes)
					for i = 1:nStds
						[mCI{i},sCI{i},~,tCI{i},~,~,~,~]  = psyfitgengaussBootstrap(stdXSub{i},cmpXSub{i},RcmpChosenSub{i},[],[],1,1.36,2,100,68,[],0,1);
						[mFit{i,Z},sFit{i,ZZ},bFit{i,ZZ},tFit(i,ZZ),pcData{i,ZZ},pcFit{i,ZZ},~] = psyfitgengauss(stdXSub{i,ZZ},cmpXSub{i,ZZ},RcmpChosenSub{i,ZZ},[],[],1,1.36,2,0);
					end
					XSS=XS;
					StCI=tCI;
					StFit=tFit;
					SpcData=pcData;
					if bPlot==1
						expPsyThreshSepPlot([],[],errorbars,fsLbl,fsTck,xtckSD,ytckSD,SUBJ,nStds,XSS,SpcData{:,Z},StFit{:,Z},StCI{:,Z},fldNames,E.stdXunqAll);
					end

					X=zeros(length(StFit),1);
					for i=1:nStds
						X(i)=XSS.(fldNames{i}).stdX(1);
						CIn(i)=StCI{i}(1);
						CIp(i)=StCI{i}(2);
					end
					[X,idx]=sort(X);
					Y=StFit;
					Y=Y(idx);
					CIn=CIn(idx);
					CIp=CIp(idx);
				end
            end
        elseif bAvgSub==0 && bAvgPrj==0
            %Get confidence intervals and fits from bootstraps
            fprintf('Please wait. \nBootstrapping ')
            for i = 1:nStds
                %confidence intervals
                [mCI{i},sCI{i},~,tCI{i},~,~,~,~]  = psyfitgengaussBootstrap(XS.(fldNames{i}).stdX,XS.(fldNames{i}).cmpX,XS.(fldNames{i}).RcmpChosen,[],[],1,1.36,2,100,68,[],0,1);
                %fits
                [mFit{i},sFit{i},bFit{i},tFit(i),pcData{i},pcFit{i},~] = psyfitgengauss(XS.(fldNames{i}).stdX,XS.(fldNames{i}).cmpX,XS.(fldNames{i}).RcmpChosen,[],[],1,1.36,2,0);
            %                 psyfitgengaussBootstrap(S.stdX,S.cmpX,S.RcmpChs,[],[],[1],1.36,1000,68,[],1);
                fprintf('.')
            end

            XSS=XS;
            StCI=tCI;
            StFit=tFit;
            SpcData=pcData;
            % -------------------------------------------------------------------------------
            if bPlot==1
                StFit
                expPsyThreshSepPlot(colors(Z),shape,errorbars,fsLbl,fsTck,xtckSD,ytckSD,subj,nStds,XSS,SpcData,StFit,StCI,fldNames,E.stdXunqAll);
            end

            X=zeros(length(StFit),1);
            for i=1:nStds
                X(i)=XSS.(fldNames{i}).stdX(1);
                CIn(i)=StCI{i}(1);
                CIp(i)=StCI{i}(2);
            end
            [X,idx]=sort(X);
            Y=StFit;
            Y=Y(idx);
            CIn=CIn(idx);
            CIp=CIp(idx);
        end

%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    end
    if bPlot==1
		if bAvgSub==1
			SUBJJ={'Avg','Avg'};
		else
			SUBJJ=repelem(SUBJ,2);
		end
		suptitle(superTitle);
		lgg.String=SUBJJ;
    end
end

if exist('files','var')~=1
    disp(newline)
    disp('No data available to calculate thresholds.');
    input('Press Return.');
end
