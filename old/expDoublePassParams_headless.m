function [D,confAll,exitflag] = expDoublePassParams_headless(prjCode,expID,subj,Dtype,confAll)
%XXX
exitflag=1;
input('Headless support not yet complete, if you are not running headless or don''t know what that means, something went wrong. Press Return.','s')
return

%XXX NEED TO HANDLE THESE
D.rhoFix
D.mu1Fix
D.mu2Fix
D.cr1Fix
D.muFrmDp
D.modelType
D.trlPerExp


%Get Stds
errorflag=0;
% -------------------------------------------------------------------------------
if length(subj)>1
    while true
        response=input('Would you like to combine data across subjects?: ','s');
        switch strYN(response)
            case 0
                D.bCombineSubj=0;
                break
            case 1
                D.bCombineSubj=1;
                break
            case 2
                return
            otherwise
                disp('Invalid Response')
        end
    end
    disp(newline)
else
    length(subj)
    D.bCombineSubj=0;
end


% -------------------------------------------------------------------------------
while true
    response=input(['Would you like to save figures automatically?' newline '(Don''t do this unless you know what you are doing.): '],'s');
    switch strYN(response)
        case 0
            D.bSave=0;
            break
        case 1
            D.bSave=1;
            break
        case 2
            return
        otherwise
            disp('Invalid Response')
    end
end


% -------------------------------------------------------------------------------
%XXX
while true
    if errorflag==1 || errorflag==2
        response=cellstr(response');
        disp('Invalid response(s): ');
        disp(['    ' response{idxInv}]);
    else
        E=loadExpStruct(prjCode,expID);
        nStdX=length(E.stdXunqAll);

        disp(newline);
        for i = 1:nStdX
            disp([ num2str(i) ' ' num2str(E.stdXunqAll(i)) ]);
        end
    end
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
end
D.idxStd=str2double(response);
D.stdXo=E.stdXunqAll(D.idxStd);
if strcmp(Dtype,'experiment')
    D.stdflds=fieldnames(E.SubjData.(subj).stdX);
elseif strcmp(Dtype,'prelim')
    D.stdflds=fieldnames(E.SubjData.(subj).prelim.stdX);
elseif strcmp(Dtype,'train')
    D.stdflds=fieldnames(E.SubjData.(subj).train.stdX);
end

% -------------------------------------------------------------------------------
disp(newline)
while true
	response=input('Swap every other trial between passes?: ','s');
	switch strYN(response)
		case 0
			D.bSwap=0;
			break
		case 1
			D.bSwap=1;
			break
		case 2
			return
		otherwise
			disp('Invalid response')
	end
end

if isempty(confAll) || confAll==-1 || confAll==0
	disp(newline)
	while true
		response=input('Plot Confidence intervals? WARNING, this can take considerably longer to compute: ','s');
		if strcmp(response,'')
			return
		end
		switch strYN(response)
			case 0
				D.bCI=0;
				break
			case 1
				D.bCI=1;
				break
			otherwise
				disp('Invalid Response.')
		end
	end
elseif confAll==2
	D.bCI=0;
elseif confAll==1
	D.bCI=1;
end

% -------------------------------------------------------------------------------
if isempty(confAll)
	while true
		response=input('Is this true for everything that you are plotting?: ','s');
		if strcmp(response,'')
			return
		end
		switch strYN(response)
			case 0
				confAll=0;
				break
			case 1
				if D.bCI==0
					confAll=2;
				elseif D.bCI==1
					confAll=1;
				end
				break
			otherwise
				disp('Invalid Response.')
		end
	end
end


% -------------------------------------------------------------------------------
if D.bCI == 1
    while true
        response=input('Enter percent confidence interval (0 < CI < 1): ','s');
        if strcmp(response,'')
            return
        end
        try
            str2double(response);
        catch
            disp('Invalid response.')
            continue
        end
        if str2double(response) <= 0 || str2double(response) >= 1
            disp('Invalid response. Value must be between 0 and 1.')
            continue
        end
        D.CI=str2double(response);
        break
    end
    while true
        response=input('How many times to resample?: ','s');
        if strcmp(response,'')
            return
        end
        try
            str2double(response);
        catch
            disp('Invalid response.')
            continue
        end
        D.nBoot=str2double(response);
        break
    end
end
% -------------------------------------------------------------------------------




% ===============================================================================
% -------------------------------------------------------------------------------
%PROMPT FOR CONFIDENCE INTERVALS
% -------------------------------------------------------------------------------
%HEADLESS

%load data

%figure(nFn)
%for i = 1:length(D.idxStd)
%    subplot(ceil(length(D.idxStd)/2),ceil(length(D.idxStd)/2),i)
%    plot(DATA{i,1},DATA{i,2},'ko','markerface','w','markersize',8);
%end
%input('','s')



%for i = 1:length(response)
%    EL=loadExpStruct(prjCodes{i},expIDs{i});
%    EL.SubjData.(subj).
%    DATA(:,i)=R
%end

%XXX
%disp(newline)
%fprintf([num2str(1) ' ' num2str(char(961)) '\n'])
%fprintf([num2str(2) ' ' num2str(char(956)) num2str(char(8321)) '\n']);
%fprintf([num2str(3) ' ' num2str(char(956)) num2str(char(8322)) '\n']);
%fprintf([num2str(4) ' ' num2str(char(946)) num2str(char(8321)) '\n']);
%fprintf([num2str(5) ' ' num2str(char(946)) num2str(char(8322)) '\n']);
%disp('Enter space seperated values to fix the above parameters. Enter ''.'' to not fix a specific value or input 1 ''.'' to not fix anything.');
%response=input(': ','s');
%if strcmp(response,'')
%    return
%elseif strcmp(response,'.')
%    D.rhoFix=[];
%    D.mu1Fix=[];
%    D.mu2Fix=[];
%    D.cr1Fix=[];
%    D.cr2Fix=[];
%end


%WINDOW
