function [D,confAll] = expDoublePassParams(prjCode,expID,subj,Dtype,confAll)

%Get Stds
errorflag=0;

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

%PROMPT FOR CONFIDENCE INTERVALS

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
        CI=str2double(response);
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
        nBoot=str2double(response);
        break
    end
end



%FIX VALUES
prompt={['Enter parameters values to fix, leave blank otherwise.' newline newline num2str(char(961)) ':'] ,[num2str(char(956)) num2str(char(8321)) ':'], [num2str(char(956)) num2str(char(8322)) ':'], [num2str(char(946)) num2str(char(8321)) ':'], [num2str(char(946)) num2str(char(8322)) ':']};
dlg_title='Fitting Parameters';
defaultans={'','','','',''};
response=inputdlg(prompt,dlg_title,1,defaultans);
response=cellfun(@str2num,response,'UniformOutput',0);
for k = 1:numel(response)
    if isnan(response{k})
        response{k} = '';
    end
end
if isempty(response)
    return
end
D.rhoFix=response{1};
D.mu1Fix=response{2};
D.mu2Fix=response{3};
D.cr1Fix=response{4};
D.cr2Fix=response{5};
try
    mu1=num2str(D.mu1Fix);
end
try
    mu2=num2str(D.mu2Fix);
end

% SELECT MODEL TYPE
str=cell(0);
%str{1}='R';
str2{1}='R';
str{1}= ...
    num2str(char(961));                          %rho
if ~strcmp(mu1,mu2) || isempty(D.mu1Fix) || isempty(D.mu2Fix)
    str2{end+1}='RMd';
    str{end+1}= ...
    [ ...
    num2str(char(961)) ', ' ...                     %rho
    num2str(char(956)) num2str(char(8321)) ' = ' ... %mu1
    num2str(char(956)) num2str(char(8322)) ' from d'''];       %mu2
end
if ~strcmp(mu1,mu2) || isempty(D.mu1Fix) || isempty(D.mu2Fix)
    str2{end+1}='RMMd';
    str{end+1}= ...
    [ ...
    num2str(char(961)) ', ' ...                     %rho
    num2str(char(956)) num2str(char(8321)) ' & ' ... %mu1
    num2str(char(956)) num2str(char(8322)) ' from d'''];       %mu2
end
if isempty(D.mu1Fix) || isempty(D.mu2Fix)
    str2{end+1}='RM';
    str{end+1}= ...
    [ ...
    num2str(char(961)) ', ' ...                     %rho
    num2str(char(956)) num2str(char(8321)) ' = ' ... %mu1
    num2str(char(956)) num2str(char(8322)) ];   ... %mu2
    %num2str(char(961)) num2str(char(8321)) ' = ' num2str(char(961)) num2str(char(8322)) ' = ... , ' ...                     %rho = rho
end
if ~strcmp(mu1,mu2) || isempty(D.mu1Fix) || isempty(D.mu2Fix)
    str2{end+1}='RMM';
    str{end+1}= ...
    [ ...
    num2str(char(961)) ', ' ...                     %rho
    num2str(char(956)) num2str(char(8321)) ', ' ... %mu1
    num2str(char(956)) num2str(char(8322)) ];       %mu2
end
if isempty(D.cr1Fix) || isempty(D.cr2Fix)
    str2{end+1}='RMMCC';
    str{end+1}= ...
    [ ...
    num2str(char(961)) ', ' ...                     %rho
    num2str(char(956)) num2str(char(8321)) ', ' ... %mu1
    num2str(char(956)) num2str(char(8322)) ', ' ... %mu2
    num2str(char(946)) num2str(char(8321)) ', ' ... %beta1
    num2str(char(946)) num2str(char(8321)) ];       %beta2
end

[response,ok]=listdlg('PromptString','Select Parameters to fit','SelectionMode','single','ListString',str);
if isempty(response) || ok==0
    return
end

D.muFrmDp=0;
D.modelType=str2(response);
if strcmp(D.modelType,'RMd')
    D.muFrmDp=1;
    D.modelType='R';
elseif strcmp(D.modelType,'RMMd')
    D.muFrmDp=2;
    D.modelType='R';
end
D.trlPerExp=E.trlPerExp;

