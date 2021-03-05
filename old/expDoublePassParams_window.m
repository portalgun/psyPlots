function [D,confAll,exitflag] = expDoublePassParams_window(prjCode,expID,subj,Dtype,confAll,bMultSubj)
E=loadExpStruct(prjCode,expID);
exitflag=0;

% -------------------------------------------------------------------------------
Title = 'MatPROMPT: Double Pass Analysis Params';
Options.Resize = 'on';
Options.Interpreter = 'tex';
Options.CancelButton = 'on';
Options.ApplyButton = 'off';
Options.ButtonNames = {'Continue','Cancel'}; %<- default names, included here just for illustration
Option.Dim = 3; % Horizontal dimension in fields

D       = struct();
Prompt  = {};
Formats = {};
DefAns  = struct([]);

% -------------------------------------------------------------------------------


%SWAP
r=1;
 Prompt(r,:)        = {'Swap responses every   ','Swap', 'trials'};
Formats(r,:).type   = 'edit';
Formats(r,1).format = 'integer';
Formats(r,1).limits = [0 10];
Formats(r,1).size   = 50;

%CI
r=2;
 Prompt(r,:)        = {'Confidence Intervals      ','CI','%'};
Formats(r,:).type   = 'edit';
Formats(r,1).format = 'float';
Formats(r,1).limits = [0 99.999];
Formats(r,1).size   = 50;
%DefAns.CI           = 0;

%nBoot
r=3;
 Prompt(r,:)        = {'    Number of bootstraps','nBoot',[]};
Formats(r,:).type   = 'edit';
Formats(r,1).format = 'integer';
Formats(r,1).limits = [0 1000000000];
Formats(r,1).span   = [1 1];
Formats(r,1).size   = 50;
%DefAns.nBoot        = 0;

% -------------------------------------------------------------------------------

%rhofix
r=4;

str                 = [num2str(char(961)) ' :'];
 Prompt(r,:)         = {str,'rhoFix',[]};
Formats(r,1).type   = 'edit';
Formats(r,1).format = 'float';
Formats(r,1).limits = [-1 1];
Formats(r,1).size   = 50;
Formats(r,1).span   = [1 1];

%mu1fix
r=5;

str=[num2str(char(956)) num2str(char(8321)) ':'];
 Prompt(r,:)         = {str,'mu1Fix',[]};
Formats(r,1).type   = 'edit';
Formats(r,1).format = 'float';
Formats(r,1).limits = [-5 5];
Formats(r,1).size   = 50;
Formats(r,1).span   = [1 1];

%mu2fix
r=6;
str=[num2str(char(956)) num2str(char(8322)) ':'];
 Prompt(r,:)         = {str,'mu2Fix',[]};
Formats(r,1).type   = 'edit';
Formats(r,1).format = 'float';
Formats(r,1).limits = [-5 5];
Formats(r,1).size   = 50;
Formats(r,1).span   = [1 1];

%mu3fix
r=7;

str=[num2str(char(946)) num2str(char(8321)) ':'];
 Prompt(r,:)         = {str,'cr1Fix',[]};
Formats(r,1).type   = 'edit';
Formats(r,1).format = 'float';
Formats(r,1).limits = [-5 5];
Formats(r,1).size   = 50;
Formats(r,1).span   = [1 1];

%mu4fix
r=8;

str=[num2str(char(946)) num2str(char(8322)) ':'];
 Prompt(r,:)         = {str,'cr2Fix',[]};
Formats(r,1).type   = 'edit';
Formats(r,1).format = 'float';
Formats(r,1).limits = [-5 5];
Formats(r,1).size   = 50;
Formats(r,1).span   = [1 1];

%standards
r=9;
vals=cellstr(num2str(E.stdXunqAll));
for i = 1:length(vals)
    strg=strsplit(vals{i});
    if length(strg)==2
        vals(i)=strg(2);
    end
end
 Prompt(r,:)        = {'Standard','stdXo',[]};
Formats(r,1).type   = 'list';
Formats(r,1).style  = 'listbox';
Formats(r,1).format = 'text'; % Answer will give value shown in items, disable to get integer
Formats(r,1).items  = vals;
Formats(r,1).labelloc= 'topleft';
Formats(r,1).limits = [0 length(vals)]; % multi-select
Formats(r,1).size   = [100 150];
DefAns(1).stdXo = vals;

%DSaveDir
info=psyComputerInfo;
r=10;
 Prompt(r,:)        = {'Save figures Automatically','bSave',[]};
Formats(r,1).type   = 'check';
if strcmp(info.localHostName,'BONICN') || strcmp(info.localHostName,'BONICH')
    DefAns.bSave        = true;
else
    DefAns.bSave        = false;
end

r=11;
 Prompt(r,:)        = {'Figure Save Directory','FigureFolder',[]};
Formats(r,:).type   = 'edit';
Formats(r,1).format = 'dir';
Formats(r,1).span   = [1 1];
Formats(r,1).size   = 400;
if strcmp(info.localHostName,'BONICN') || strcmp(info.localHostName,'BONICH')
    if isdir(['/home/dambam/FWork/' E.expID])
        DefAns.FigureFolder = ['/home/dambam/FWork/' E.expID];
    else
        DefAns.FigureFolder = '/home/dambam/FWork/tmp';
    end
else
    DefAns.FigureFolder = pwd;
end

[D,Cancelled] = inputsdlg(Prompt,Title,Formats,DefAns,Options);
if Cancelled==1
    exitflag=1;
    return
end
% -------------------------------------------------------------------------------
try
    D.mu1=num2str(D.mu1Fix);
    mu1=D.mu1;
end
try
   D.mu2=num2str(D.mu2Fix);
   mu2=D.mu2;
end

if isempty(D.Swap)
    D.bSwap=0;
else
    D.bSwap=1;
end
if isempty(D.CI)
    D.bCI=0; %AUTO
else
    D.bCI=1; %AUTO
end

D.idxStd=find(ismember(vals,D.stdXo));
D.stdXo=str2double(D.stdXo);

% -------------------------------------------------------------------------------

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

Title = 'MatPROMPT: Double Pass model type';
Options.Resize = 'on';
Options.Interpreter  = 'tex';
Options.CancelButton = 'on';
Options.ApplyButton = 'off';
Options.ButtonNames  = {'Continue','Cancel'}; %<- default names, included here just for illustration
Prompt  = {};
Formats = {};
DefAns  = struct([]);
r=1;

 Prompt(r,:)        = {'Parameters to fit','mType',[]};
Formats(r,:).type   = 'list';
Formats(r,1).style  = 'listbox';
Formats(r,1).format = 'text'; % Answer will give value shown in items, disable to get integer
Formats(r,1).items  = str;
Formats(r,1).limits = [1 1]; % multi-select
Formats(r,1).size   = [150 250];
Formats(r,1).labelloc= 'topleft';
DefAns(1).mType = str(3);

r=2;
 Prompt(r,:)        = {'Apply all params to remaining projects','confAll',[]};
Formats(r,1).type   = 'check';
DefAns(1).confAll = false;

%bCombineSubj
if bMultSubj==1 && isempty(confAll)
    r=3;
     Prompt(r,:)        = {'Combine subject data across experiments','bCombineSubj',[]};
    Formats(r,1).type   = 'check';
    DefAns.bCombineSubj = false;
end

[D2,Cancelled] = inputsdlg(Prompt,Title,Formats,DefAns,Options);
confAll=D2.confAll;
if isfield(D2,'bCombineSubj')
    D.bCombineSubj=D2.bCombineSubj;
else
    D.bCombineSubj=0;
end

if Cancelled==1
    exitflag=1;
    return
end

response=find(ismember(str,D2.mType));

%[response,ok]=listdlg('PromptString','Select Parameters to fit','SelectionMode','single','ListString',str);

D.muFrmDp=0;
D.modelType=str2(response);
if strcmp(D.modelType,'RMd')
    D.muFrmDp=1;
    D.modelType='R';
elseif strcmp(D.modelType,'RMMd')
    D.muFrmDp=2;
    D.modelType='R';
end

% -------------------------------------------------------------------------------
D.trlPerExp=E.trlPerExp;

if strcmp(Dtype,'experiment')
    D.stdflds=fieldnames(E.SubjData.(subj).stdX);
elseif strcmp(Dtype,'prelim')
    D.stdflds=fieldnames(E.SubjData.(subj).prelim.stdX);
elseif strcmp(Dtype,'train')
    D.stdflds=fieldnames(E.SubjData.(subj).train.stdX);
end
