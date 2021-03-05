function [prjCodes,expIDs] = expDoublePassList(prjCode,expID)

E=loadExpStruct(prjCode,expID);
if isfield(E,'link') && isfield(E.link,'pass') && ~isempty(E.link.pass)
    passList=E.link.pass;
    passList{end+1}=[prjCode ' ' expID];
    passList=sortrows(passList,2);
else
    disp('The selected experiment is not linked to any other passes')
    input('Press Return: ','s')
    return
end

%DISPLAY PASSES
disp(newline)
for i = 1:length(passList)
    disp([num2str(i) '. ' passList{i}]);
end

%SELECT PASSES
while true
    response=input('Enter spaced values corresponding to the passes you want to analyze: ','s');
    if strcmp(response,'')
        return
    end
    response=strsplit(response,' ');
    response(strcmp('',response))=[];
    idxInvalid=valChars(response,numV);
    if ~isempty(idxInvalid)
        disp(['Invalid response ' response(idxInvalid)])
        continue
    end
    response=str2double(response);
    if length(response)<2
        disp('This query requires two or more values.')
        continue
    end
    break
end
[prjCodes,expIDs]=cellfun(@(s) strtok(s,' '),passList(response), 'UniformOutput',false);
for i = 1:length(expIDs)
    expIDs{i}=strrep(expIDs{i},' ','');
end
