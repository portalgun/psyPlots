function [exitflag] = expDoublePassSetup(prjCode,expID,subj,Dtype,prjCode2,expID2)
    exitflag=0;
    % -------------------------------------------------------------------------------
    %IGNORE DUP PRJCODES

    disp(newline)
    if length(unique(expID2))<length(expID2) && length(expID2)==2
        while true
            disp(newline)
            response=input('Ignore duplicate prjCodes?: ','s');
            switch strYN(response)
                case 2
                    return
                case 1
                    break
                case 0
                    prjCode=prjCode2;
                    expID=expID2;
                    break
                otherwise
                    disp('Invalid Response');
            end
        end
    elseif length(expID2)>2
        prjCode=prjCode2;
        expID=expID2;
    end
    % -------------------------------------------------------------------------------
    % PARAMS

    if ~iscell(prjCode); prjCode={prjCode}; end
    if ~iscell(expID); expID={expID}; end
    if ~iscell(subj); subj={subj}; end

    D=cell(length(prjCode),1);
    prjcodes=cell(length(prjCode),1);
    expids=cell(length(prjCode),1);
    confAll=[];
    if length(subj)>1
        bMultSubj=1;
    else
        bMultSubj=0;
    end
    for W= 1:length(prjCode)
        [prjcodes{W,:},expids{W}]=expDoublePassList(prjCode{W},expID{W});
        if isempty(confAll) || confAll==0
            if isWnMngr
                [D{W},confAll,exitflag]=expDoublePassParams_window(prjCode{W},expID{W},subj{1},Dtype,confAll,bMultSubj);
            else
                [D{W},confAll,exitflag]=expDoublePassParams_headless(prjCode{W},expID{W},subj{1},Dtype,confAll);
            end
            if exitflag==1
                return
            end
        elseif confAll==1
            D{W}=D{W-1};
        end
    end
    % -------------------------------------------------------------------------------
    % DOUBLE PASS ANALYSIS

    if D{1}.bCombineSubj==0
        rho=cell(length(prjCode),length(subj));
        DP=cell(length(prjCode),length(subj));
        y=cell(length(prjCode),length(subj));
        signal=cell(length(prjCode),length(subj));
        Data=cell(length(prjCode),length(subj));
        yMenSub=cell(length(prjCode),length(subj));
        rhoMenSub=cell(length(prjCode),length(subj));
        bCombineCount=0;
        for W = 1:length(prjCode)
            for Z = 1:length(subj)
                [DP{W,Z},rho{W,Z},ncmp,sutitle,signal{W,Z},Data{W,Z},rhoMenSub{W,Z},yMenSub{W,Z},]=expDoublePass(D{W},prjcodes{W,:},expids{W,:},subj{Z},D{W}.bSave);
                disp(newline)
                bCombineCount=bCombineCount+D{W}.bCombineSubj;
            end
        end

        assignin('base','D',D)
        assignin('base','DP',DP)
        assignin('base','rho',rho)
        assignin('base','ncmp',ncmp)
        assignin('base','sutitle',sutitle)
        assignin('base','signal',signal)
        assignin('base','prjCode',prjCode)
        assignin('base','subj',subj)
        assignin('base','Data',Data)

        if length(subj)>1 && bCombineCount<length(subj) && length(D{1}.stdXo)>1
            %RATIO
            yMedTmp=zeros(length(prjCode),length(D{1}.stdXo),length(subj),ncmp);
            yMenTmp=zeros(length(prjCode),length(D{1}.stdXo),length(subj),ncmp);
            DPTmp=zeros(length(prjCode),length(D{1}.stdXo),length(subj),ncmp);
            signalTmp=zeros(length(prjCode),length(D{1}.stdXo),length(subj),ncmp);
            sI=zeros(length(prjCode),length(D{1}.stdXo),length(subj),ncmp);
            sE=zeros(length(prjCode),length(D{1}.stdXo),length(subj),ncmp);
            for W = 1:length(prjCode)
                for Z = 1:length(subj)
                    for s = 1:length(D{1}.stdXo)
                        for c = 1:ncmp
                            rhoTmp=abs(real(rho{W,Z}(s,c)));
                            rhoTmp2=abs(imag(rho{W,Z}(s,c)));
                            if rhoTmp<rhoTmp2
                                rhoTmp=rhoTmp2;
                            end
                            if rhoTmp==0
                                y{W,Z}(s,c)=10;
                            else
                                y{W,Z}(s,c)=sqrt((1-rhoTmp)./(rhoTmp));
                            end
                            yMedTmp(W,s,Z,c)=y{W,Z}(s,c);
                            yMenTmp(W,s,Z,c)=y{W,Z}(s,c);
                            DPTmp(W,s,Z,c)=mean(DP{W,Z}(:,s,c),1);
                            signalTmp(W,s,Z,c)=mean(DP{W,Z}(:,s,c),1);
                            sI(W,s,Z,c)=(signalTmp(W,s,Z,c)*(y{W,Z}(s,c)-1))/(y{W,Z}(s,c)*DPTmp(W,s,Z,c));
                            sE(W,s,Z,c)=(signalTmp(W,s,Z,c))/(y{W,Z}(s,c)*DPTmp(W,s,Z,c));
                        end
                    end
                end
            end

            yMedTmp=reshape(yMedTmp,length(prjCode),length(D{1}.stdXo),length(subj)*ncmp);
            yMenTmp=reshape(yMenTmp,length(prjCode),length(D{1}.stdXo),length(subj)*ncmp);

            %geometric mean of kappa accross subjects
            yMen=zeros(length(prjCode),length(D{1}.stdXo));
            yMed=zeros(length(prjCode),length(D{1}.stdXo));
            for W=1:length(prjCode)
                for s = 1:length(D{1}.stdXo)
                    yMed(W,s)=median(yMedTmp(W,s,:));
                    yMen(W,s)=geomean(yMenTmp(W,s,:));
                end
            end

            CIsz=68;
            numBoots=10000;

            means=zeros(length(prjCode),length(D{1}.stdXo),numBoots);
            for W=1:length(prjCode)
                for s = 1:length(D{1}.stdXo)
                for i = 1:numBoots
                    [~,indSmp]=datasample(squeeze(yMenTmp(W,s,:)),size(yMenTmp(W,s,:),3));
                    means(W,s,i)=geomean(squeeze(yMenTmp(W,s,indSmp)));
                end
                end
            end

            CIlohi = 0.5*(1-CIsz/100) + [0 CIsz/100];
            mCI=zeros(length(prjCode),length(D{1}.stdXo),2);
            for W=1:length(prjCode)
                for s = 1:length(D{1}.stdXo)
                    mCI(W,s,:) = quantile(means(W,s,:), CIlohi);
                end
            end

            %DP{W,Z}
            %rho{W,Z}
            %y{W,z}
            %ym{W,z}
            %signal{W,Z}
            %Data{W,Z}(s,c,:)

            %assignin('base','y',y)
            %assignin('base','yMed',yMed)
            %assignin('base','yMen',yMed)
        %save('2DsortedData','D','DP','rho','y','yM','ncmp','sutitle','signal','prjCode','subj','yMed','Data','yMen')
        %save('1DsortedData','D','DP','rho','y','yM','ncmp','sutitle','signal','prjCode','subj','yMed','Data','yMen')

        %PLOT COMBINED DATA PLOTS &&
            expDoublePassPopulation(signal,D,prjCode,subj,DP,rho,y,yMed,yMen,ncmp,sutitle,sI,sE,rhoMenSub,yMenSub,mCI)
        end

    elseif D{1}.bCombineSubj==1
        for W = 1:length(prjCode)
            expDoublePass(D{W},prjcodes{W},expids{W},subj,D{W}.bSave);
            disp(newline)
        end
    end
    return
end
