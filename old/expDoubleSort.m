function [DATA,DP,cmpXsrt,SS,E,signal,subj] = expDoubleSort(D,prjCodes,expIDs,SUBJ)

DATA    = zeros(length(D.idxStd),D.trlPerExp*length(SUBJ),length(prjCodes));
cmpXsrt = zeros(length(D.idxStd),D.trlPerExp*length(SUBJ),length(prjCodes));
stdXAll = zeros(length(prjCodes)*D.trlPerExp*length(SUBJ),length(D.idxStd));
cmpXAll = zeros(length(prjCodes)*D.trlPerExp*length(SUBJ),length(D.idxStd));
RAll    = zeros(length(prjCodes)*D.trlPerExp*length(SUBJ),length(D.idxStd));
stdXpss = zeros(length(prjCodes),length(D.idxStd),D.trlPerExp*length(SUBJ));
cmpXpss = zeros(length(prjCodes),length(D.idxStd),D.trlPerExp*length(SUBJ));
Rpss    = zeros(length(prjCodes),length(D.idxStd),D.trlPerExp*length(SUBJ));
SSS     =  cell(length(prjCodes),length(D.idxStd));

for p = 1:length(prjCodes)
    E=loadExpStruct(prjCodes{p},expIDs{p});
    if p == 1
        signal  = zeros(length(prjCodes),length(D.idxStd),E.ncmp);
    end
    for s = 1:length(D.idxStd)
        cmpIind=[];
        stdIind=[];
        R=[];
        cmpX=[];
        stdX=[];
        SS=struct();
        for Z=1:length(SUBJ)
            subj=SUBJ{Z};
            for b = 1:E.blkPerExp
                S=load([E.serverDataDir filesep subj filesep E.SubjData.(subj).stdX.(D.stdflds{D.idxStd(s)}).data{b}]);
                cmpIind=[cmpIind; S.cmpIind];
                stdIind=[stdIind; S.stdIind];
                R=[R; S.R==S.cmpIntrvl];
                cmpX=[cmpX; S.cmpX];
                stdX=[stdX; S.stdX];
                if isempty(fieldnames(SS))
                    SS=S;
                else
                    SS=structmerge(SS,S);
                end
            end
            if length(SUBJ)>1
                subj='ALL';
            end
            nPssInd=ones(size(stdIind)).*p;
            [DATA(s,(1:E.trlPerExp)+(Z-1)*E.trlPerExp,p), cmpXsrt(s,(1:E.trlPerExp)+(Z-1)*E.trlPerExp,p)] = psyTrialSort(stdIind,cmpIind,nPssInd,R,cmpX);

        end
        SS.nPssInd=ones(size(stdIind)).*p;
        SS.cmpChs=R;
        SS.cmpX=cmpX;
        SSS{p,s}=SS;
        signal(p,s,:)=S.cmpXunq-S.stdXunq;

        stdXAll((1:E.trlPerExp*length(SUBJ))*p,s)=stdX;
        cmpXAll((1:E.trlPerExp*length(SUBJ))*p,s)=cmpX;
        RAll((1:E.trlPerExp*length(SUBJ))*p,s)=R;
        stdXpss(p,s,:)=stdX;
        cmpXpss(p,s,:)=cmpX;
        Rpss(p,s,:)=R;
    end
end

%SWAP RESPONSES BETWEEN TRIALS IN CIRCLE, WORKS FOR MORE THAN 2 PASSES
if D.bSwap==1
	disp('SWAPPING')
	P=length(prjCodes);
	for s = 1:length(D.idxStd)
		K=0;
		for i = 1:size(Rpss,3)
			if K>P-1; K=0; end
			if K~=0
                Rpss(:,s,i)=circshift(Rpss(:,s,i)',K)';
                DATA(s,i,:)=circshift(DATA(s,i,:),K);
            end
			K=K+1;
		end
	end
end


%PERCENT CORRECT to D-prime
if D.muFrmDp==1
    PCi=zeros(length(D.idxStd),E.ncmp);
    DP=zeros(length(D.idxStd),E.ncmp);
    for s = 1:length(D.idxStd)
        i=psyPercentChosen(stdXAll(:,s),cmpXAll(:,s),RAll(:,s));
        PCi(s,:)=i(:,1);
        for c = 1:E.ncmp
            DP(s,c)=percentCorrect2dprime(PCi(s,c),2,1);
        end
    end
elseif D.muFrmDp==2
    PCi=zeros(length(prjCodes),length(D.idxStd),E.ncmp);
    DP=zeros(length(prjCodes),length(D.idxStd),E.ncmp);
    for p = 1:length(prjCodes)
    for s = 1:length(D.idxStd)
        i=psyPercentChosen(stdXpss(p,s,:),cmpXpss(p,s,:),Rpss(p,s,:));
        PCi(p,s,:)=i(:,1);
        for c = 1:E.ncmp
            DP(p,s,c)=percentCorrect2dprime(PCi(p,s,c),2,1);
        end
    end
    end
end

%MERGE STRUCT
SS=cell(length(D.idxStd),1);
for s = 1:length(D.idxStd)
for p = 1:(length(prjCodes))
    if isempty(SS{s})
        SS{s}=SSS{p,s};
    else
        SS{s}=structmerge(SS{s},SSS{p,s});
    end
end
end
