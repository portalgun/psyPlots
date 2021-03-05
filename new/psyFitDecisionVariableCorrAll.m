
% NORMALIZED MODEL MEANS
if size(mu1Fix,1) == size(stdX,1)
    %%
    Rideal = mu1Fix;
    stdXunq = unique(stdX(:))';
    [~,~,~,~,DPfit,DPdta]=psyfitgengaussAll(stdX,cmpX,Rideal,stdXunq,[],[1],1.36,0);
    clear mu1Fix;
    for s = 1:length(unique(abs(stdX)))
    disp('psyFitDecisionVariableCorrAll: WARNING! mu1Fix is klooged!!! BE CAREFUL!!!');
    mu1Fix(:,s) = DPfit(:,s).*sqrt(.43)./sqrt(2);
    end
end

% PARAMETER BOUNDS
[LB,UB] = psyFitDecisionVariableCorrAllParamBounds(modelType,mu1Fix);

% INITIALIZE PARAMETERS
param0  = psyFitDecisionVariableCorrAllInit(LB,UB);

[paramF,negLL] = fmincon(@(param) psyFitDecisionVariableCorrAllNegLL(param,modelType,stdX,cmpX,RcmpChsn,nPssInd,rhoFix,mu1Fix,mu2Fix,cr1Fix,cr2Fix),param0,[],[],[],[],LB,UB,[],opts);

% FINAL PARAMETER VALUES
[rho,mu1,mu2,cr1,cr2] = psyFitDecisionVariableCorrAllParamUnpack(paramF,modelType,rhoFix,mu1Fix,mu2Fix,cr1Fix,cr2Fix);


function [mu1CI,mu2CI,rhoCI,cr1CI,cr2CI,mu1,mu2,rho,cr1,cr2] = psyFitDecisionVariableCorrAllBootstrap(modelType,stdX,cmpX,RcmpChsn,nPssInd,rhoFix,mu1Fix,mu2Fix,cr1Fix,cr2Fix,CIsz,nBoot,bPLOT,bQUIET)

    parfor i = 1:nBoot
        % PROGRESS REPORT
        if bQUIET~=1
            progressreport(i,10,nBoot)
        end
        % RESAMPLE DATA W. REPLACEMENT
        [~,indSmp] = datasample(DATA,size(DATA,1));
        %%%%%%%%%%%%%%%%%%%%%
        % BOOTSTRAP FITS!!! %
        %%%%%%%%%%%%%%%%%%%%%
        [mu1(i,1),mu2(i,1),rho(i,1),cr1(i,1),cr2(i,1)] = psyFitDecisionVariableAllCorr(modelType,DATA(indSmp,:),rhoFix,mu1Fix,mu2Fix,cr1Fix,cr2Fix,0,1);
    end

    % CONFIDENCE INTERVAL: LO & HI BOUNDS
    CIlohi = 0.5*(1-CIsz/100) + [0 CIsz/100];

    % BOOSTRAPPED CONFIDENCE INTERVALS
    mu1CI = quantile(mu1, CIlohi);
    mu2CI = quantile(mu2, CIlohi);
    rhoCI = quantile(rho, CIlohi);
    cr1CI = quantile(cr1, CIlohi);
    cr2CI = quantile(cr2, CIlohi);

    % BOOSTRAPPED STD ERR OF STATISTIC
    mu1SE = std(mu1(~isnan(mu1)));
    mu2SE = std(mu2(~isnan(mu2)));
    rhoSE = std(rho(~isnan(rho)));
    cr1SE = std(cr1(~isnan(cr1)));
    cr2SE = std(cr2(~isnan(cr2)));
end

function param0  = psyFitDecisionVariableCorrAllInit(LB,UB)

    param0  = randInterval([LB; UB]).*1/8;
end

function [negLLall,Npp,Nnn,Npn,Nnp,Ppp,Pnn,Ppn,Pnp] = psyFitDecisionVariableCorrAllNegLL(param,modelType,stdX,cmpX,RcmpChsn,nPssInd,rhoFix,mu1Fix,mu2Fix,cr1Fix,cr2Fix,bPLOT)


    % SORT TRIALS (TO BE SURE THEY ARE SORTED %
    bSorted = psyDataTrialSortCheck(stdX,cmpX,nPssInd);
    if bSorted == 0
        disp(['psyFitDecisionVariableCorrAllNegLL: WARNING! data not sorted!!!']);
    end

    % UNIQUE STANDARD VALUES
    stdXunq = unique(stdX);
    nStdUnq = length(stdXunq);

    % UNPACK PARAMETERS (RETURNS MATRIX OF PARAMETERS GIVEN SIZE OF mu1Fix
    [rho,mu1,mu2,cr1,cr2]=psyFitDecisionVariableCorrAllParamUnpack(param,modelType,rhoFix,mu1Fix,mu2Fix,cr1Fix,cr2Fix);

    % LOOP OVER STANDARDS
    for s = 1:length(stdXunq)
        % UNIQUE COMPARISON VALUES FOR EACH STANDARD
        cmpXunq = unique(cmpX(stdX==stdXunq(s)));
        nCmpUnq = length(cmpXunq);

        % LOOP OVER COMPARISONS
        for c = 1:length(cmpXunq)
            % INPUT CHECKING
            if ~isempty(mu1Fix) && ( size(mu1Fix,1) ~= nCmpUnq || size(mu1Fix,2) ~= nStdUnq )
                error(['psyFitDecisionVariableCorrAllNegLL: WARNING! mu1Fix is an inappropriate size! Fix it!!! UPDATE CHECK 4 FUTURE MODELS']);
            end
            if ~isempty(mu2Fix) && ( size(mu2Fix,1) ~= nCmpUnq || size(mu2Fix,2) ~= nStdUnq )
                error(['psyFitDecisionVariableCorrAllNegLL: WARNING! mu2Fix is an inappropriate size! Fix it!!! UPDATE CHECK 4 FUTURE MODELS']);
            end

            % FIND EACH STD/CMP CONDITION FOR EACH PASS
            indCndPss1 = cmpX==cmpXunq(c) & stdX==stdXunq(s) & nPssInd==1;
            indCndPss2 = cmpX==cmpXunq(c) & stdX==stdXunq(s) & nPssInd==2;

            % PROBABILITY OF AGREEMENTS AND DISAGREEMENTS GIVEN THE MODEL
            warning off; % IN CASE mu1Fix is -Inf OR Inf
            [Ppp(c,s),Pnn(c,s),Ppn(c,s),Pnp(c,s)]=psyFitDecisionVariableCorrFunc(rho(c,s),mu1(c,s),mu2(c,s),cr1(c,s),cr2(c,s));
            warning on;

            % NUMBER OF EACH TYPE OF RESPONSE
            Npp(c,s) = sum(RcmpChsn(indCndPss1) == 1 & RcmpChsn(indCndPss2) == 1);
            Nnn(c,s) = sum(RcmpChsn(indCndPss1) == 0 & RcmpChsn(indCndPss2) == 0);
            Npn(c,s) = sum(RcmpChsn(indCndPss1) == 1 & RcmpChsn(indCndPss2) == 0);
            Nnp(c,s) = sum(RcmpChsn(indCndPss1) == 0 & RcmpChsn(indCndPss2) == 1);
        end
    end

    % TOTAL NUMBER OF RESPONSES IN EACH CONDITION
    Ntrl = Npp + Nnn + Npn + Nnp;

    % TOTAL NUMBER OF TRIALS
    NtrlAll = sum(Ntrl(:));

    % NEGATIVE LOG-LIKELIHOOD PER CONDITION
    negLL = -( Npp.*log(Ppp) + Nnn.*log(Pnn) + Npn.*log(Ppn) + Nnp.*log(Pnp) );

    % MEAN NEGATIVE LOG-LIKELIHOOD PER TRIAL
    negLLall = sum(negLL(~isnan(negLL(:))));

end

function [mFit,sFit,bFit,tFit,DPfit,DPdta,negLL] = psyfitgengaussAll(Xstd,Xcmp,RcmpChsn,mFix,sFix,bFix,DPcrt,nIntrvl,bPLOT,nBoot,CIsz,prcntUse,xLbl,yLbl,color,shape,bPLOTindi)
    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % FIT PSYCHOMETRIC FUNCTIONS %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i = 1:length(XstdUnq)
        % FIT EACH PSYCHOMETRIC
        indCnd  = Xstd==XstdUnq(i);
        mInd = intersect(i,1:length(mFix));
        sInd = intersect(i,1:length(sFix));
        bInd = intersect(i,1:length(bFix));
        [mFit(i),sFit(i),bFit(i),tFit(i),PCdta(:,i),PCfit(:,i),negLL(i)]  = psyfitgengauss(         Xstd(indCnd),Xcmp(indCnd),RcmpChsn(indCnd),mFix(mInd),sFix(sInd),bFix(bInd),DPcrt,nIntrvl,bPLOTindi,xLbl,[],color,shape);
        % BOOTSTRAP EACH PSYCHOMETRIC
        if nBoot > 1
        [~,sCI(i,:),bCI(i,:),tCI(i,:),~,sDstb(:,i),bDstb(:,i),Tdstb(:,i)] = psyfitgengaussBootstrap(Xstd(indCnd),Xcmp(indCnd),RcmpChsn(indCnd),mFix(mInd),sFix(sInd),bFix(bInd),DPcrt,nIntrvl,nBoot,CIsz,prcntUse,0);
        end
        DPdta(:,i)  = percentCorrect2dprime(PCdta(:,i),nIntrvl);
        DPfit(:,i)  = percentCorrect2dprime(PCfit(:,i),nIntrvl);
    end
    %%
end
