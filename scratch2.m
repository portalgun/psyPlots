ym=cell(length(prjCode),length(D{1}.idxStd));
for W= 1:2
for Z = 1:3
shape=markerDefs(subj{Z});
for s = 1:length(D{1}.idxStd)
    if Z==1
        ym{W,s}=yM{W,Z}(s);
    else
        %yM{W,Z}(s)
        %ym(W,s)
        ym{W,s}=[yM{W,Z}(s);ym{W,s}];
    end
    Y=y{W,Z}(s,:).*DP{W,Z}(s,:)
    scatter(D{W}.stdXo(s)*ones(ncmp,1),Y,['k' shape],'MarkerFaceColor',[1 1 1]);
    hold on
end
end
end


    function [pA,pC] = psyFitDecisionVariableCorrRspAgree(symblType,bSmooth,bPlotManual)%,symType)

        %SMOOTHS CURVES AND EXTENDS CURVES TO 100% AGREEMENT
        MU1=-3:.01:obj.mu1(1,1);
        MU2=-3:.01:obj.mu2(1,1);
        for u=1:length(mu1(1,:))-1
            MU1=[MU1, obj.mu1(1,u):.01:obj.mu1(1,u+1)];
            MU2=[MU2, obj.mu2(1,u):.01:obj.mu2(1,u+1)];
        end
        MU1=[MU1, obj.mu1(1,end):.01:3];
        MU2=[MU2, obj.mu2(1,end):.01:3];

        % FIND PERCENT AGREEMENT AND PERCENT COMPARISON CHOSEN FOR EACH SET OF MU
        pA=zeros(length(MU1),1);
        pC=zeros(length(MU1),1);
        for i = 1:length(MU1)
            % XXX
            [Ppp,Pnn,Ppn,Pnp]=psyFitDecisionVariableCorrFunc(rho,MU1(i),MU2(i),obj.cr1(1),obj.cr2(1),0);
            pA(i)=Pnn+Ppp;
            pC(i)=Ppp+Ppn/2+Pnp/2;
            [~, index] = min(abs(pC-.5));
        end
    end
