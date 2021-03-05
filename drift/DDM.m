N=5;
uB= 100;
lB=-100;
uBias=3;
lBias=2;
subtract=3;
DV=0;
lambda=1000000;
mu=.2;
sigma=1;
%sequential

tt=zeros(N,1);
for i = 1:N
    eT=0;
    eC=0;
    t=1;
    while true
        t=t+1;
        eT(t)=normrnd(mu,sigma,1);

        if eT(t)>0
            eT(t)=eT(t)*uBias;
        else
            eT(t)=eT(t)*lBias;
        end

        if lambda ==0
            eC(t)=eC(t-1)-subtract+eT(t);
        else
            eC(t)=eC(t-1)*2^(-t/lambda)+eT(t);
        end
        if eC(t)>=uB || eC(t)<=lB; break; end
    end
    tt(i)=t;
    plot(1:t,eC);hold on
    plot([0 t],[0 0],'k')
    ylim([lB, uB])
end
hold off

