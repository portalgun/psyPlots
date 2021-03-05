function [fPr,tPr] = dave_HitFA(criterion,DV_L,DV_R,X)
    tP=sum(DV_R(X>=criterion));
    fP=sum(DV_L(X>=criterion));

    fN=sum(DV_R(X<criterion));
    tN=sum(DV_L(X<criterion));

    tPr=tP./(tP+fN);
    fPr=fP./(fP+tN);
end
