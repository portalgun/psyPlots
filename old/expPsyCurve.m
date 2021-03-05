function [] = expPsyCurve(prjCode,expID,subj,Dtype)

while true
    response=input('Would you like to plot individual psychometric curves for a each block (0), or would you like to pool across all blocks (1)?: ','s');
    switch strYN(response)
        case 0
            expPsyCurveInd(prjCode,expID,subj,Dtype);
            return
        case 1
            while true
                response=input('Would you like to plot psychometric curves on one plot (1) or each to their own subplot (0)?: ','s');
                switch strYN(response)
                    case 0
                        bSubPlot=1;
                        break
                    case 1
                        bSubPlot=0;
                        break
                    otherwise
                        disp('Invalid response')
                end
            end
            expPsyCurveAll(prjCode,expID,subj,Dtype,bSubPlot);
            return
        case 2
            return
        otherwise
            disp('Invalid response')
    end
end
