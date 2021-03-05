function[shape,color,title]=def_Exp_plotParamsByExp(exp,subs)
    switch exp
    case {'DPA','DPA_shapeFunc'}
            switch numel(subs)
            case 2
                e =subs(1);
                se=subs(2);

            otherwise
                e=subs(2);
                se=subs(3);
            end
            [shape,color,title]=DPA_shapeFunc(e,se);
    end
