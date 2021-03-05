function [line,color]=def_Exp_subj(subj)
switch subj
       case 'DNW'
           line='-';
           color='k';
       case 'JDB'
           line='--';
           color='r';
       case 'THD'
           line=':';
           color='b';
       otherwise
           line='-';
           color='k';
end
