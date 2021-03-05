function [shape, color, titl]=DPA_shapeFunc(e,se)
if     e==1 && se==1 %full
    shape='o';
    color='w';
    titl='Full Cues';
elseif e==1 && se==2 %no
    shape='s';
    color='w';
    titl='No Cues';
elseif e==2 && se==1 %no mono
    shape='^';
    color='r';
    titl='No Mono BG';
elseif e==2 && se==2 %no bino bg
    shape='v';
    color='b';
    titl='No Bino BG';
elseif e==2 && se==3 %no bino
    shape='d';
    color='c';
    titl='No BG';
elseif e==2 && se==4 %no bino, black tex
    shape='d';
    color='w';
    titl ='No BG, No Texture';
else
    shape='o';
    color='w';
    titl='';
end
