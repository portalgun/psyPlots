function fig=Exp_plot_histogramsBootBatch(obj,opts,figNum)
%BOOTSTRAP ALL DATA
M  =zeros(size(obj.subsAllLimit,1),1);
SEM=zeros(size(obj.subsAllLimit,1),1);
CI =zeros(size(obj.subsAllLimit,1),2);
opts.shape=cell(size(obj.subsAllLimit,1),1);
opts.color=cell(size(obj.subsAllLimit,1),1);
opts.lColor=cell(size(obj.subsAllLimit,1),1);
opts.titl=cell(size(obj.subsAllLimit,1),1);
X=zeros(size(obj.subsAllLimit,1),1);
opts.Axis='square';

for i = 1:size(obj.subsAllLimit,1)
    sub=obj.subsAllLimit(i,:);
    x=[sub(2) sub(3)];
    switch num2strSane(x)
          case '1,1'
              X(i)=1;
          case '1,2'
              X(i)=6;
          case '2,1'
              X(i)=2;
          case '2,2'
              X(i)=3;
          case '2,3'
              X(i)=4;
          case '2,4'
              X(i)=5;
          otherwise
    end

    SS=obj.select(sub);
    data=SS.Data.(opts.field);

    [~,opts.lColor{i}]=def_Exp_subj(SS.Subjs);
    [opts.shape{i}, opts.color{i}, opts.titl{i}]=def_Exp_plotParamsByExp(obj.shapeFunc,sub);
    [M(i),SEM(i),CI(i,:)]=bootMeanBasic(data,opts.nBoot,opts.prcntChoose);
end
ind=cellfun(@(x) isempty(x),opts.titl);
opts.titl(ind)=[];
opts.titl=unique(opts.titl);
opts.titl=join(opts.titl,', ');
opts.titl=opts.titl{1};

fig=figure(figNum);
plot_histogramsBoot(M,SEM,CI,opts,figNum,X);

opts.figNum=figNum;
gFigure(opts);
save_plot(fig,opts);
