function val=def_plot(param,n)
    switch param
        case 'color'
            map=colormap(purplegreen);
            ind=linspace(1,size(map,1),n);
            val=map(ind,:);
        case 'shape'
            map=['<','>','^','v','d','s','o','h','*','x','+','.'];
            l=length(map);
            ind=repmat(1:l,1,ceil(n/l));
            ind=ind(1:n);
            val=map(ind);
        case 'xtitles'
            val=cell(n,1);
            for i = 1:n
                val{i}=['Condition ' num2str(i)];
            end
    end
