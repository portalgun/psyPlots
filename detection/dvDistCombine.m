classdef dvDistCombine < handle
properties
    names
    INname
    dvDists

    LLRa
    LLRb
    LLR_IN

    LLRaHist
    LLRbHist

    LLRaCtrs
    LLRbCtrs
end
properties(Hidden = true)
    n
    ninds
    bIN
end
methods
    function obj=dvDistCombine(varargin)
        varargin(cellfun(@isempty,varargin))=[];
        obj.n=length(varargin);
        obj.ninds=varargin{1}.n;
        obj.names=varargin{1}.names;
        n=length(varargin);
        if isfield(varargin{1}.OUT,'LLR_IN')
            obj.bIN=1;
        else
            obj.bIN=0;
        end
        for i = 1:n
            obj.LLRa{i}=[varargin{i}.OUT.LLRa];
            obj.LLRb{i}=[varargin{i}.OUT.LLRb];
            if obj.bIN
                obj.LLR_IN{i}=[varargin{i}.OUT.LLR_IN];
                obj.INname{i}=[varargin{i}.INname];
            end
        end

        obj.get_LLRs;
    end

    function obj=get_LLRs(obj)
        LLRa=cell(obj.ninds,1);
        LLRb=cell(obj.ninds,1);
        for j = 1:obj.ninds
            bNew=1;
            for ind = 1:obj.n
                llra=obj.LLRa{ind}{j};
                llrb=obj.LLRb{ind}{j};
                if bNew==1
                    bNew=0;
                    LLRa{j}=exp(llra);
                    LLRb{j}=exp(llrb);
                else
                    LLRa{j}=LLRa{j}+llra;
                    LLRb{j}=LLRb{j}+llrb;
                end
            end
        end
        obj.LLRa=LLRa;
        obj.LLRb=LLRb;

        obj.LLRaHist=cell(obj.ninds,1);
        obj.LLRbHist=cell(obj.ninds,1);
        obj.LLRaCtrs=cell(obj.ninds,1);
        obj.LLRbCtrs=cell(obj.ninds,1);
        for j = 1:obj.ninds
            LLRs=[obj.LLRa{j}; obj.LLRb{j}];
            ctrs=bin_widths_FD(LLRs);
            if numel(ctrs) > 1000
                [~,ctrs]=hist(LLRs,1000);
            end
            [A,obj.LLRaCtrs{j}]=hist(obj.LLRa{j},ctrs);
            [B,obj.LLRbCtrs{j}]=hist(obj.LLRb{j},ctrs);
            obj.LLRaHist{j}=A./sum(A(:));
            obj.LLRbHist{j}=B./sum(B(:));
        end

    end
    function DV=convert(obj)
        DV=dvDists();

        DV.OUT.LLRa=obj.LLRa;
        DV.OUT.LLRb=obj.LLRb;

        DV.LLRa=obj.LLRaHist;
        DV.LLRb=obj.LLRbHist;
        DV.LLRaCtrs=obj.LLRaCtrs;
        DV.LLRbCtrs=obj.LLRbCtrs;

        DV.names=obj.names;
        DV.n=length(DV.LLRa);
        if obj.bIN
            DV.OUT.LLR_IN=obj.LLR_IN;
            DV.INname=obj.INname;
        end

        DV.get_ROCs();
        DV.get_AUCs();
    end
end
end
