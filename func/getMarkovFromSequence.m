function back=getMarkovFromSequence(seq,alphabet, order)

pseudo=1;

seqH=horzcat(seq{:});

alphSize=length(alphabet);

back0=sum(seqH(:)==(1:alphSize));

pseudoFrac=pseudo/alphSize;
back0=back0+pseudoFrac;

back0=back0(:)/sum(back0);

back=cell(order+1,1);
back{1}=back0;

numSeqs=length(seq);
seqLens=zeros(numSeqs, 1);

for iSq=1:numSeqs
    seqLens(iSq)=length(seq{iSq});
end


for orderL=1:order

    chainLen=alphSize^(orderL+1);

    pseudoFrac=pseudo/chainLen;

    kMers=convertNMer(seqH,seqLens, orderL+1, true);

    kMers=kMers(:, 1:orderL+1);
    
%     totalCount=size(kMers, 1);
    
    [kMersU, ~, jU]=unique(kMers, "rows");
    
    counts=accumarray(jU, 1);
    
    indexes=sum((kMersU-1).*(alphSize.^(orderL:-1:0)), 2)+1;
    
    countsAll=zeros(chainLen, 1);
    
    countsAll(indexes)=counts;
    
    countsAll=reshape(countsAll, alphSize, []);

    countsAll=countsAll+pseudoFrac;

    probsAll=countsAll./sum(countsAll);
    
    probsAll=probsAll(:);
    
    back{orderL+1}=probsAll.';
end







