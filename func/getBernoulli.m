function    bernoulli=getBernoulli(posSeq, negSeq, W)
    numPSeqs=length(posSeq);
    numNSeqs=length(negSeq);
    pLens=zeros(numPSeqs, 1);
    nLens=zeros(numNSeqs, 1);
    
    for iP=1:numPSeqs
        pLens(iP)=length(posSeq{iP});
    end
    
    for iN=1:numNSeqs
        nLens(iN)=length(negSeq{iN});
    end
    
    
    Lp=mean(pLens);
    Lc=mean(nLens);
    
    if Lp~=Lc
        Sp=numPSeqs*(Lp-W+1);
        Sc=numNSeqs*(Lc-W+1);
        bernoulli=Sp/(Sp+Sc);
    else
        bernoulli=-1;
    end