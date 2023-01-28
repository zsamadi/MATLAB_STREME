function scoreBack=getScoreBack(seeds,backSpecs)    

backgound=backSpecs.back;
order=backSpecs.mkvOrder;

[numSeeds, W]=size(seeds);
    
    
    scoreBack=zeros(numSeeds, W);
    
    scoreBack(:,1)=log2(backgound{1}(seeds(:,1)));
    
    
    
    for i=2:W
        li=max(1, i-order);
    
        nWi=i-li;
    
        wrdi=sum((seeds(:, li:i)-1).*(numCells.^(nWi:-1:0)), 2)+1;
        scoreBack(:, i)=log2(backgound{nWi+1}(wrdi));
    end