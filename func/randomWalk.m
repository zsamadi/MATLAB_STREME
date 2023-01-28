
function genSeqCellType=randomWalk(G, gStruct, RWParams)

SubtypeVec=G.Nodes.label(:,1);

numPaths=RWParams.numPaths;
pathLength=RWParams.pathLength;
alpha=RWParams.alpha;

plotGraph=RWParams.plotGraph;
h=RWParams.figHandle;

    numNodes = length(SubtypeVec);

    genSeq = zeros(pathLength, numNodes, numPaths);

    parfor i =1:numPaths
        for j =1:numNodes
            indexList = generateSequence(j, gStruct, pathLength, alpha);

            genSeq(:,j, i)=indexList;

        end
    end
    genSeq=(reshape(genSeq,pathLength, [])).';



    genSeqCellType=SubtypeVec(genSeq);

    numSeqs=numNodes*numPaths;


    if plotGraph
        ncolor=rand(numSeqs, 3);
    
    % Highlight the Selected RW Seqs in the Graph Data
    for seqID = 1: numSeqs
        
         highlight(h,genSeq(:,seqID),'NodeColor',ncolor(seqID, :),'EdgeColor',ncolor(seqID, :))
    end
    end
    
end



function result=generateSequence(startIndex, gStruct, pathLength, alpha)
    result=zeros(pathLength, 1);
    result(1) = startIndex;
    current = startIndex;

    randPathLength=rand(pathLength, 1);

    for i =2: pathLength
        if randPathLength(i) < alpha
            nextIndex = startIndex;
        else

            indexes=gStruct.neighs{current};
            probs = gStruct.weights{current};
            
            nextIndex =randsample(length(probs),1,true,probs) ;

            nextIndex=indexes(nextIndex);
            
        end


        result(i)= nextIndex;
        current = nextIndex;
    end
end




