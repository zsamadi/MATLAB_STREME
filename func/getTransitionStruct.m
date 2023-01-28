function gStruct=getTransitionStruct(network)


    endNodes=network.Edges.EndNodes;
    weights=network.Edges.Weight;
    numNodes=length(unique(endNodes));

    neighsCell = cell(numNodes,1);
    weightsCell= cell(numNodes,1);




    for i =1: numNodes
        neighsFlag= any(endNodes==i, 2);
        neighs =endNodes(neighsFlag, :);
        neighs=neighs.';
        neighs=neighs(:);
        neighs=neighs(neighs~=i);
        neighsCell{i}=neighs;
        weightsCell{i}=weights(neighsFlag);

    end

%     for i =1: numNodes
%         neighsWeights=weights(neighsFlagCell{i});
%         weightsCell{i}=neighsWeights/sum(neighsWeights);
%     end

    gStruct.neighs=neighsCell;
    gStruct.weights=weightsCell;