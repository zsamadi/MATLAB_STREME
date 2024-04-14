function [seedsEvaled, seedsWmax]=evaluateInitialSeeds(unqMersCell,options)

% sorts WMers of length W  by pvalues
% input arguments:
%   - infoMat: a matrix containing 
%     - uniWMers: all seeds of length W
%     - posCnts: counts of Wmer seeds in primary data
%     - negCnts: counts of Wmer seeds in control data
%  Output arguments:
%   - seedsEvaled: struct for sorted input WMers
%      - seedsEvaled.seeds: sorted WMers NEVAL of which to be passed to
%      refinement step
%      - seedsEvaled.pvalues: significance obtained by fisher exact test
NEVAL=options.NEVAL;
Nt=options.numSeqs;
numCells=options.numCells;



nPos=Nt(1);
nNeg=Nt(2);

numWs=length(unqMersCell);

unqWMersCell=cell(numWs, 1);
countsCell=cell(numWs, 1);
% scaleCell=cell(numWs, 1);

wSeeds=zeros(numWs, 1);
maxWidthCell=cell(numWs, 1);

wMax=size(unqMersCell{end}.seeds, 2);
merLength=zeros(numWs, 1);

for iInfomat=1:numWs
    unqMers=unqMersCell{iInfomat};

    nuMers=size(unqMers.seeds,1);


    wSeeds(iInfomat,1)=unqMers.w;
    maxWidthCell{iInfomat}=unqMers.maxWidth;


    unqWMersCell{iInfomat}=unqMers.seeds;
    countsCell{iInfomat}=unqMers.counts;
%     scaleCell{iInfomat}=unqMers.countScale;

    merLength(iInfomat)=nuMers;


end





unqWMersAll=cell2mat(unqWMersCell);
wSeeds=repelem(wSeeds,merLength,1);
maxWidths=cell2mat(maxWidthCell);

countsMat=cell2mat(countsCell);
% scaleMat=cell2mat(scaleCell);


% [unqWMersAll, iu]=unique(unqWMersAll, 'stable','rows');
% 
% countsMat=countsMat(iu, :);
% wSeeds=wSeeds(iu);



% compute pvaluse using fisher exact test, computation performed on GPU

fisherTest=computePvalue(countsMat, [nPos, nNeg], options.bernoulli);


% fisherTest=round(fisherTest, 10);



pvaluesWMax=fisherTest(end-merLength(end)+1:end);

countsWMax=unqMersCell{end}.counts;
uniWMersWMax= unqMersCell{end}.seeds;
uniWMersMWidth= unqMersCell{end}.maxWidth;

[pvaluesWMax, idxWmax]=sortrows([pvaluesWMax,-uniWMersMWidth, -countsWMax(:,1)], 'ascend');

% [pvaluesWMax, idxWmax]=sortrows([pvaluesWMax,-countsWMax(:,1), uniWMersWMax], 'ascend');




seedsWmax.pvalues=pvaluesWMax(:, 1);
seedsWmax.seeds=uniWMersWMax(idxWmax, :);
seedsWmax.counts=countsWMax(idxWmax, :);
seedsWmax.maxWidth=uniWMersMWidth(idxWmax);

seedsWmax.siteCell=unqMersCell{end}.siteCell(idxWmax);
seedsWmax.seedsIndex=unqMersCell{end}.seedsIndex(idxWmax);


% sorting with significance test results

[fisherTestSort, fisherPerm] = sortrows([fisherTest, -wSeeds,-countsMat(:, 1)], 'ascend');

fisherTestSort=fisherTestSort(:,1);

unqWMersAll=unqWMersAll(fisherPerm, :);

wSeeds=wSeeds(fisherPerm);
maxWidths=maxWidths(fisherPerm);

%discarding rv symmetric seeds

if (options.rvp)

    unqWMersAllS=unqWMersAll;
    swFlag=unqWMersAllS(:, 1)>unqWMersAllS(:, end);
    unqWMersAllS(swFlag, :)=revCmp(unqWMersAllS(swFlag, :), numCells);
    
    [~, iSW]=unique(unqWMersAllS, "rows", 'stable');
    unqWMersAll=unqWMersAll(iSW, :);
    
    maxWidths=maxWidths(iSW);
    wSeeds=wSeeds(iSW);
    fisherTestSort=fisherTestSort(iSW);

end




seedsNEVAL=NEVAL;

if numWs>1
    
for wi=wMax-numWs+1:wMax
    idxs=(1:length(wSeeds));
    extFlag=unqWMersAll(:, end)<numCells+1;
    idxs=idxs(extFlag);
    wSeedsExt=wSeeds(extFlag);
    wiFlad=(wSeedsExt==wi);

    idxs=idxs(wiFlad);
    if ~isempty(idxs)
        seedsNEVALN=idxs(min(NEVAL, length(idxs)));
        if seedsNEVALN>seedsNEVAL
            seedsNEVAL=seedsNEVALN;
        end
    end
end



    
    unqWMersAll=unqWMersAll(1:seedsNEVAL, :);
    fisherTestSort=fisherTestSort(1:seedsNEVAL, :);
    wSeeds=wSeeds(1:seedsNEVAL);
    maxWidths=maxWidths(1:seedsNEVAL);

%     [unqWMersAll, iu]=unique(unqWMersAll, 'stable','rows');
% 
%     fisherTestSort=fisherTestSort(iu);
%     wSeeds=wSeeds(iu);







else
     unqWMersAll=unqWMersAll(1:NEVAL, :);
    fisherTestSort=fisherTestSort(1:NEVAL, :);
    wSeeds=wSeeds(1:NEVAL);
    maxWidths=maxWidths(1:NEVAL);


end
 



% sorted seeds by pvalue
seedsEvaled.seeds=unqWMersAll;

seedsEvaled.wSeeds=wSeeds;
seedsEvaled.pvalues=fisherTestSort;
seedsEvaled.maxWidth=maxWidths;


