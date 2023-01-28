function [zoopsCnt,totalCnt] =seedEnrichCount(seedsWMax,xSeedIndex,nPos)

% from initial NEVAL motifs, chooses NREF motifs by counting all matches
% and approximate matches in primary and control data to the motifs and computing their pvalue
% input arguments:
%   - XnMerZoops: Cell containing zoops WMers of primary and control
%   sequences
%   - appxSeeds: Nm Seeds of length W that shoud be counted inside the zoops sites 
%   - initial: true or false,  if true, only total counts are calculated,
%   if false, accumulated counts are calculated
%  Output arguments:
%   - seedCount: zoops count of the appxSeeds in the input WMers

seedsIndex=seedsWMax.seedsIndex;
seedsSites=seedsWMax.siteCell;
seedsCounts=seedsWMax.counts;

seedsCountsSum=sum(seedsCounts, 2);


seedIndexg=gpuArray(seedsIndex);
xSeedIndexg=gpuArray(xSeedIndex);
[~,seedStFlit] = ismember(xSeedIndexg,seedIndexg);
seedStFlit=gather(seedStFlit);
seedsCountsFLit=seedsCountsSum(seedStFlit);
sitesFlag=repelem((1:length(seedStFlit)).', seedsCountsFLit,1);

seedStIdT=cell2mat(seedsSites(seedStFlit));

seedStIdT=seedStIdT(:, 1);
[~, iSit]=unique(seedStIdT, 'stable');
sitesFlagZp=sitesFlag(iSit);
seedStIdTz=seedStIdT(iSit);

[sitesU, iz, jz]=uniqueSorted(sitesFlagZp);


jzP=jz(seedStIdTz<=nPos);
jzN=jz(seedStIdTz>nPos);



nuMers=length(iz);
zoopsCntNP=zeros(nuMers, 1);
zoopsCntNN=zeros(nuMers, 1);
if ~isempty(jzP)
    pCountst = accumarray(jzP,1);
    zoopsCntNP(1:length(pCountst))=pCountst;
end
if ~isempty(jzN)   
    nCountst = accumarray(jzN,1);
    zoopsCntNN(1:length(nCountst))=nCountst;
end




zoopsCntT=[zoopsCntNP, zoopsCntNN];
zoopsCnt=zeros(length(xSeedIndex), 2);
zoopsCnt(sitesU, :)=zoopsCntT;
totalCnt=seedsCounts(seedStFlit,:);









