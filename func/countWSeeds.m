function [unqMers, nMersAndSites]=countWSeeds(nMers, W,nMersMax, options)
numCells=options.numCells;
wMax=options.wMax;
numSeqs=options.numSeqs;

nMersW=[nMers(:, 1:W), nMers(:, end-1:end)];

nMersW(:, end)=-nMersW(:, end);
[~, iStW]=sortrows(nMersW );

nMers=nMers(iStW, :);

nMersAndSites=[nMers(:, 1:W), nMers(:, end-1)];


[nMerZoops, iZ] =unique(nMersAndSites, 'rows');


XNmerZoops=nMerZoops(:, 1:W);

% [unqEMers,iU, jU]=unique(XNmerZoops, 'rows');


[unqEMers,iU, jU]=uniqueSorted(XNmerZoops, numCells+1);
if nargout>1
    nMersAndSites=nMers(iZ, :);
    nMersAndSites=nMersAndSites(iU, :);
end

XNmerzMaxWidth=nMers(iZ, end);





jUP=jU(nMerZoops(:,end)<=numSeqs);
jUN=jU(nMerZoops(:,end)>numSeqs);
nuMers=size(unqEMers, 1);
pCounts=zeros(nuMers, 1);
nCounts=zeros(nuMers, 1);

if ~isempty(jUP)
pCountst = accumarray(jUP,1);
pCounts(1:length(pCountst))=pCountst;

end
if ~isempty(jUN)

nCountst = accumarray(jUN,1);
nCounts(1:length(nCountst))=nCountst;

end


uniqueMersV=(unqEMers(:, end)<numCells+1);

pCounts=pCounts(uniqueMersV);
nCounts=nCounts(uniqueMersV);
unqMers.counts=[pCounts, nCounts];


mWDTSort=sortrows([jU, XNmerzMaxWidth]);
XNmerzWidth=mWDTSort(:,2);
XNmerzMaxWidth=XNmerzWidth(iU);
XNmerzMaxWidth=XNmerzMaxWidth(uniqueMersV);
% STREME goes to this depth
% unqMers.maxWidth=min(XNmerzMaxWidth, wMax+1);
unqMers.maxWidth=XNmerzMaxWidth;


% siteCell is later used in initial and nested seed enrichment, and modify counts
sitesAll=nMerZoops(:, end);

iUD=[iU(2:end);size(nMerZoops,1)+1]-iU;
% this will be used in case of modify count
% sitesCell=mat2cell([sitesAll, XNmerzWidth],iUD);
sitesCell=mat2cell(sitesAll,iUD);

sitesCell=sitesCell(uniqueMersV);
unqMers.siteCell=sitesCell;
unqMers.w=W;



if W<wMax

    nMersMaxV=nMersMax(:, W)~=numCells+1;
    nMersMax=nMersMax(nMersMaxV, :);


    [nMersMaxT, iSt]=sortrows([nMersMax(:, 1:W),nMersMax(:, end-1), -nMersMax(:, end)]);

    [~,iT]=uniqueSorted(nMersMaxT(:, 1:W), numCells+1);
    nMersMaxSt=nMersMax(iSt, 1:end-1);

    nMersMaxEx=nMersMaxSt(:, 1:end-1);

    iTN=min(iT+1, [iT(2:end);length(iSt)]);

     unqMers.seeds=nMersMaxEx(iTN, :);
     unqMers.stExt=nMersMaxSt(iTN,end);

else
    unqEMers=unqEMers(uniqueMersV, :);
    unqMers.seeds=unqEMers;
    unqMers.seedsIndex=sum(unqEMers.*(numCells.^(W-1:-1:0)), 2);


end






