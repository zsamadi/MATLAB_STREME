function wMerInfo= scoreSeq(dSeq, PWMS, scrSpecs)
% computes scores of all Zoops Wmers in dseq
% input arguments:
%     - dSeq: dLen*sLen input sequence
%     - PWMS: 15*W  PWM scores

%  Output arguments:
%      -wMerInfo: 
%           - wMerInfo.seq  : Zoops Wmers
%           - wMerInfo.score: scores of all possible Zoops WMers

% used in scoreModelPssm

[numCells, W]=size(PWMS);
Ns=length(dSeq);
dSeqLens=zeros(Ns, 1);

for i=1:length(dSeq)
    dSeqLens(i)=length(dSeq{i});
end

dSeqLens=dSeqLens(dSeqLens>0);

dSeqSt=horzcat(dSeq{:});

if (scrSpecs.rvp)
    cnvOptions.numPSeqs=Ns/2;
else
    cnvOptions.numPSeqs=Ns;
end


cnvOptions.W=W;
cnvOptions.rvp=scrSpecs.rvp;

cnvOptions.allMers=false;
cnvOptions.pnSeq=false;


nMers=convertNMer(dSeqSt,dSeqLens, cnvOptions);






nMersAndSites=nMers(:, 1:end-1);
nMerZoops =unique(nMersAndSites, 'rows');


XNmerZoops=nMerZoops(:, 1:W);

% [unqEMers,iU, jU]=unique(XNmerZoops, 'rows');


[XZoopsU, ~, jx]=uniqueSorted(XNmerZoops, numCells+1);




score=scoreWords(XZoopsU, PWMS, scrSpecs);
score=score(jx);
wMerInfo.seq=nMerZoops;
wMerInfo.score=score(:);


