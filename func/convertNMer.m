function  [XNmerso, merEyson]=convertNMer(seq,seqLens, options)

% Sorts input sequence matrix in NMers and Zoops NMers
% Input arguments:
%   - seq   : input matrix of sequences
%   - W     : length of mers
%   - zpsOut: true or false, if true, outputs zoops nmenrs as well
% Output arguments:
%   - XZoops: Zoops WMers, 0 if zpsOut is false
%   - XNmers: output WMers
%used in countSeed and scoreSeq and seqFilterNew


numPSeqs=options.numPSeqs;
W=options.W;
rvcSeq=options.rvc;
allMers=options.allMers;

N=length(seqLens);
numShifts=seqLens-W+1;

seqLensC=cumsum(seqLens);
numShiftsRN=repelem([0;seqLensC(1:end-1)], numShifts);


totShifts=sum(numShifts);
numShiftsRi=(1:totShifts).';

numShiftsRit=[1;1+numShiftsRi(cumsum(numShifts(1:end-1)))];

numShiftsR=repelem(numShiftsRit, numShifts);

seqLensRi=numShiftsRi-numShiftsR;


tmp=repelem(seqLensRi.', W,1);
tmp=tmp(:);

tmp2=repelem((1:W), totShifts, 1);
tmp2=tmp2.';
tmp2=tmp2(:);

merEyso=(tmp2+tmp);


sitest=numShiftsRN;
sitesE=repelem(sitest, W, 1);
merEyson=merEyso+sitesE;


dSeqa=seq(merEyson);

XNmer=(reshape(dSeqa, W, totShifts)).';

if allMers
    XNmerErased=false(totShifts, 1);
else
    XNmerErased=any(XNmer==0, 2);
end
XNmer=XNmer(~XNmerErased, :);

if (rvcSeq)
    XNmerSitesP=repelem((1:numPSeqs).', numShifts(1:numPSeqs), 1);
    XNmerSitesP=XNmerSitesP(:);
    if options.pnSeq
        XNmerSitesN=repelem((numPSeqs+1:N/2).', numShifts(1:numPSeqs), 1);
        XNmerSitesN=XNmerSitesN(:);
        XNmerSites=[XNmerSitesP;XNmerSitesP;XNmerSitesN;XNmerSitesN];

    else
        XNmerSites=[XNmerSitesP;XNmerSitesP];
    end


else

    XNmerSites=repelem((1:N).', numShifts, 1);
    
    XNmerSites=XNmerSites(:);
end



XNmerMaxWidth=repelem(seqLens, numShifts)-seqLensRi;
XNmerMaxWidth=XNmerMaxWidth(~XNmerErased);



XNmerSites=XNmerSites(~XNmerErased);

XNmers=[XNmer, XNmerSites];
XNmerso=[XNmers, XNmerMaxWidth];


