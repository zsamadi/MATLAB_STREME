function [seqData, sitesErased, seqDataCopy]=eraseMotif(outMotif,  seqData)

% Erasing final motif obtained from enrichment step
% Input arguments:
%  - outMotif:
%       - outMotif.PWMSE        : Final refined PWM score matrix
%       - outMotif.consensusSeed: Final consensus Seed
%       - outMotif.scoreThr     : score threshold used for erasing
%       - outMotif.enrichPvalue : Enrichment pvalue of the output motif
%   - seqData: 
%       - seqData.pSeq  : primary sequence
%       - seqData.nSeq  : control sequence
%       - seqData.pHSeq : primary holdoput sequence
%       - seqData.nHSeq : control holdoput sequence
%       - seqData.PWM0  : backgroung model
%       - seqData.lens  : data lengths [size(posSeq90, 1),size(negSeq90, 1),size(pHoldSeq, 1),size(nHoldSeq, 1)];
% Output arguments:
%   - seqData: output seqData structure with input motif erased



% Filter the input sequence with threshold score, letters only erased if
% they're part of a motis and their score is greater than zero



filterSpecs.PWMS=outMotif.PWMSE;
filterSpecs.scoreThr=outMotif.scoreThr;
filterSpecs.back=seqData.back;
filterSpecs.mkvOrder=0;

filterSpecs.rvc=seqData.rvc;


[pSeqFiltered, pSeqErased]=seqFilterNew(seqData.pSeq, filterSpecs);
[nSeqFiltered, nSeqErased]=seqFilterNew(seqData.nSeq, filterSpecs);
[pHSeqFiltered,pHSeqErased]=seqFilterNew(seqData.pHSeq, filterSpecs);
[nHSeqFiltered,nHSeqErased]=seqFilterNew(seqData.nHSeq, filterSpecs);


seqDataCopy.pSeq=seqData.pSeq;
seqDataCopy.nSeq=seqData.nSeq;
seqDataCopy.pHSeq=seqData.pHSeq;
seqDataCopy.nHSeq=seqData.nHSeq;



seqData.pSeq=pSeqFiltered;
seqData.nSeq=nSeqFiltered;
seqData.pHSeq=pHSeqFiltered;
seqData.nHSeq=nHSeqFiltered;

sitesErased.pSeq=pSeqErased;
sitesErased.nSeq=nSeqErased;
sitesErased.pHSeq=pHSeqErased;
sitesErased.nHSeq=nHSeqErased;
