function outMotif=scoreModelPssm(seedsRefined, seqData, options)

% Testing refined motifs on the holdout data 
% input arguments:
%   - seedsRefined: struct for sorted input WMers
%       - seedsRefined.seeds   :  Refined NREF consensus motifs
%       - seedsRefined.pvalues : significance obtained by optimum
%       thresholding
%       - seedsRefined.PWMSCell: PWMS motif to be used for enrichment test
%       - seedsRefined.thresholdOptimum  : optimum thresholds
%   - seqData: 
%       - seqData.pSeq  : primary sequence
%       - seqData.nSeq  : control sequence
%       - seqData.pHSeq : primary holdoput sequence
%       - seqData.nHSeq : control holdoput sequence
%       - seqData.PWM0  : backgroung model
%       - seqData.lens  : data lengths [size(posSeq90, 1),size(negSeq90, 1),size(pHoldSeq, 1),size(nHoldSeq, 1)];
%  Output arguments:
%   - outMotif: struct for sorted input WMers
%       - outMotif.consensusSeed   :  output consensus motif
%       - outMotif.testPvalue   : significance obtained from hold-out
%       - outMotif.trainPvalue  : significance obtained from training
%       - outMotif.PWMSE           : out PWMS motif to be erased
%       - outMotif.scoreThr        : Score thresholds




posHoldSeq=seqData.pHSeq;
negHoldSeq=seqData.nHSeq;
PWMSCell=seedsRefined.PWMSCell;
PWMOutCell=seedsRefined.PWMOutCell;

thrOpt=seedsRefined.thresholdOptimum;

if options.rvp
    nPos=length(posHoldSeq)/2;
    nNeg=length(negHoldSeq)/2;
else
    
    nPos=length(posHoldSeq); 
    nNeg=length(negHoldSeq);
end

scrSpecs.back=seqData.back;
scrSpecs.mkvOrder=seqData.mkvOrder;
scrSpecs.rvp=options.rvp;


if nPos>0

    if scrSpecs.mkvOrder>0
        PMSi=log2(PWMSCell{1});
    else
        PMSi=log2(PWMSCell{1})-log2(scrSpecs.back{1});
    end        
        scoreThreshold=thrOpt(1);

        % Compute Zoops WMers and their scores
        posWmerInfo=scoreSeq(posHoldSeq, PMSi, scrSpecs);
        pZoopsSites=posWmerInfo.seq(:,end);
    
        negWmerInfo=scoreSeq(negHoldSeq, PMSi, scrSpecs);
        nZoopsSites=negWmerInfo.seq(:,end);
        
        % Finds sites above threshold
        posSeqSites=pZoopsSites(posWmerInfo.score>=scoreThreshold);
        negSeqSites=nZoopsSites(negWmerInfo.score>=scoreThreshold);

        % only one best site from each sequence
        posCnt=length(unique(posSeqSites));
        negCnt=length(unique(negSeqSites));

        PEnrich=computePvalue([posCnt, negCnt], [nPos, nNeg], options.bernoulli);

else
    PEnrich=seedsRefined.pvalues(1);
end



outMotif.PWMSE=log2(PWMSCell{1})-log2(scrSpecs.back{1});
outMotif.PWMOut=log2(PWMOutCell{1})-log2(scrSpecs.back{1});

outMotif.cSeed=seedsRefined.seeds(1, :);
outMotif.scoreThr=thrOpt(1);
outMotif.trainPvalue=seedsRefined.pvalues(1);

outMotif.testPvalue=PEnrich;
outMotif.secSeeds=seedsRefined.seeds(2:end, :);

