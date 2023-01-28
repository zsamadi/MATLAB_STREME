function outMotif=scoreModelPssm(seedsRefined, seqData)

% Testing refined motifs on the holdout data 
% input arguments:
%   - seedsRefined: struct for sorted input WMers
%       - seedsRefined.seeds   :  Refined NREF consensus motifs
%       - seedsRefined.pvalues : significance obtained by optimum
%       thresholding
%       - seedsRefined.PWMSCell: PWMS motif to be used for enrichment test
%       - seedsRefined.thrOpt  : optimum thresholds
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
%       - outMotif.testLogPvalue   : significance obtained from hold-out
%       - outMotif.trainLogPvalue  : significance obtained from training
%       - outMotif.PWMSE           : out PWMS motif to be erased
%       - outMotif.scoreThr        : Score thresholds




posHoldSeq=seqData.pHSeq;
negHoldSeq=seqData.nHSeq;
PWMSCell=seedsRefined.PWMSCell;
thrOpt=seedsRefined.thrOpt;



nPos=size(posHoldSeq, 1);

nNeg=size(negHoldSeq,1);

if nPos>0

    
        PMSi=PWMSCell{1};

        % Compute Zoops WMers and their scores
        posWmerInfo=scoreWmer(posHoldSeq, PMSi);
        pZoopsSites=posWmerInfo.seq(:,end);
    
        negWmerInfo=scoreWmer(negHoldSeq, PMSi);
        nZoopsSites=negWmerInfo.seq(:,end);
        
        % Finds sites above threshold
        posSeqSites=pZoopsSites(posWmerInfo.score>=0);
        negSeqSites=nZoopsSites(negWmerInfo.score>=0);

        % only one best site from each sequence
        posCnt=length(unique(posSeqSites));
        negCnt=length(unique(negSeqSites));

        pEnrich=computePvalue([posCnt, negCnt], [nPos, nNeg], 'fisher');

else
    pEnrich=seedsRefined.pvalues(1);
end



outMotif.PWMSE=PWMSCell{1};
outMotif.consensusSeed=seedsRefined.seeds(1, :);
outMotif.scoreThr=thrOpt(1);
outMotif.trainLogPvalue=log(seedsRefined.pvalues(1));

outMotif.testLogPvalue=log(pEnrich);

