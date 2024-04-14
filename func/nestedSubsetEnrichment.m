
function motifsRefined=nestedSubsetEnrichment(seedsNREF, seedsInfo, seqDataSpecs, options)

% further refined NREF motifs by optimizing score threshold
% input arguments:
%   - seedsNREF : struct for sorted input WMers
%      - seedsNREF.seeds    :  NREF seeds to be further refined
%      - seedsNREF.pvalues  : significance obtained by fisher exact test
%   - seedsInfo:
%       - seedsInfo.seeds   : sorted seeds by pvalue
%       - seedsInfo.pvalues : sorted pvalues
%   - XnMerZoops: Cell containing zoops WMers of primary and control
%   - seqData:
%       - seqData.pSeq  : primary sequence
%       - seqData.nSeq  : control sequence
%       - seqData.pHSeq : primary holdoput sequence
%       - seqData.nHSeq : control holdoput sequence
%       - seqData.PWM0  : backgroung model
%       - seqData.lens  : data lengths [size(posSeq90, 1),size(negSeq90, 1),size(pHoldSeq, 1),size(nHoldSeq, 1)];
%   - options
%       - options.NEVAL     : Number of motifs for initial evaluation
%       - options.NREF      : Number of motifs for refinement
%       - options.prior     : Dirichlet prior for initial evaluation
%       - options.REFPrior  : Dirichlet prior for refinement step
%  Output arguments:
%   - motifsRefined: struct for sorted input WMers
%       - motifsRefined.seeds   :  Refined NREF consensus motifs
%       - motifsRefined.pvalues : significance obtained by optimum
%       thresholding
%       - motifsRefined.PWMSCell: PWMS motif to be used for enrichment test
%       on holdout sequences
%       - motifsRefined.thresholdOptimum  : optimum thresholds



NREF=options.NREF;
prior=options.REFPrior;


PWM0=seqDataSpecs.back{1};

% wSeeds=seedsNREF.wSeeds(1:NREF, :);
seedsNREF=seedsNREF.seeds(1:NREF, :);
% PWMSThr=nodesNREF.thresholdOptimumimum;




% numCells=length(PWM0);

priorBg=prior*PWM0/(1+prior);
PWM0L=log2(PWM0);



pvalueOptimum=ones(NREF,1);
thrOptVec=zeros(NREF,1);
iterNumbers=zeros(NREF,1);
PWMSCell=cell(NREF, 1);
PWMOutCell=cell(NREF, 1);

% motifsNREF=motifsNREF(2:end, :);
consensusSeedMat=seedsNREF;

for imotifSeq=1:NREF

    seedSeq =seedsNREF(imotifSeq, :);
    





    % initialize PWM1 with motif



    PWM1=(seedSeq(:)==(1:length(PWM0)))/(1+prior)+priorBg.';
    PWM1L=log2(PWM1.');

    if options.mkvOrder>0
        PWMS=PWM1L;
    else
        PWMS=(PWM1L-PWM0L);
    end



    % Continue refining while pvalue improves upto nRefIter of iterations

    [seedEnriched, iterNumber]=nestedSeedEnrichment(PWMS, seedsInfo, seqDataSpecs, options);


    PWMSCell{imotifSeq}=seedEnriched.PWMS;
    PWMOutCell{imotifSeq}=seedEnriched.PWMOut;


    thrOptVec(imotifSeq)=seedEnriched.thresholdOptimum;
    pvalueOptimum(imotifSeq)=seedEnriched.pvalue;
    consensusSeedMat(imotifSeq, :)=seedEnriched.seed;
    iterNumbers(imotifSeq)=iterNumber;


end

[pvalueOptimum,idx]=sort(pvalueOptimum, 'ascend');
PWMSCell=PWMSCell(idx);
thrOptVec=thrOptVec(idx);
seedsNREF=seedsNREF(idx, :);
consensusSeedMat=consensusSeedMat(idx, :);

motifsRefined.seeds=consensusSeedMat;
motifsRefined.pvalues=pvalueOptimum;
motifsRefined.PWMSCell=PWMSCell;
motifsRefined.thresholdOptimum=thrOptVec;
motifsRefined.PWMOutCell=PWMOutCell;


motifsRefined.seedsChar=[char(seedsNREF+64), num2str([pvalueOptimum, thrOptVec])];
motifsRefined.isPal=false;

if (options.rvp)
    
    [palModel, edOptimum]=getPalMode(PWMSCell{1}, PWM0);
    
    % Testing to see if palindromic model works better
    
    if seqDataSpecs.mkvOrder>0
        PWMSPal=log2(palModel);
    else
    
        PWMSPal=log2(palModel)-PWM0L;
    end
    
    options.isPal=true;
    
    motifPal=nestedSeedEnrichment(PWMSPal, seedsInfo, seqDataSpecs, options);
    
    logPvalueRatio=motifPal.pvalue/pvalueOptimum(1);
    
    if edOptimum<options.maxPalED && logPvalueRatio>=options.minPalRatio
        
        motifsRefined.seeds(1, :)=motifPal.seed;
        motifsRefined.pvalues(1)=motifPal.pvalue;
        motifsRefined.PWMSCell{1}=motifPal.PWMS;
        motifsRefined.PWMOutCell{1}=motifPal.PWMOut;
        motifsRefined.thresholdOptimum(1)=motifPal.thresholdOptimum;
        motifsRefined.isPal=true;
    end


end






% if iterNumbers(1)>=nRefIter
% 
%     motiFinal=nestedSeedEnrichment(PWMSCell(1), seedsInfo, wMerZoops,seqData, options);
% 
%     motifsRefined.seeds(1, :)=    motiFinal.seed;
%     motifsRefined.pvalues(1)=motiFinal.pvalue;
%     motifsRefined.PWMSCell(1)=motiFinal.PWMS;
%     motifsRefined.thresholdOptimum(1)=motiFinal.thresholdOptimum;
%     seedFinal=seedsNREF(1, :);
%     motiFinal.seedChar=[char(seedFinal+64), num2str([pvalueOptimum(1), thrOptVec(1)])];
%     motifsRefined.seedsChar(1, :)=motiFinal.seedChar;
% end




