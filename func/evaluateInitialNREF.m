
function seedsRefined=evaluateInitialNREF(seedsInfo, XnMerZoops,seqData, options)

% from initial NEVAL motifs, chooses NREF motifs by counting all matches
% and approximate matches in primary and control data to the motifs and computing their pvalue
% input arguments:
%   - seedsInfo
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
%   - seedsRefined: struct for sorted input WMers
%      - seedsRefined.seeds: sorted NEVAL seeds, NREF of which to be further refined
%      - seedsRefined.pvalues: significance obtained by fisher exact test


NEVAL=options.NEVAL;
prior=options.prior;



posSeeds=seedsInfo.seeds;
motifsNEVAL=posSeeds(1:NEVAL, :);


nPos=seqData.lens(1);
nNeg=seqData.lens(2);


PWM0=seqData.PWM0;
priorBg=prior*PWM0/(1+prior);
PWM0L=log2(PWM0);


pvalueMin=zeros(NEVAL, 1);
thrOptimum=zeros(NEVAL, 1);
refCntCell=cell(NEVAL, 1);
posPWMS0Cell=cell(NEVAL, 1);
refCntLen=zeros(NEVAL, 1);
for imt=1:NEVAL

    if imt==25
        check=1;
    end

    motifSeq=motifsNEVAL(imt, :);

    % initialize PWM1 with motif
    PWM1=(motifSeq(:)==(1:length(PWM0)))/(1+prior)+priorBg;
    PWM1L=log2(PWM1);
    PWMS=(PWM1L-PWM0L).';


    % compute PWM score and find approximate matches


    posPWMS=scoreWords(posSeeds, PWMS);

    % threshold vector to be used for finding optimum threshold
%     PWMSTV=PSSMThr(PWMS, max(posPWMS), 1);

    % Count pos and neg seeds above each threshold

    appxSeeds=posSeeds(posPWMS>=0, :);
    posPWMS=posPWMS(posPWMS>=0);
    [posPWMS, idxP]=sort(posPWMS, 'descend');
    
    appxSeeds=appxSeeds(idxP, :);
    [PWMSTV,~, j]=unique(posPWMS);
    PWMSTV=PWMSTV(end:-1: 1);

    thrPositions = accumarray(j,1);
    thrPositions=cumsum(thrPositions(end:-1:1));



    refCnt=motifCnt(XnMerZoops, appxSeeds, thrPositions, true);
    zoopsCnt=refCnt.zoopsCnt;


    posPWMS0Cell{imt}=PWMSTV;

    seedCntsThr=cumsum(zoopsCnt);


    refCntCell{imt}=seedCntsThr;
    refCntLen(imt)=size(seedCntsThr, 1);

    % find best threshold by minimizing pvale

end


refCntZoopsCnt = cell2mat(refCntCell);


posInfoGPU=gpuArray(refCntZoopsCnt);
lenInfoGPU=gpuArray([nPos, nNeg]);
fisherTestMotif=computePvalue(posInfoGPU, lenInfoGPU, 'fisher');
fisherTestMotif = gather(fisherTestMotif);

ind0=1;
for imotifSeq=1:NEVAL
    if imotifSeq==25
        check=1;
    end
    ind1=refCntLen(imotifSeq);
    
    [minPvalue, idx]=min(fisherTestMotif(ind0:ind0+ind1-1));
    pvalueMin(imotifSeq)=minPvalue;
    posPWMS0i=posPWMS0Cell{imotifSeq};
    thrOptimum(imotifSeq)=posPWMS0i(idx);

    ind0=ind0+ind1;


end


[pvalueNREFSort, sindex]=sort(pvalueMin);
motifsNEVALNREF=motifsNEVAL(sindex, :);
seedsRefined.seeds=motifsNEVALNREF;
seedsRefined.pvalues=pvalueNREFSort;
seedsRefined.thrOptimum=thrOptimum;
seedsRefined.seedsChar=[char(motifsNEVALNREF+64), num2str([log(pvalueNREFSort), thrOptimum])];

check=1;
