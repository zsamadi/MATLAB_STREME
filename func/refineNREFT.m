
function motifsRefined=refineNREFT(seedsNREF, seedsInfo, XnMerZoops,seqData, options)

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
%       - motifsRefined.thrOpt  : optimum thresholds



NREF=options.NREF;
prior=options.REFPrior;
nRefIter=options.nRefIter;




posSeeds=seedsInfo.seeds;
posPvalues=seedsInfo.pvalues;
PWM0=seqData.PWM0;
nPos=seqData.lens(1);
nNeg=seqData.lens(2);

motifsNREF=seedsNREF.seeds(1:NREF, :);
% PWMSThr=nodesNREF.thrOptimum;

W=size(motifsNREF, 2);

cellCharsLen=length(PWM0);

priorBg=prior*PWM0/(1+prior);
PWM0L=log2(PWM0);

cellCharsDouble=(1:cellCharsLen);

imotifSeq=1;

pvalueOptimum=ones(NREF,1);
thrOptVec=zeros(NREF,1);
PWMSCell=cell(NREF, 1);

% motifsNREF=motifsNREF(2:end, :);
consensusSeedMat=motifsNREF;

for motifSeq =motifsNREF.'
    logPvalue=-1;
    logPvalue0=0;
    refIter=1;



    % initialize PWM1 with motif

    PWM1=(motifSeq(:)==cellCharsDouble)/(1+prior)+priorBg;
    PWM1L=log2(PWM1);

    % compute PWM score and find approximate matches
    PWMS=(PWM1L-PWM0L).';

    thrOptimum=0;
    consensusSeed=motifSeq;
    % Continue refining while pvalue improves upto nRefIter of iterations

    while (refIter<nRefIter && logPvalue<logPvalue0)

        posPWMS=scoreWords(posSeeds, PWMS);
        appxSeeds=posSeeds(posPWMS>=0, :);
        appxPosPvalues=posPvalues(posPWMS>=0);

        posPWMS=posPWMS(posPWMS>=0);
        [posPWMS, iN]=sort(posPWMS, 'descend');

        appxSeeds=appxSeeds(iN, :);
        consensusSeed0=consensusSeed;
        consensusSeed=appxSeeds(1, :);
        appxPosPvalues=appxPosPvalues(iN);

        appxSeedsChar=char(appxSeeds+64);

        [PWMSTV,~, j]=unique(posPWMS);
        PWMSTV=PWMSTV(end:-1: 1);

        thrPositions = accumarray(j,1);
        thrPositions=cumsum(thrPositions(end:-1:1));
    
        refCntOptimum=motifCnt(XnMerZoops, appxSeeds,thrPositions,  true);
        zoopsCnt=refCntOptimum.zoopsCnt;
%         totalCnt=refCntOptimum.totalCnt-[0,0;refCntOptimum.totalCnt(1:end-1, :)];

%         zoopsCnt=zoopsCnt(iN, :);




%         PWMSTV=PSSMThr(PWMS,max(posPWMS), 1);

%         PWMSTV=sort(posPWMS(PWMSTV>=0),'descend');
%         [PWMSTV,~, j]=unique(posPWMS);
%         PWMSTV=PWMSTV(end:-1:1);

% 
%         seedCntsPos = accumarray(j,zoopsCnt(:,1));
%         seedCntsNeg = accumarray(j,zoopsCnt(:,2));
%         seedCntsThr=[seedCntsPos, seedCntsNeg];
        seedCntsThr=cumsum(zoopsCnt);



        % compute pvalues for the thresholdings and obtain optimum pvalue
        pvalues=computePvalue(seedCntsThr, [nPos, nNeg], 'fisher');
        logPvalues=log(pvalues);

        [logPvaluemin, idx]=min(logPvalues);

        logPvalue0=logPvalue;
        logPvalue=logPvaluemin;

        thrOptimum0=thrOptimum;

        thrOptimum=PWMSTV(idx);

        selIndexes=(posPWMS>=thrOptimum);
        appxSeeds=appxSeeds(selIndexes, :);
        appxSeedsLogPvalues=log(appxPosPvalues(selIndexes));
        % update motif with the sites above the threshold, only if pvalue
        % is improved
        if logPvalue<logPvalue0
            if sum(selIndexes)>1
% 
%                 posPWMScores=posPWMS(posPWMS>=thrOptimum);
%                 [~, idx]=sort(posPWMScores, 'descend');
%                 appxSeeds=appxSeeds(idx, :);
%                 appxSeedsLogPvals=appxSeedsLogPvals(idx);
%                 refCntOptimum=motifCnt(XnMerZoops, appxSeeds, false);
%                 zoopsCnt=refCntOptimum.zoopsCnt(:,1)-[0;refCntOptimum.zoopsCnt(1:end-1, 1)];
                refCntOptimum=motifCnt(XnMerZoops, appxSeeds,thrPositions,  false);
                zoopsCnt=refCntOptimum.zoopsCnt(:, 1);
                totalCnt=refCntOptimum.totalCnt(:,1);

                wgt_pos_count=-zoopsCnt./totalCnt.*appxSeedsLogPvalues;
                seedsToPWM=appxSeeds;
                numSelIndex=size(seedsToPWM,1);

                seedCellTi=seedsToPWM(:)==cellCharsDouble;
                PWM10=seedCellTi.*repmat(wgt_pos_count(:), W,1);
                PWM10=reshape(PWM10,numSelIndex,W, cellCharsLen);
                PWM1=squeeze(sum(PWM10, 1));
                PWM1=PWM1+prior*PWM0;
                PWM1=PWM1./sum(PWM1,2);
                PWM1L=log2(PWM1);
                PWMS=(PWM1L-PWM0L).';


            else
                seedsToPWM=appxSeeds;
                wgt_pos_count=-appxSeedsLogPvalues;
                PWM1=wgt_pos_count*(seedsToPWM(:)==cellCharsDouble);
                PWM1=PWM1+prior*PWM0;
                PWM1=PWM1./sum(PWM1,2);
                PWM1L=log2(PWM1);
                PWMS=(PWM1L-PWM0L).';
            end
        end



        refIter=refIter+1;
    end
    PWMSCell{imotifSeq}=PWMS;


    thrOptVec(imotifSeq)=thrOptimum0;
    pvalueOptimum(imotifSeq)=logPvalue0;
    consensusSeedMat(imotifSeq, :)=consensusSeed0;


    imotifSeq=imotifSeq+1;


end

motifsRefined.seeds=consensusSeedMat;
motifsRefined.pvalues=pvalueOptimum;
motifsRefined.PWMSCell=PWMSCell;
motifsRefined.thrOpt=thrOptVec;
