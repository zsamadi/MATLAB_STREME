
function motiFinal=nestedFinalEnrichment(motifsRefined, seedsInfo, XnMerZoops,seqData, options)

% further refined NREF motifs by optimizing score threshold
% input arguments:
%   - motifsRefined: struct for sorted input WMers
%       - motifsRefined.seeds   :  Refined NREF consensus motifs
%       - motifsRefined.pvalues : significance obtained by optimum
%       thresholding
%       - motifsRefined.PWMSCell: PWMS motif to be used for enrichment test
%       on holdout sequences
%       - motifsRefined.thresholdOptimum  : optimum thresholds
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
%   - motiFinal: struct for sorted input WMers
%       - motiFinal.seed   :  Refined consensus motif
%       - motiFinal.pvalue : significance obtained by optimum
%       thresholding
%       - motiFinal.PWMS: PWMS motif to be used for enrichment test
%       on holdout sequences
%       - motiFinal.thresholdOptimum  : optimum thresholds


% motifsRefined.seeds=consensusSeedMat(idx, :);
% motifsRefined.pvalues=pvalueOptimum;
% motifsRefined.PWMSCell=PWMSCell(idx);
% motifsRefined.thresholdOptimum=thrOptVec(idx);

PWMS=motifsRefined.PWMSCell{1};
PWMS0=PWMS;
consensusSeed=motifsRefined.seeds(1,:);
thrOpt=motifsRefined.thresholdOptimum(1);
prior=options.REFPrior;
nRefIter=options.nRefIter;




posSeeds=seedsInfo.seeds;
posPvalues=seedsInfo.pvalues;
PWM0=seqData.PWM0;
nPos=seqData.lens(1);
nNeg=seqData.lens(2);

seedRefined=motifsRefined.seeds(1, :);
% PWMSThr=nodesNREF.thresholdOptimumimum;

W=size(seedRefined, 2);

cellCharsLen=length(PWM0);

PWM0L=log2(PWM0);

cellCharsDouble=(1:cellCharsLen);


    logPvalue=motifsRefined.logPvalues(1);
    logPvalue0=0;
    refIter=1;


    % Continue refining while pvalue improves upto nRefIter of iterations

    while (refIter<nRefIter && logPvalue<logPvalue0)

        posPWMS=scoreWords(posSeeds, PWMS);
        appxSeeds=posSeeds(posPWMS>=0, :);
        appxPosPvalues=posPvalues(posPWMS>=0);

        posPWMS=posPWMS(posPWMS>=0);
        [posPWMS, iN]=sort(posPWMS, 'descend');

        appxSeeds=appxSeeds(iN, :);
        
        appxPosPvalues=appxPosPvalues(iN);

        appxSeedsChar=char(appxSeeds+64);

        [PWMSTV,~, j]=unique(posPWMS);
        PWMSTV=PWMSTV(end:-1: 1);

        thrPositions = accumarray(j,1);
        thrPositions=cumsum(thrPositions(end:-1:1));
    
        refCntOptimum=motifCnt(XnMerZoops, appxSeeds,thrPositions,  true);
        zoopsCnt=refCntOptimum.zoopsCnt;

        seedCntsThr=cumsum(zoopsCnt);



        % compute pvalues for the thresholdings and obtain optimum pvalue
        pvalues=computePvalue(seedCntsThr, [nPos, nNeg], 'fisher');
        logPvalues=log(pvalues);

        [logPvaluemin, idx]=min(logPvalues);

        logPvalue0=logPvalue;
        logPvalue=logPvaluemin;

        thrOpt0=thrOpt;

        thrOpt=PWMSTV(idx);

        selIndexes=(posPWMS>=thrOpt);
        appxSeeds=appxSeeds(selIndexes, :);
        appxSeedsLogPvalues=log(appxPosPvalues(selIndexes));
        
        % update motif with the sites above the threshold, only if pvalue
        % is improved
        if logPvalue<logPvalue0
            PWMS0=PWMS;
            if sum(selIndexes)>1

                refCntOptimum=motifCnt(XnMerZoops, appxSeeds,thrPositions,  false);
                zoopsCnt=refCntOptimum.zoopsCnt(:, 1);
                totalCnt=refCntOptimum.totalCnt(:,1);
                totalCnt(totalCnt<1)=1;

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
            [~, consensusSeed]=max(PWM1, [], 2);
            consensusSeed=consensusSeed.';
        end



        refIter=refIter+1;
    end


    thrOptVec=thrOpt0;
    logPvalueOptimum=logPvalue0;





motiFinal.seed=consensusSeed;
motiFinal.logPvalue=logPvalueOptimum;
motiFinal.PWMS=PWMS0;
motiFinal.thresholdOptimum=thrOptVec;

motiFinal.seedsChar=[char(seedRefined+64), num2str([logPvalueOptimum, thrOptVec])];


