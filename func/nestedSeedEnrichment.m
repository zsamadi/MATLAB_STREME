function [seedEnriched, iterNumber]=nestedSeedEnrichment(PWMS, seedsWMax, seqDataSpecs, options)


refIter=1;
allPvalues=seedsWMax.pvalues;

% allCounts=seedsInfo.counts;

% posSites=seedsInfo.sites;
allSeeds=seedsWMax.seeds;

% posmaxWidth=seedsInfo.maxWidth;



PWM0=seqDataSpecs.back{1};
PWM0L=log2(PWM0);

nPos=seqDataSpecs.lens(1);
nNeg=seqDataSpecs.lens(2);

prior=options.REFPrior;
nRefIter=options.nRefIter;

thrOptimum=0;
pvalue=2;
[numCells,wMax]=size(PWMS);
cellCharsDouble=(1:numCells);

% [~, consensusSeed]=max(PWMS, [], 2);
% consensusSeed=consensusSeed.';
PWMS0=PWMS;

if seqDataSpecs.mkvOrder>0
    PWM1=2.^(PWMS);
else

    PWM1=2.^(PWMS+PWM0L);
end
PWM10=PWM1;


D15=numCells.^(wMax-1:-1:0);


% lenInfoGPU=gpuArray([nPos, nNeg]);



while (refIter<=nRefIter)

    posPWMS=scoreWords(allSeeds,PWMS, seqDataSpecs);
    positiveFlags=posPWMS>=0;


    appxPvalues=allPvalues(positiveFlags);
    possibleSeeds=allSeeds(positiveFlags, :);





    positivePWMS=posPWMS(positiveFlags);

    [appxPWMS, iN]=sort(positivePWMS, 'descend');

    appxSeeds=possibleSeeds(iN, :);

    appxPvalues=appxPvalues(iN);
    


    appxSeedIndex=(sum(appxSeeds.*D15,2));

    [zoopsCntAll,totalCntAll] =seedEnrichCount(seedsWMax,appxSeedIndex,nPos);

%     tic
% 
%     [pZoopsCntAll0GPU, pTotalCntAll0GPU]=motifCntGPU(pXNmerSeqd,pWMerSites, appxSeedIndex,thrPositions);
%     [nZoopsCntAll0GPU, nTotalCntAll0GPU]=motifCntGPU(nXNmerSeqd,nWMerSites, appxSeedIndex,thrPositions);
% 
% toc

% 
% 
%     pRefCntAll0N=[pZoopsCntAll0GPU, pTotalCntAll0GPU];
%     nRefCntAll0N=[nZoopsCntAll0GPU, nTotalCntAll0GPU];
%     zoopsCntAll0=[pRefCntAll0N(:,1), nRefCntAll0N(:,1)];
%     totalCntAll0=[pRefCntAll0N(:,2), nRefCntAll0N(:,2)];







    effectIndex=sum(zoopsCntAll, 2)>0;

    appxSeeds=appxSeeds(effectIndex, :);
    appxPvalues=appxPvalues(effectIndex);
    zoopsCntAll=zoopsCntAll(effectIndex, :);
    totalCntAll=totalCntAll(effectIndex, :);

    appxPWMS=appxPWMS(effectIndex);


    [PWMSTV,~, j]=unique(appxPWMS);

    % We could just do the counting at the thresholds, but later we need
    % these counts to compute wgtPosCount values

    seedCntsThr= [accumarray(j,zoopsCntAll(:,1)),accumarray(j,zoopsCntAll(:,2))] ;
    seedCntsThr=cumsum(seedCntsThr(end:-1:1, :), 1);






    % compute pvalues for the thresholdings and obtain optimum pvalue



     pvalues=computePvalue(seedCntsThr, [nPos, nNeg], options.bernoulli);



    [pvaluemin, idx]=min(pvalues);


    pvalue0=pvalue;
    pvalue=pvaluemin;

    thrOptimum0=thrOptimum;

    thrOptimum=PWMSTV(end-idx+1);

    selIndexes=(appxPWMS>=thrOptimum);
    sumSelIndexes=sum(selIndexes);
    appxSeeds=appxSeeds(1:sumSelIndexes, :);
    appxSeedsLogPvalues=appxPvalues(1:sumSelIndexes);
    appxSeedsLogPvalues=min(appxSeedsLogPvalues, 0);
    PWM1Out=PWM1;


    % update motif with the sites above the threshold, only if pvalue
    % is improved
    if (pvalue<pvalue0 && refIter<nRefIter)
        PWM10=PWM1;

        if sumSelIndexes>1

            zoopsCnt=zoopsCntAll(1:sumSelIndexes,1);
            totalCnt=totalCntAll(1:sumSelIndexes,1);


            totalCnt(totalCnt<1)=1;

            wgtPosCount=-zoopsCnt./totalCnt.*appxSeedsLogPvalues;
            seedsToPWM=appxSeeds;
            numSelIndex=size(seedsToPWM,1);

            seedCellTi=seedsToPWM(:)==cellCharsDouble;
            PWM1T=seedCellTi.*repmat(wgtPosCount(:), wMax,1);
            PWM1T=reshape(PWM1T,numSelIndex,wMax, numCells);
            PWM1=squeeze(sum(PWM1T, 1));
            PWM1=PWM1.'+prior*PWM0;
            PWM1=PWM1./sum(PWM1);

            if (options.isPal)
                PWM1=(PWM1+PWM1(:, end:-1:1))/2;
            end



            PWM1L=log2(PWM1);
            if options.mkvOrder>0
                PWMS=PWM1L;
            else
                PWMS=(PWM1L-PWM0L);
            end


        else
            seedsToPWM=appxSeeds;
            wgtPosCount=-appxSeedsLogPvalues;
            PWM1=wgtPosCount*(seedsToPWM(:)==cellCharsDouble);
            PWM1=PWM1.'+prior*PWM0;
            PWM1=PWM1./sum(PWM1);

            if (options.isPal)
                PWM1=(PWM1+PWM1(:, end:-1:1))/2;
            end
            
            PWM1L=log2(PWM1);
            if options.mkvOrder>0
                PWMS=PWM1L;
            else
                PWMS=(PWM1L-PWM0L);
            end
        end
    refIter0=refIter;
    refIter=refIter+1;

    else
        refIter0=refIter;
    
        refIter=nRefIter+1;

    end



end

if refIter0==nRefIter
    [~, consensusSeed]=max(PWM1);
    seedEnriched.seed=consensusSeed;
    iterNumber=refIter;
    seedEnriched.pvalue=pvalue;
    seedEnriched.PWMS=PWM1;
    seedEnriched.thresholdOptimum=thrOptimum;
    seedEnriched.PWMOut=PWM1;
else

    [~, consensusSeed]=max(PWM10);
    seedEnriched.seed=consensusSeed;
    iterNumber=refIter0;
    seedEnriched.pvalue=pvalue0;
    seedEnriched.PWMS=PWM10;
    seedEnriched.thresholdOptimum=thrOptimum0;
    seedEnriched.PWMOut=PWM1Out;
end

end
