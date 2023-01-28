function [PWMSCell, thrOptimumVec]=refineNREF(pSeq, nSeq, motifsNREF, PWM0, prior)

% Refines NREF Motifs

% input arguments:
%     - pSeq: primary sequence
%     - nSeq: control sequence
%     - motifsNREF: input NREF motifs
%     - PWM0: backgroung frequncies
%     - prior: Dirichlet prior
% % Output arguments:
%         -PWMSCell: Refined PWMS scores of NREF motifs
%      - thrOptimumVec: optimum thresholds to minimize pvalue

priorBg=prior*PWM0/(1+prior);
cellCharsLen=length(PWM0);
cellCharsDouble=(1:cellCharsLen);
pLen=size(pSeq, 1);
nLen=size(nSeq, 1);
NREF=size(motifsNREF, 1);
PWM0L=log(PWM0);


% initialize PWM1 with motif
PWMSCell=cell(NREF,1);
thrOptimumVec=zeros(NREF,1);
motifSeqi=1;
for motifSeq=motifsNREF.'
    PWM1=(motifSeq(:)==cellCharsDouble)/(1+prior)+priorBg;
    
    % compute PWM score and find approximate matches
    pvalueThrmin=0.5;
    pvalueThrmin0=1;


    while(pvalueThrmin<pvalueThrmin0)
        PWM1L=log(PWM1);
        PWMS=(PWM1L-PWM0L).';
        pvalueThrmin0=pvalueThrmin;


        pscore= scoreWmer(pSeq, PWMS);
        [bestScoresPseq,bestIndexPseq]=max(pscore, [], 2);

        
        nscore= scoreWmer(nSeq, PWMS);
        bestScoresNseq=max(nscore, [], 2);

%         [bestScoresPseqSort, bestScoresPseqSortIndex]=sort(bestScoresPseq, "descend");
%         bestScoresNseqSort=sort(bestScoresNseq, "descend");

        bestScoresPN=[bestScoresPseq;bestScoresNseq];
        bestScoresPN=bestScoresPN(bestScoresPN>=0);
        bestScoresPNUnique=sort(unique(bestScoresPN),"descend");
%         bestScoresPNUnique=[bestScoresPNUnique;0];

        iScore=1;
        pvalueThrVec=zeros(length(bestScoresPNUnique),1);
        pvalueThrmin=1;

        for scoreThr=bestScoresPNUnique.'
            nPos=sum(bestScoresPseq>=scoreThr);

            nNeg=sum(bestScoresNseq>=scoreThr);

            if (nPos>0) && (nNeg>0)
                pvalueThr=1-binocdf(nPos,pLen,nNeg/nLen);
            else
                pvalueThr=1;
            end


            pvalueThrVec(iScore)=pvalueThr;

            if pvalueThr<pvalueThrmin
                pvalueThrmin=pvalueThr;
                thrOptimum=scoreThr;
            end
            iScore=iScore+1;
        end
        selIndex=(bestScoresPseq>=thrOptimum);
        numSelIndex=sum(selIndex);
        pSeqSel=pSeq(selIndex, :);


        bestIndexPseqSelFirst=bestIndexPseq(selIndex);
        W=length(motifSeq);
        indexIncr=(1:W*numSelIndex).';

        bestIndexPseqSel=numSelIndex*repmat(bestIndexPseqSelFirst-1, W,1)+indexIncr;
        pSeqSelW=pSeqSel(bestIndexPseqSel);

        PWM10=(pSeqSelW(:)==cellCharsDouble)/(1+prior)+repmat(priorBg, numSelIndex*W, 1);
        %     PWM10=(pSeqSelW(:)==cellCharsDouble);

        PWM11=reshape(PWM10,numSelIndex,W, cellCharsLen);
        PWM1N=squeeze(sum(PWM11, 1));
        PWM1N=PWM1N./sum(PWM1N, 2);
        if norm(PWM1N-PWM1)<1e-5 || pvalueThrmin==0
            pvalueThrmin0=0;
        end

    end
    PWM1L=log(PWM1);
    PWMS=(PWM1L-PWM0L).';

    PWMSCell{motifSeqi}=PWMS;

    thrOptimumVec(motifSeqi)=thrOptimum;

    motifSeqi=motifSeqi+1;



end

