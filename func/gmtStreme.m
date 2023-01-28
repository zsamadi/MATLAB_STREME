function [outMotifCell,textOut,commandText, background]= gmtStreme(seqData, options)

% main function for applying streme on the input data
% Input arguments:
%   - seqData:
%       - seqData.pSeq  : primary sequence
%       - seqData.nSeq  : control sequence
%       - seqData.pHSeq : primary holdoput sequence
%       - seqData.nHSeq : control holdoput sequence
%       - seqData.lens  : data lengths [size(posSeq, 1),size(negSeq, 1),size(pHoldSeq, 1),size(nHoldSeq, 1)];
%   - options
%       - options.NEVAL     : Number of motifs for initial evaluation
%       - options.NREF      : Number of motifs for refinement
%       - options.prior     : Dirichlet prior for initial evaluation
%       - options.REFPrior  : Dirichlet prior for refinement step
% Output arguments:
%   - outMotifCell: Cell of outMotifs
%       - outMotif.PWMSE        : Final refined PWM score matrix
%       - outMotif.consensusSeed: Final consensus Seed
%       - outMotif.scoreThr     : score threshold used for erasing
%       - outMotif.testPvalue   : significance obtained from hold-out
%       - outMotif.trainPvalue  : significance obtained from training
%       - textOut               : output text indicating the reason for finishing the algorithm 
%       - commandText               : the command that generated the output results





commandText= sprintf("mStreme(pSeq={input data}  ,nSeq=={shuffled data},pHSeq={positive holdout seq}, nHSeq={control holdout seq},...\n PWM0={bg frequency}, wMin=%d, wMax=%d,threshold=%10.10f, alphabet={input alphabet})",...
     options.wMin, options.wMax, options.threshold);



background=getMarkovFromSequence(vertcat(seqData.pSeq, seqData.pHSeq),options.alphabet, options.mkvOrder);


options.numCells=length(options.alphabet);

if options.rvp
    numPSeqs=length(seqData.pSeq)/2;
    numNSeqs=length(seqData.nSeq)/2;
    numPHSeqs=length(seqData.pHSeq)/2;
    numNHSeqs=length(seqData.nHSeq)/2;
else
    numPSeqs=length(seqData.pSeq);
    numNSeqs=length(seqData.nSeq);
    numPHSeqs=length(seqData.pHSeq);
    numNHSeqs=length(seqData.nHSeq);
end

seqData.lens=[numPSeqs,numNSeqs,numPHSeqs,numNHSeqs];
seqData.mkvOrder=options.mkvOrder;
seqData.back=background;
seqData.rvp=options.rvp;


wMin=options.wMin;
wMax=options.wMax;

isEraseOccured=true;
iErase=1;
outMotifCell=cell(100,wMax-wMin+1);
cSeeds=zeros(100, wMax);
trainPvalues=zeros(100, wMax-wMin+1);
testPvalues=zeros(100, wMax-wMin+1);


validThreshold= log(options.threshold);
numPateince=0;


coptions.numCells=options.numCells;
coptions.wMax=options.wMax;
coptions.MinSeedWidth=options.MinSeedWidth;

coptions.numSeqs=numPSeqs;
  
coptions.rvp=options.rvp;


eoptions.NEVAL=options.NEVAL;
eoptions.numCells=options.numCells;
eoptions.numSeqs=seqData.lens(1:2);
eoptions.rvp=options.rvp;
eoptions.bernoulli=options.bernoulli;


seqDataSpecs.back=background;
seqDataSpecs.lens=seqData.lens;
seqDataSpecs.mkvOrder=seqData.mkvOrder;

while(isEraseOccured && numPateince<options.patience)

    for iW=1:wMax-wMin+1

        % Finding first NEVAL significant seeds
        fprintf('Counting Seeds         (Motif #%d)... ', iErase)
        ticCS=tic;
        coptions.modCount=false;
        unqMersCell=countSeeds(seqData, coptions);
        elapsedTime=toc(ticCS);
        fprintf('Elapsed Time %3.3f seconds \n' , elapsedTime)

        fprintf('Evaluate Initial Seeds (Motif #%d)... ', iErase)


        ticEV=tic;

        
        [seedsNEVAL, seedsWmax]=evaluateInitialSeeds(unqMersCell, eoptions);
        elapsedTime=toc(ticEV);
        fprintf('Elapsed Time %3.3f seconds \n' , elapsedTime)


        % refining first NEVAL seeds to find first NREF seeds
%         tic
        fprintf('Refine Initial Seeds   (Motif #%d)... ', iErase)

        ticRI=tic;
        seedsNREF=refineInitialSeeds(seedsNEVAL,seedsWmax, seqDataSpecs, options);
        elapsedTime=toc(ticRI);
        fprintf('Elapsed Time %3.3f seconds \n' , elapsedTime)
%         toc



        % Further refine first NREF seeds
        % uses following data from seqData
        %  PWM0=seqData.PWM0;
        %  nPos=seqData.lens(1);
        %  nNeg=seqData.lens(2);
       fprintf('Nested Subset Enrich   (Motif #%d)... ', iErase)

        ticNS=tic;
        seedsRefined=nestedSubsetEnrichment(seedsNREF,seedsWmax,seqDataSpecs,options);
        elapsedTime=toc(ticNS);
        fprintf('Elapsed Time %3.3f seconds \n' , elapsedTime)


        % Enrichment test on first NREF motifs
        fprintf('Hold-out Scoring       (Motif #%d)... ', iErase)

        ticHS=tic;
        outMotif=scoreModelPssm(seedsRefined, seqData, options);
        elapsedTime=toc(ticHS);
        fprintf('Elapsed Time %3.3f seconds \n' , elapsedTime)

        fprintf('Erase Found Motif      (Motif #%d)... ', iErase)
         ticEM=tic;
        [seqData, sitesErased]=eraseMotif(outMotif,  seqData);
        elapsedTime=toc(ticEM);
        fprintf('Elapsed Time %3.3f seconds \n' , elapsedTime)






        outMotif.nsites=sum(sitesErased.pSeq)+sum(sitesErased.pHSeq);

        PWMSE=outMotif.PWMSE;
        if options.mkvOrder>0
           PWM1=2.^(PWMSE);
        else
            PWM1=2.^(PWMSE).*background{1};
        end

        PWM1=round(PWM1, 6);

        outMotif = rmfield(outMotif,'PWMSE');

        outMotif.PWM=PWM1;


        outMotifCell{iErase, iW}=outMotif;
        cSeeds(iErase, :)=outMotif.cSeed;
        trainPvalues(iErase, :)=outMotif.trainPvalue;
        testPvalues(iErase, :)=outMotif.testPvalue;
        
        % Check if erasure has been occured

        isEraseOccured=any(sitesErased.pSeq|sitesErased.nSeq);

%         if isEraseOccured
%         coptions.modCount=true;
% 
%             unqMersCellN=modifyCountSeeds(seqDataCopy, seqData,sitesErased,unqMersCell, coptions);
%         end


        if (options.evalue)
            validThreshold= log(options.threshold)-log(iErase);
        end

        if outMotif.testPvalue>validThreshold
            numPateince=numPateince+1;
        else
            numPateince=0;
        end

        if (options.nmotifs>0 && iErase==options.nmotifs)
            numPateince=options.patience +1;
        end


    end


    iErase=iErase+1;


end

if ~isEraseOccured
    textOut='Motif Analysis Finished Becasue No Site Erased\n';
else

    if iErase==options.nmotifs
           textOut=sprintf('Motif Analysis Finished Becasue %d Motifs Were Found\n', options.nmotifs);
    
    else
          textOut=sprintf('Motif Analysis Finished Becasue %d Consecutive Nonsignificant Motifs Were Found\n', options.patience);
    
    end


end
  fprintf(textOut);
mumFoundMotifs=iErase-1;
% cSeeds=cSeeds(1:mumFoundMotifs, :);
% trainPvalues=trainPvalues(1:mumFoundMotifs, :);
testPvalues=testPvalues(1:mumFoundMotifs, :);

outMotifCell=outMotifCell(1:iErase-1, :);
% validMotifs=testPvalues<=validThreshold;

validMotifs=testPvalues<=0;

testPvalues=testPvalues(validMotifs);
outMotifCell=outMotifCell(validMotifs);
[~, idx]=sort(testPvalues);

% idx=(1:length(testPvalues));

outMotifCell=outMotifCell(idx);
% trainPvalues=trainPvalues(idx);
% consensusSeeds=consensusSeeds(idx, :);

mumFoundMotifs=length(testPvalues);

for im=1:mumFoundMotifs

    outMotifCell{im}.testPvalue=log(mumFoundMotifs)+outMotifCell{im}.testPvalue;
end
% check=1;

end



