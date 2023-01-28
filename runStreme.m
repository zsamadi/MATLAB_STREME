function runStreme(filename, options)


arguments
    filename='empty.txt'
    options.NEVAL=25;
    options.NREF=4;
    options.prior=0.01;
    options.REFPrior=0.1;
    options.ENRCHPrior=0.1;
    options.W=4;
    options.nRefIter=20;
    options.threshold=0.05;
    options.patience=3;
    options.evalue=false;
    options.nmotifs=0;
    options.alphabet='ABCDEFGHIJKL';
    options.mkvOrder=0;
    options.rvp=false;
    options.hFrac=0.1;
    options.isPWM=false;
    options.patternFile='patternList.txt';

end

if(~isdeployed)
    cd(fileparts(which(mfilename)));
end
addpath(genpath(pwd))



timerValue=tic;

sqOptions.hFrac=options.hFrac;
sqOptions.mkvOrder=options.mkvOrder;
sqOptions.rvp=options.rvp;
sqOptions.isChar=true;
sqOptions.aLen=length(options.alphabet);


[posSeq, negSeq, pHoldSeq, nHoldSeq]=inputFastaRead(filename, sqOptions);


bernoulli=getBernoulli(posSeq, negSeq, options.W);

options.prior=0.01;
options.REFPrior=0.1;
options.ENRCHPrior=0.1;
options.wMin=options.W;
options.wMax=options.W;
options.MinSeedWidth=3;

options.maxPalED=5;
options.minPalRatio=0.85;
options.isPal=false;
options.bernoulli=bernoulli;


seqData.pSeq=posSeq;
seqData.nSeq=negSeq;
seqData.pHSeq=pHoldSeq;
seqData.nHSeq=nHoldSeq;


[extMotif,textOut, commandText, background]=gmtStreme(seqData, options);






elapsedTime=toc(timerValue);


bkg.PWM0=background{1};
bkg.mkvOrder=options.mkvOrder;
filename='output/MatSTREME.txt';
writeTextOutput(filename, extMotif,bkg,  textOut,commandText, elapsedTime);


%% Tomtom Test


rpOtions.W=options.W;
rpOtions.isPWM=options.isPWM;
patternFile=options.patternFile;


ptMotif=readPattern(patternFile, background, rpOtions);
numPtMotif=length(ptMotif);

toptions.background=background;
toptions.rvp=true;

outMatchMP=TTest(extMotif, ptMotif, toptions);
outMatchPM=TTest(ptMotif,extMotif,toptions);


minPvalueMP=zeros(numPtMotif,1);
enrchPvalueMP=zeros(numPtMotif,1);
minPvalueIndexMP=zeros(numPtMotif,2);



for iPt=1:numPtMotif

    minPvalueMP(iPt)=outMatchMP{iPt}.stats(1,1);
    minPvalueIndexMP(iPt, :)=[iPt, outMatchMP{iPt}.stats(1,3)*(-1)^outMatchMP{iPt}.stats(1,4)];

    enrchPvalueMP(iPt)=outMatchMP{iPt}.mtPvalue;


end

[minPvalueMP, iSt]=sort(log10(minPvalueMP), 'descend');
enrchPvalueMP=enrchPvalueMP(iSt);
minPvalueIndexMP=minPvalueIndexMP(iSt, :);




numExtMotif=length(extMotif);
minPvaluePM=zeros(numExtMotif,1);
minPvalueIndexPM=zeros(numExtMotif,2);
enrchPvaluePM=zeros(numExtMotif,1);

for iMt=1:numExtMotif



    minPvaluePM(iMt)=outMatchPM{iMt}.stats(1,1);
    minPvalueIndexPM(iMt, :)=[iMt, outMatchPM{iMt}.stats(1,3)*(-1)^outMatchPM{iMt}.stats(1,4)];

    enrchPvaluePM(iMt)=extMotif{iMt}.testPvalue;

end

[minPvaluePM, iSt]=sort(log10(minPvaluePM), 'descend');
enrchPvaluePM=enrchPvaluePM(iSt);
minPvalueIndexPM=minPvalueIndexPM(iSt, :);





%% Plot the results

figure
hold on

[~, enrchPvalueMPi]=sort(enrchPvalueMP);
[~, enrchPvalueMPi]=sort(enrchPvalueMPi);

colorMP=colormap(jet(length(enrchPvalueMPi)));
colorMP=colorMP(enrchPvalueMPi, :);


scatter((1:length(minPvalueMP)),minPvalueMP,[],colorMP,'filled','v')
text((1:length(minPvalueMP)),minPvalueMP+0.2,num2str(minPvalueIndexMP))


line((1:length(minPvalueMP)),minPvalueMP,'LineStyle','--', 'LineWidth',1.1)


grid on

title('Minimum LogPValues(Best Match) of Pattern Motifs in Extracted Motifs')
xlabel('Motif index')
xlim([0 numPtMotif+1])
ylabel('Log pvalue')
c=colorbar;
c.Ticks=[0, 1];
c.TickLabels = {'low evalue','high evalue'};



% showing how many false positives MATLAB has

figure
hold on



[~, enrchPvaluePMi]=sort(enrchPvaluePM);
[~, enrchPvaluePMi]=sort(enrchPvaluePMi);

colorPM = colormap(jet(length(enrchPvaluePMi)));
colorPM=colorPM(enrchPvaluePMi, :);
scatter((1:length(minPvaluePM)),minPvaluePM,[],colorPM,'filled','v')
text((1:length(minPvaluePM)),minPvaluePM+0.2,num2str(minPvalueIndexPM))

line((1:length(minPvaluePM)),minPvaluePM,'LineStyle','--', 'LineWidth',1.1)

grid on
title('Minimum LogPValues(Best Match) of Extracted Motifs in Pattern Motifs')
xlabel('Motif index')
xlim([0 numExtMotif+1])
ylabel('Log pvalue')

c=colorbar;
c.Ticks=[0, 1];
c.TickLabels = {'low evalue','high evalue'};