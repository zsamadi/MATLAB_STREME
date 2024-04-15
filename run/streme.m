function streme(filename, options)


arguments

    filename='exData.fasta';
    options.outFolder='output';
    options.NEVAL=25;
    options.NREF=4;
    options.prior=0.01;
    options.REFPrior=0.1;
    options.ENRCHPrior=0.1;
    options.W=6;
    options.nRefIter=20;
    options.threshold=0.01;
    options.patience=3;
    options.evalue=false;
    options.nmotifs=0;
    options.alphabet='ACGT';
    options.mkvOrder=0;
    options.rvc=true;
    options.hFrac=0.1;


end

outputFolderName=strcat(options.outFolder, '/');

cd(fileparts(which(mfilename)));
cd('..\')





addpath(genpath(pwd))



timerValue=tic;

sqOptions.hFrac=options.hFrac;
sqOptions.mkvOrder=options.mkvOrder;
sqOptions.rvc=options.rvc;
sqOptions.isChar=true;
sqOptions.aLen=length(options.alphabet);
sqOptions.alphabet=options.alphabet;

fprintf('Creating Negative Control and Hold-Out Data... \n ');

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
filename=strcat(outputFolderName, 'MatSTREME.txt');
writeTextOutput(filename, extMotif,bkg,  textOut,commandText, elapsedTime, options.alphabet);


%% Seqlogo
numExtMotifs=length(extMotif);

for iexMotif=1:numExtMotifs

    exMotifi=extMotif{iexMotif};
    exMotifi.background=background;

    [~, hfig] = seqlogoGen(exMotifi, options.alphabet);
    [~, exMotifi.cSeed]=max(exMotifi.PWM);

    % cSeedStr = strjoin(string(exMotifi.cSeed), '_');

    figname=strcat(outputFolderName, 'mLogo',num2str(iexMotif), '.jpg');
    saveas(hfig,figname)

end

%% PValues

pvalueVec=zeros(numExtMotifs, 1);
scoreThrVec=zeros(numExtMotifs, 1);

for iexMotif=1:numExtMotifs

    exMotifi=extMotif{iexMotif};
    pvalueVec(iexMotif)=exMotifi.testPvalue;
    scoreThrVec(iexMotif)=exMotifi.scoreThr;

end

figure('visible','off');
plot(pvalueVec, '-kv', 'LineWidth',1)
grid on

xlabel('motif index')
ylabel('evalue')

figname=strcat(outputFolderName,'pvalueW', num2str(options.W), '.jpg');
saveas(gcf,figname)



