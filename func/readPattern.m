function ptMotif=readPattern(patternFile, background, options)
    prior=0.1;
if ~options.isPWM
  PatternListChar=importdata(patternFile);
  
    PatternListChar=PatternListChar(2:end);
    PatternListChar=vertcat(PatternListChar{:});
    
    PatternList=double(PatternListChar)-64;
    
    

    priorBg=prior*background{1}/(1+prior);
    ptMotif=cell(size(PatternList,1), 1);
    for iPt=1:size(PatternList,1)
        motifSeq=PatternList(iPt, :);
        PWM=(motifSeq(:)==(1:length(background{1})))/(1+prior)+priorBg.';
        ptLMotif.PWM=PWM;
    
        ptLMotif.cSeed=motifSeq;
        ptLMotif.testPvalue=1;
    
        ptMotif{iPt}=ptLMotif;
    end
else




    opts = detectImportOptions(patternFile); 
    W=options.W;
    opts.DataLines = [2 2+W-1];
    opts.Delimiter=',';
    C = readmatrix(patternFile,opts);
        
    ptMotif=cell(100,1);
    lIdx=1;
    while ~isempty(C)
    
        try
            PWMC=str2num(vertcat(C{:}));
        catch
            break
        end

        PWMC=PWMC./sum(PWMC, 2);
    
        PWM=reshape(PWMC, W, [])+(prior*background{1}).';
        PWM=PWM./(sum(PWM, 2));
        ptLMotif.PWM=PWM;
        [~, cSeed]=max(PWM, [],2);
        cSeed=cSeed.';
        ptLMotif.cSeed=cSeed;
        ptLMotif.testPvalue=1;

        
        ptMotif{lIdx}=ptLMotif;
        opts.DataLines = opts.DataLines+W+1; 
        C = readmatrix(patternFile,opts);

    
        lIdx=lIdx+1;
    end
    
    ptMotif=ptMotif(1:lIdx-1);
end