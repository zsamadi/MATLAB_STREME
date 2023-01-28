function outMatchCell=TTest(Target, Querry, options)

motifPWMs=getMotifPWMS(Target);


outMatchCell=cell(length(Querry), 1);
outMatchIndexCell=cell(length(Querry), 1);
minPvalues=zeros(length(Querry), 4);
tomThr=0.05;

for iPt=1:length(Querry)
    motifSeq=Querry{iPt};
    PWM1=motifSeq.PWM;


    [outMatch, outIdx]=tomtom(PWM1, motifPWMs);

    if options.rvp
        [outMatchInv, outIdxInv]=tomtom(PWM1(end:-1:1, :), motifPWMs);
        invFlag=outMatchInv(:,1)<outMatch(:,1);

        outMatch(invFlag, :)=outMatchInv(invFlag, :);
        outIdx(invFlag)=outIdxInv(invFlag);


    end

    outMatchValid=outMatch(:,1)<=tomThr;

    matchResults.stats=[outMatch, outIdx, invFlag];
    matchResults.cSeed=Target{outIdx(1)}.cSeed;
    matchResults.matchs=outIdx(outMatchValid);
    matchResults.mtPvalue=Target{outIdx(1)}.testPvalue;



    outMatchCell{iPt}=matchResults;
    outMatchIndexCell{iPt}= outIdx(outMatchValid);


    minPvalues(iPt, :)=[outMatch(1,:), outIdx(1), invFlag(1)];

    


end


