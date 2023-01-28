function motifPWMs=getMotifPWMS(Target)

motifPWMCell=cell(length(Target), 1);



for iMEx=1:length(Target)
    motifPWM=(Target{iMEx}.PWM);

    motifPWMCell{iMEx}=motifPWM.';
end

motifPWMs=cell2mat(motifPWMCell);