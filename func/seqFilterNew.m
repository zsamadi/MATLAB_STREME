function [seqFiltered, erasedSites]=seqFilterNew(seq, filterSpecs)

% Erasing PWMS motif obtained from enrichment step with score above threshold
% Input arguments:
%   - seq       : input sequence
%   - PWMS      : Refined PWM Score matrix
%   - threshold : score threshod for erasing
% Output arguments:
%  - seqFiltered:input sequence with motif erased

% 
% filterSpecs.back=seqData.back;
% filterSpecs.order=seqData.order;

PWMS=filterSpecs.PWMS;
threshold=filterSpecs.scoreThr;

if (filterSpecs.rvc)
    numSeqs=length(seq)/2;


else
    numSeqs=length(seq);
end

seqH=horzcat(seq{1:numSeqs});




seqLens=zeros(numSeqs, 1);
for iSl=1:numSeqs
    seqLens(iSl)=length(seq{iSl});
end
W=size(PWMS, 2);


cnvOptions.numPSeqs=numSeqs;
cnvOptions.W=W;
% cnvOptions.rvc=filterSpecs.rvc;
cnvOptions.rvc=false;

cnvOptions.allMers=true;
cnvOptions.pnSeq=false;

[XNmerso , merEyesOn]=convertNMer(seqH,seqLens, cnvOptions);

% sites=XNmerso(:, end-1);

% merEyesOnt=merEyesOn;
% XNmersot=XNmerso;

% if (filterSpecs.rvc)
%     rvcFlag=XNmerso(:, 1)>XNmerso(:, W);
%     XNmerso(rvcFlag, 1:W)=XNmerso(rvcFlag, W:-1:1);
%     merEyesOnW=(reshape(merEyesOn, W, [])).';
%     merEyesOnW(rvcFlag, :)=merEyesOnW(rvcFlag, end:-1:1);
%     merEyesOnW=merEyesOnW.';
%     merEyesOn=merEyesOnW(:);
% 
% end




[XNmerUo, ~, jxo]=unique(XNmerso(:, 1:end-2), 'rows');

[PWMScoreo, PWMScoreLettero]=scoreWords(XNmerUo, PWMS, filterSpecs);


eraseFlago=(PWMScoreo>=threshold);

eraseFlago=eraseFlago(jxo);
%     eraseFlagMat(:,iErFlag)=eraseFlag;

PWMScoreLettero=PWMScoreLettero(jxo, :);



eraseFlagLettero=(PWMScoreLettero>0);
eraseFlagLettero=eraseFlagLettero&eraseFlago;

eraseFlagLetteron=eraseFlagLettero.';

eraseFlagLetteron=eraseFlagLetteron(:);

%Also erase reverse path
if filterSpecs.rvc
    [PWMScoreorvc, PWMScoreLetterorvc]=scoreWords(XNmerUo, PWMS(end:-1:1, end:-1:1), filterSpecs);
    
    eraseFlagorvc=(PWMScoreorvc>=threshold);
    
    eraseFlagorvc=eraseFlagorvc(jxo);
    PWMScoreLetterorvc=PWMScoreLetterorvc(jxo, :);
    eraseFlagLetterorvc=(PWMScoreLetterorvc>0);
    eraseFlagLetterorvc=eraseFlagLetterorvc&eraseFlagorvc;
    eraseFlagLetteronrvc=eraseFlagLetterorvc.';
    
    eraseFlagLetteronrvc=eraseFlagLetteronrvc(:);
    eraseFlagLetteron=eraseFlagLetteron|eraseFlagLetteronrvc;
end





eraseFlagLetteront=zeros(numel(seqH), 1);
temp=accumarray(merEyesOn, eraseFlagLetteron);

eraseFlagLetteront(1:length(temp))=temp;

% erasedSites=reshape(eraseFlagLetteront, size(seq, 2), []);
% 
% erasedSites=erasedSites.';
eraseLetters=eraseFlagLetteront>0;

seqFilteredERS=seqH;
seqFilteredERS(eraseLetters)=0;
seqFiltered=mat2cell(seqFilteredERS, 1,seqLens);

if (filterSpecs.rvc)
    seqFilteredrvc=mat2cell(revCmp(seqFilteredERS, size(PWMS, 1)), 1,seqLens(end:-1:1));
    seqFilteredrvc=seqFilteredrvc(end:-1:1);
    seqFiltered=(horzcat(seqFiltered, seqFilteredrvc)).';


end


seqNum=repelem((1:length(seqLens)).', seqLens-W+1);

erasedSites=accumarray(seqNum, eraseFlago);


