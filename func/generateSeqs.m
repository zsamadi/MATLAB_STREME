function [posSeq, negSeq, pHoldSeq, nHoldSeq]=generateSeqs(dstream,options)
% from input sequence, generates control sequence by shuffling, and
% separates 10% of the data as hold out data

% input arguments:
    % - dstream: input sequence
% % Output arguments:
    % - posSeq90: output primary sequence(90% of input sequence)
    % - negSeq90: output control sequence(90% of shuffled input sequence)
    % - pHoldSeq: output primary hold out sequence
    % - nHoldSeq: output control hold out sequence

hFrac=options.hFrac;
order=options.mkvOrder;

rvpath=options.rvp;

numSeqs=length(dstream);

shStream=cell(numSeqs, 1);

 if (rvpath)
    dstreamInv=cell(numSeqs, 1);
    shStreamInv=cell(numSeqs, 1);
 end

 dstreamDouble=cell(numSeqs, 1);
for iSp=1:numSeqs

    if (options.isChar)
        dstreamChar=dstream{iSp};
        [~, dstreamDouble{iSp}]=ismember(dstreamChar, options.alphabet);

    else
        dstreamChar=char(dstream{iSp}+64);
        dstreamDouble{iSp}=dstream{iSp};
        

    end



    shuffled=ushuffle(dstreamChar, order+1);
  
    [~,shStream{iSp}]=ismember(shuffled, options.alphabet);

   if (rvpath)
       dstreamInv{iSp}=revCmp(dstreamDouble{iSp}, options.aLen); 
       shStreamInv{iSp}=revCmp(shStream{iSp}, options.aLen);
    end
            




end

% 
% [~, sID]=sort(dstreamID);
% 
% dstreamDouble=dstreamDouble(sID);
% shStream=shStream(sID);
% rng('default');
% 
% sIDR=randperm(numSeqs);
% 
% dstreamDouble=dstreamDouble(sIDR);
% shStream=shStream(sIDR);




numHSeq=floor(numSeqs*hFrac);

pHoldSeq=dstreamDouble(1:numHSeq);
posSeq=dstreamDouble(numHSeq+1:end);

nHoldSeq=shStream(1:numHSeq);
negSeq=shStream(numHSeq+1:end);


 if (rvpath)
    posSeq=vertcat(posSeq, dstreamInv(numHSeq+1:end));
    negSeq=vertcat(negSeq, shStreamInv(numHSeq+1:end));
    pHoldSeq=vertcat(pHoldSeq, dstreamInv(1:numHSeq));
    nHoldSeq=vertcat(nHoldSeq, shStreamInv(1:numHSeq));
 end





