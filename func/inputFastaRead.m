function [posSeq, negSeq, pHoldSeq, nHoldSeq]=inputFastaRead(filename, options)

% hFrac=options.hFrac;
% order=options.mkvOrder;
% rvpath=options.rvp;

fastaStruct = fastaread(filename);
numData=length(fastaStruct);

dstream=cell(numData, 1);
dstreamID=zeros(numData, 1);

aLen=options.aLen;

for iStrct=1:numData
    datai=fastaStruct(iStrct).Sequence;
    dstream{iStrct}=datai;
    [~, dataidb]=ismember(datai, options.alphabet);
    dataID=sum(dataidb.*(aLen.^(length(dataidb)-1:-1:0)));
    dstreamID(iStrct)=dataID;
    
end

[~, sID]=sort(dstreamID);

dstream=dstream(sID);

rng('default');

sIDR=randperm(numData);

dstream=dstream(sIDR);


[posSeq, negSeq, pHoldSeq, nHoldSeq]=generateSeqs(dstream,options);
