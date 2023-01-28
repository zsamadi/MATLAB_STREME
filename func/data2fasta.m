function result=data2fasta(posSeq90, negSeq90)

% Write Fasta file to be analyzed by official STREME 
numDStreamFasta=size(posSeq90, 1);
posFasta=char(posSeq90+64);
negFasta=char(negSeq90+64);

% numDStreamFasta=5;
rng('default');
% dstreamFasta=dstream(randi(numDStream, numDStreamFasta, 1), :);
% load('dataSaved.mat')
% 
% OIndex=randi(numel(dstreamFasta)-2, floor(numel(dstreamFasta)), 1);
% dstreamFasta(OIndex)='O';
% dstreamFasta(OIndex+1)='O';
% dstreamFasta(OIndex+2)='O';

posData(numDStreamFasta) = struct();
negData(numDStreamFasta) = struct();

posWriteFileName=strcat('posSeq',num2str(numDStreamFasta),'.fa');
negWriteFileName=strcat('negSeq',num2str(numDStreamFasta),'.fa');


posWriteFile=strcat('C:\cygwin64\home\askar\', posWriteFileName);
negWriteFile=strcat('C:\cygwin64\home\askar\', negWriteFileName);

for idat=1:size(posFasta, 1)
    posData(idat).Sequence = posFasta(idat, 1:8);
    posData(idat).Header =strcat('sequence ', num2str(idat));
    negData(idat).Sequence = negFasta(idat, 1:8);
    negData(idat).Header =strcat('sequence ', num2str(idat));
end
if exist(posWriteFile, 'file')==2
  delete(posWriteFile);
end

fastawrite(posWriteFile, posData)


if exist(negWriteFile, 'file')==2
  delete(negWriteFile);
end

fastawrite(negWriteFile, negData)
result=0;
%%%
