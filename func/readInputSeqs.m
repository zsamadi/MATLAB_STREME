function [posSeq, negSeq, pHoldSeq,nHoldSeq, dLengs]=readInputSeqs(filename)

opts = detectImportOptions(filename); 


opts.DataLines = [1 4]; % fist block of data
C = readcell(filename,opts);
C=C(end-1:end, :);
pSeqs=cell(2,1);
nSeqs=cell(2,1);
dLengst=cell(2,1);
for iC=1:2
    C1=C{iC,1};
    C1S=find(C1==':');
    if ~isempty(C1S)
        C1=C1(C1S+2:end);
    else
        C1=C{iC,2};
    end


    
    dpattern1=double(C1)-64;
    dpattern1=dpattern1(dpattern1>0);
    
    dSepLoci=find(dpattern1>15);
    if dSepLoci(end)<length(dpattern1)
        dSepLocin=[dSepLoci, length(dpattern1)+1];
        dSepLoci=dSepLocin;
    end
    
    dpattern1(dpattern1>12)=[];
    
    dLengs=dSepLoci-[0, dSepLoci(1:end-1)]-1;
    dLengs=dLengs(dLengs>0);
    
    dStream=mat2cell(dpattern1, 1, dLengs);
    
    numSeqs=length(dStream);
    pSeqs{iC}=dStream(1:numSeqs/2);
    nSeqs{iC}=dStream(numSeqs/2+1:end);
    dLengst{iC}=dLengs;
end

posSeq=pSeqs{1};
pHoldSeq=pSeqs{2};

negSeq=nSeqs{1};
nHoldSeq=nSeqs{2};

dLengs=dLengst{1};


