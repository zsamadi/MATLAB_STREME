function rcData=revCmp(idata, aLen)
if any(size(idata)==1)
    rcData=revCmps(idata, aLen);
else
    rcData=zeros(size(idata));
    for i=1:size(idata, 1)
        rcData(i, :)=revCmps(idata(i, :), aLen);
    end
end
