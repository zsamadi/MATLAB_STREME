function logPvalue=computePvalue(cntInfo, lenInfo, bernoulli)

% computes significance by using either fisher exact test or binomial test
% input arguments:
%   - cntInfo    : vector of posCount and negCount
%   - lenInfo    : vector of Np and NC, posLavg, negLavg average length of sequences, and W 
%   - testMethod : 'fisher' or 'binomial'
%  Output arguments:  
%   - pvalue     : computed pvalue

% [~, idx]=sort(cntInfo(:,1), 'descend');
% 
% cntInfo=cntInfo(idx, :);

% Used in evaluateInitialSeeds, nestedSeedEnrichment, scoreModelPssm





if bernoulli<0

    nPosVec=repmat(lenInfo(1), length(cntInfo(:,1)), 1);

    nNegVec=repmat(lenInfo(2), length(cntInfo(:,1)), 1);


    logPvalue=getFisherLogPvalues(nPosVec, cntInfo(:,1), nNegVec, cntInfo(:,2));


else

    logPvalue=getBinomLogPvalues(cntInfo(:,1), cntInfo(:,2), bernoulli);

    
%     posCount=cntInfo(:,1);
%     negCount=cntInfo(:,2);
%     nPos=lenInfo(1);
%     nNeg=lenInfo(2);
%     posLavg=lenInfo(3);
%     negLavg=lenInfo(4);
%     W=lenInfo(5);
%     Pb=nPos*(posLavg-W+1)/(nPos*(posLavg-W+1)+nNeg*(negLavg-W+1));
%     logPvalue=1- binocdf(posCount-1,posCount+negCount,Pb);  
%     logPvalue=log(max(logPvalue, realmin("double")));
end
% pvalue=gather(pvalue);

