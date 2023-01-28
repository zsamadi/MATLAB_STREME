function seedsEvaled=evaluateInitialSeeds(infoMat,Nt)

% sorts WMers of length W  by pvalues
% input arguments:
%   - infoMat: a matrix containing 
%     - uniWMers: all seeds of length W
%     - posCnts: counts of Wmer seeds in primary data
%     - negCnts: counts of Wmer seeds in control data
%  Output arguments:
%   - seedsEvaled: struct for sorted input WMers
%      - seedsEvaled.seeds: sorted WMers NEVAL of which to be passed to
%      refinement step
%      - seedsEvaled.pvalues: significance obtained by fisher exact test

nPos=Nt(1);
nNeg=Nt(2);

W=size(infoMat, 2)-2;
uniWMers=infoMat(:, 1:W);
posCnts=infoMat(:, W+1);
negCnts=infoMat(:, W+2);

% compute pvaluse using fisher exact test, computation performed on GPU
posInfoGPU=gpuArray([posCnts, negCnts]);
lenInfoGPU=gpuArray([nPos, nNeg]);
fisherTest=computePvalue(posInfoGPU, lenInfoGPU, 'fisher');
fisherTest = gather(fisherTest);


% sorting with significance test results

[fisherTestSort, fisherPerm] = sort(fisherTest);
uniWMers=uniWMers(fisherPerm, :);
posCnts=posCnts(fisherPerm);
negCnts=negCnts(fisherPerm);


% sorted seeds by pvalue
seedsEvaled.seeds=uniWMers;
seedsEvaled.seedsChar=[char(uniWMers+64), num2str([posCnts, negCnts, log(fisherTestSort)])];

seedsEvaled.Cnts=[posCnts, negCnts, log(fisherTestSort)];
seedsEvaled.pvalues=fisherTestSort;
