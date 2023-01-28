function [positivePWMS, iN]=sortNodeScores(positivePWMS)

%     positivePWMS=round(positivePWMS,10);
    [positivePWMS, iN]=sort(positivePWMS, 'descend');
