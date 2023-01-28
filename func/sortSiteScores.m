function [positivePWMS, iN]=sortSiteScores(appxSites, positivePWMS, is_positive)
    

%     positivePWMS=round(positivePWMS,10);
    [positivePWMS, iN]=sortrows([positivePWMS, -is_positive,-appxSites], 'descend');

    positivePWMS=positivePWMS(:, 1);


