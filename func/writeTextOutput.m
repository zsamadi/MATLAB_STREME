function done=writeTextOutput(filename, motifsCell,bkg, textIn,commandText, elapsedTime, alphabet)

PWM0=bkg.PWM0;


fileID = fopen(filename,'w');
if ~isempty(motifsCell)
    [alength, w]=size(motifsCell{1}.PWM);
else
    alength=1;
    w=1;
end

starLine=repelem('*', 1, 50);
starLine=strcat(starLine, '\n');

fprintf(fileID,starLine);
fprintf(fileID,'MATLAB STREME - Sensitive, Thorough, Rapid, Enriched Motif Elicitation\n');
fprintf(fileID,starLine);
fprintf(fileID,'MEME version 5.5.0 (Release date: Wed Sep 7 14:18:26 2022 -0700)\n');
fprintf(fileID,starLine);
fprintf(fileID,'\n');

fprintf(fileID,'\n');
fprintf(fileID,starLine);

fprintf(fileID,'REFERENCE\n');
fprintf(fileID,starLine);

fprintf(fileID,'Timothy L. Bailey,\n');
fprintf(fileID,'"STREME: accurate and versatile sequence motif discovery",\n');

fprintf(fileID,'Bioinformatics, Mar. 24, 2021.\n');
fprintf(fileID,starLine);

fprintf(fileID,'Zainalabedin Samadi, Amjad Askary\n');
fprintf(fileID,'"Spatial motifs reveal patterns in cellular architecture of complex tissues",\n');

fprintf(fileID,'biorxiv, Apr. 21, 2024.\n');
fprintf(fileID,starLine);

fprintf(fileID,'\n');

fprintf(fileID,'\n');
fprintf(fileID,starLine);
fprintf(fileID,sprintf('ALPHABET "%s"\n', alphabet));
fprintf(fileID,starLine);

fprintf(fileID,'Background letter frequencies\n');

for ipw=1:length(PWM0)

fprintf(fileID,'%s:%1.5f,', alphabet(ipw), PWM0(ipw));
end
fprintf(fileID,'\n');
fprintf(fileID,'\n');

for im=1:length(motifsCell)
    outMotif=motifsCell{im};
    seedi=outMotif.cSeed;
    seediChar=alphabet(seedi);

%     secSeedi=outMotif.secSeeds(1, :);
%     secSeediChar=char(secSeedi+64);


    strText=strcat("MOTIF  ", num2str(im), '-', seediChar,' MSTREME-', num2str(im), '\n');
    fprintf(fileID,strText);

    testPvalue=outMotif.testPvalue/log(10);

    ep=floor(testPvalue);
    ei=10^(testPvalue-ep);
    eText=strcat(num2str(ei), 'e', num2str(ep));

    fprintf(fileID,'letter-probability matrix: alength= %d w= %d nsites= %d E= %s\n', alength, w, outMotif.nsites, eText);



        PWM1=outMotif.PWM;


        for ii = 1:size(PWM1,2)
            fprintf(fileID,'%1.7f\t',(PWM1(:,ii)).');
            fprintf(fileID,'\n');
        end
fprintf(fileID,'\n');

%     writematrix(PWMSE,filename,'WriteMode','append','Delimiter','tab')

end
fprintf(fileID,starLine);

fprintf(fileID,textIn);
fprintf(fileID,starLine);

% commandText=strcat("COMMAND:    ",commandText, '\n');
% fprintf(fileID,commandText);
% fprintf(fileID,starLine);

pcName = getenv('COMPUTERNAME');

cpuRext=strcat("CPU:		",pcName, '\n');

fprintf(fileID,cpuRext);

fprintf(fileID,starLine);
fprintf(fileID,'FINALTIME:	%5.2f seconds\n', elapsedTime);

fprintf(fileID,starLine);


fclose(fileID);

done=0;