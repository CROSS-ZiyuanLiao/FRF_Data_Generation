clc, clear;
fclose('all');

% write first n1 FRFs as training data set
n1 = 5E4;
% total number of FRFs
n2 = 6E4;

folderName = '191220_LeakScale';
% folderName = '191220_MultipleLeaks';
% folderName = '191220_SpeedAndFric';

WD(folderName, n1, n2);

function WD(folderName, n1, n2)
%-----------------------------------------------------------------------------------------
% get all names of data.mat
file = dir(fullfile('.', folderName, '*.mat'));
fileName = {file(:).name}';

nFile = length(fileName);
if nFile == 0
    return;
end

for i = 1: nFile

    % load data
    tempName0 = fileName{i};                               % M1.mat
    tempName1 = fullfile('.', folderName, tempName0);      % ./191220_LeakScale/M1.mat
    data = load(tempName1);
    fprintf([tempName0, ' loaded.\n']);

    tempName2 = tempName0(1: length(tempName0)-4);         % M1
    tempName3 = fullfile('.', folderName, tempName2, 'train');        % ./191220_LeakScale/M1/train
    tempName4 = fullfile('.', folderName, tempName2, 'test');         % ./191220_LeakScale/M1/test

    mkdir(tempName3);
    mkdir(tempName4);

    % write train data
    waitbar1 = waitbar(0, ['Writing training data in ', tempName0, ' into txt files: ', num2str(0), '%...']);
    n = n1;
    for j = 1: n1
        tempName5 = fullfile(tempName3, [num2str(j), '.txt']);
        fileID = fopen(tempName5, 'w');
        leakPipeID = data.leakPipeID(j, :);
        txtLeakPipeID = [];
        for k = 1: length(leakPipeID)
            txtLeakPipeID = [txtLeakPipeID, ',', num2str(leakPipeID(k))];
        end
        txtLeakPipeID = txtLeakPipeID(2: end);
        fprintf(fileID, txtLeakPipeID);
        l = length(data.normalizedFRF(1, :));
        for k = 1: l
            fprintf(fileID, [',', num2str(data.normalizedFRF(j, k))]);
        end
        fclose(fileID);
        waitbar(j/n, waitbar1, ['Writing training data in ', tempName0, ' into txt files: ', num2str(j/n*100), '%...']);
    end
    close(waitbar1);

    % writing test data
    waitbar1 = waitbar(0, ['Writing test data in ', tempName0, ' into txt files: ', num2str(0), '%...']);
    n = n2 - n1;
    for j = (n1+1): n2
        tempName6 = fullfile(tempName4, [num2str(j), '.txt']);
        fileID = fopen(tempName6, 'w');
        leakPipeID = data.leakPipeID(j, :);
        txtLeakPipeID = [];
        for k = 1: length(leakPipeID)
            txtLeakPipeID = [txtLeakPipeID, ',', num2str(leakPipeID(k))];
        end
        txtLeakPipeID = txtLeakPipeID(2: end);
        fprintf(fileID, txtLeakPipeID);
        l = length(data.normalizedFRF(1, :));
        for k = 1: l
            fprintf(fileID, [',', num2str(data.normalizedFRF(j, k))]);
        end
        fclose(fileID);
        waitbar((j-n1)/n, waitbar1, ['Writing test data in ', tempName0, ' into txt files: ', num2str((j-n1)/n*100), '%...']);
    end
    close(waitbar1);

end
%-----------------------------------------------------------------------------------------
end
