% *** function addMultiLeaks ***

function addMultiLeaks(obj, pipeID, lFromStart, coeff, pressure)
%-----------------------------------------------------------------------------------------
M = [pipeID, lFromStart, coeff, pressure];
M = sortrows(M, [1 2]);
pipeIndex = unique(M(:, 1));
nLeakPipe = length(pipeIndex);
c = cell(nLeakPipe, 1);
for i = 1: nLeakPipe
    c{i} = M(M(:, 1)==pipeIndex(i), :);
    c{i}(:, 5) = 0;
end
for i = 1: nLeakPipe
    np = obj.nPipe;
    c{i}(1, 5) = c{i}(1, 2);
    n = size(c{i}, 1);
    if n >= 2
        for j = 2: n
            np = np + 1;
            c{i}(j, 1) = np;
            c{i}(j, 5) = c{i}(j, 2) - c{i}(j-1, 2);
        end
    end
    pipeID = c{i}(:, 1);
    lFromStart = c{i}(:, 5);
    coeff = c{i}(:, 3);
    pressure = c{i}(:, 4);
    for j = 1: size(c{i}, 1)
        obj.addOneLeak(pipeID(j), lFromStart(j), coeff(j), pressure(j));
    end
end
%-----------------------------------------------------------------------------------------
end

