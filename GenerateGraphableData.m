% Note that the functions V2VGraphableData and V2IGraphableData both rely
% on raw data files stored on the V2V harddrives, likely stored in the lab
% computer
stepSize = 0.1;
averageWidth = 1;

for testPoint = ['F' 'G' 'H' 'I' 'J']
    fprintf('Working on test point %s', testPoint);
    tic
    V2VGraphableData(testPoint, stepSize, averageWidth);
    toc
end

for testPoint = ['A']
    fprintf('Working on test point %s', testPoint);
    tic
    V2IGraphableData(testPoint, stepSize, averageWidth);
    toc
end
