% Set up y-z grid
% Unitless!

clear all
close all

startingY = 1;
endingY = 100;
startingZ = 1;
endingZ = 100;
spacing = .5;

[yCoords, zCoords] = meshgrid(startingY:spacing:endingY, startingZ:spacing:endingZ);

% Matrix representing what is being deposited/built up
% Get the size for later reference
sizeZCoords = size(zCoords);
sedimentationMatrix = zeros(sizeZCoords);
% Note that z coords should be thought of as depth, not height!
%% Establish parameters
description = {'empty', 'stromatolite', 'carbonate', 'shale'};
descriptionCodes = [0, 1, 2, 3];
% What percent of the time is the sediment shale?
percentShale = .25;
% Average sediment thickness, in unitless units
avgSedThickness = .5;
% Given number of cells in length domain, how many do we have to fill to
% hit avgSedThickness?
cellsOfSed = round(avgSedThickness * sizeZCoords(2));
% Sed makeup -- for splitting each pulse of sed into different lithologies
% numShale = round(percentShale * cellsOfSed);
% sedMakeup = [repmat(3, numShale, 1); repmat(2, cellsOfSed - numShale, 1)];
% Local or absolute minima for deposition?
whereDeposited = 'local';
% What is the initial condition?
initialCondition = 'flat';
% Relative scale of noise and elevation
noiseScale = 1;
elevationScale = 1;
% And the growth rate?
avgGrowthRate = .1;
%% Initialize sedimentationMatrix
% We begin by setting the initial condition @ max depth
switch initialCondition
    case 'flat'
        % Set a flat surface, made up entirely of stromatolites
        sedimentationMatrix(zCoords == endingZ) = 1;
    case 'roughStroms'
        % Set a rough initial surface
        % We establish a maximum surface height of pct * the max depth
        maxHeight = .025 * endingZ;
        % Now, we want to know how many y cells we'll need to get a height
        % for
        heights = 0 + (maxHeight-0).*rand(1, sizeZCoords(2));
        sedimentationMatrix(zCoords > endingZ - heights) = 1;
        % Note: think about creating a smooth surface...
    case 'disc'
        % Set a small disc of strom in the center @ the max depth
        % 50% of the total width (i.e., centered)
        sedimentationMatrix(zCoords == endingZ & ...
            yCoords > (.5 * endingY) / 2 & ...
            yCoords < endingY - ((.5 * endingY) / 2)) = 1;
end
%% Run the model
keepGrowing = true;
shalePulse = 0;
sedPulse = 0;
while keepGrowing
    % First, we pulse in some sediment
    % Note that, since coordinates are in depth, we look for local or
    % absolute maxima, not minima when deciding where to put sediment
    % Drop a sed grain, then find maxima, then drop in another grain, and
    % so forth
    % First, draw grains
    % drawnGrains = randsample(sedMakeup, cellsOfSed, false);
    if (shalePulse + 1) / (sedPulse + 1) <= percentShale
        sed = 3;
        shalePulse = shalePulse + 1;
    else
        sed = 2;
    end
    sedPulse = sedPulse + 1;
    for x = 1:cellsOfSed
        % Begin by defining the current surface
        toTest = zCoords;
        toTest(sedimentationMatrix == 0) = nan;
        [currentSurface, currentSurfaceIdxs] = nanmin(toTest, [], 1);
        % Any remaining nans are set to the max z coord
        currentSurface(isnan(currentSurface)) = endingZ;
        switch whereDeposited
            case 'local'
                % In a local situation, we look for local *maxima*
                localMax = islocalmax(currentSurface);
                if sum(localMax) == 0
                    % Note... to handle situations where a local maxima is not
                    % found (i.e., in the case of a flat or single stepped
                    % surface), we just find the max points
                    validMaxIdxs = find(currentSurface == max(currentSurface(:)));
                else
                    validMaxIdxs = find(localMax);
                end
            case 'global'
                % In a global situation, we find just find the *maxima*
                validMaxIdxs = find(currentSurface == max(currentSurface(:)));
        end
        % Draw one valid location
        if length(validMaxIdxs) > 1
            whichMax = randsample(validMaxIdxs, 1);
        else
            whichMax = validMaxIdxs;
        end
        % Now, update
        sedimentationMatrix(max(currentSurfaceIdxs(whichMax) - 1, 1), whichMax) = sed;
    end
    % Now, we want to check if there are any exposed stroms to grow
    surfaceValues = zeros(1, sizeZCoords(2));
    surfaceIdxs = surfaceValues;
    surfaceDepths = surfaceValues;
    for x = 1:sizeZCoords(2)
        thisCol = sedimentationMatrix(:, x);
        findNonZero = find(cumsum(thisCol ~= 0) == 1);
        if isempty(findNonZero)
            surfaceIdxs(x) = sizeZCoords(1);
        else
            surfaceIdxs(x) = findNonZero;
        end
        surfaceValues(x) = thisCol(surfaceIdxs(x));
        surfaceDepths(x) = zCoords(surfaceIdxs(x), x);
    end
    validGrowthLocations = surfaceValues == 1;
    if sum(validGrowthLocations) > 0
        % Now, we grow!
        toGrowFrom = find(validGrowthLocations); 
        for x = 1:length(toGrowFrom)
            % Height on left (if on edge, we just repeat this value
            leftHeight = endingZ - surfaceDepths(min(max(toGrowFrom(x) - 1, 1), sizeZCoords(2)));
            % Height on right
            rightHeight = endingZ - surfaceDepths(min(max(toGrowFrom(x) + 1, 1), sizeZCoords(2)));
            % This height
            thisHeight = endingZ - surfaceDepths(toGrowFrom(x));
            % We need to grow this strom
            growHeight = max(avgGrowthRate + ...
                (noiseScale * randn(1)) + ...
                (elevationScale * (leftHeight + rightHeight - (2*thisHeight))), avgGrowthRate);
            % Now, add to sed matrix
            % Put grow height in terms of spacing
            growHeightCells = floor(growHeight / spacing);
            toFillIn = max(surfaceIdxs(toGrowFrom(x)) - growHeightCells, 1) :...
                max(surfaceIdxs(toGrowFrom(x)), 1);
            sedimentationMatrix(toFillIn, toGrowFrom(x)) = 1;
        end
        % Now, some lateral growth
        preLateralGrowth = sedimentationMatrix;
        carbonateLocs = find(surfaceValues == 2);
        for x = 1:length(carbonateLocs)
            % Okay, you want to look at the cell above this one
            aboveCell = max(surfaceIdxs(carbonateLocs(x)) - 1, 1);
            % Value to left of this above cell?
            leftValue = preLateralGrowth(aboveCell, min(max(carbonateLocs(x) - 1, 1),...
                sizeZCoords(2)));
            % Value on right
            rightValue = preLateralGrowth(aboveCell, min(max(carbonateLocs(x) + 1, 1),...
                sizeZCoords(2)));
            if leftValue == 1 || rightValue == 1
                % If adjacent to strom, add a strom in the cell above
                sedimentationMatrix(max(surfaceIdxs(carbonateLocs(x)) - 1, 1),...
                    carbonateLocs(x)) = 1;
            end
        end
    else
        keepGrowing = false;
    end
    if any(sedimentationMatrix(zCoords == min(zCoords(:))) > 0)
        keepGrowing = false;
    end
end

imagesc(sedimentationMatrix);
axis image;