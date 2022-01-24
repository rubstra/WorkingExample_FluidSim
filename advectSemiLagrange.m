function [ advectedSource ] = advectSemiLagrange( ...
    cellsize,...
    timestep,...
    dissipation,...
    velocityField,...
    source )


%ADVECTSEMILAGRANGE Apply Semi-Lagrangian advection, by use of the 
%   velocityfield, on source data.
%       
%OUTPUT:
%       advectedSource:  A matrix(yNodes, xNodes, [dim]) expressing the  
%       advected state of the source. The source is either the 
%       distancefield or the velocityfield. 
%INPUT:
%       cellSize: grid cell size in m
%       timestep: time per step, deltatime
%       dissipation: friction within fluid
%       velocityField: matrix(yNodes, xNodes, dim) in m/s
%       source: matrix(yNodes, xNodes, [dim]), velocity or distance field,
%       as described elsewhere.
    
%   mSize and nSize are numbers representing number of rows and columns
%   respectively.
    mSize = size(velocityField, 1); % y-direction
    nSize = size(velocityField, 2); % x-direction

    % create two matrices, in which number one, coords(:, :, 1), contains 
    % coordinates from 1 to width, 1:nSize (x-dir), and the other, 
    % coords(:, :, 2) 1 to height, 1:mSize (y-dir).
    coords = zeros(mSize, nSize, 2);
    [coords(:, :, 1), coords(:, :, 2)] = meshgrid(1 : nSize, 1 : mSize);
    
    % track coordinates back to their origin in the last step. Remember
    % that the involved matrices consists of two matrices, one for each
    % dimension (x and y).
    cellVelocity = velocityField / cellsize; % velocity in cells per sec
    previouscoords = coords - cellVelocity * timestep;
    
    %The interpolation is kept as an implemented function, to better keep
    %track of the overall advection process, it is given below. The 
    %interpolatedsource matrix consist of the advected source.
    advectedSource = interpolateSource(previouscoords, source);
    
    % add dissipation (friction)
    advectedSource =  dissipation .*advectedSource;

end


function [ advectedSource ] = interpolateSource( ...
                              previouscoords,...
                              source)
%INTERPOLATESOURCE This function interpolate source values based on the
%coordinates from the previous timestep. Semilagrange advection may
%cause coordinates(particles) to occur outside of the domain. There are 
%many ways to treat this problem of advection, here we choose to force 
%coordinates to within the domain.
%       
%OUTPUT:
%       A matrix(yNodes, xNodes, [dim]) containing the advected state of 
%       the source matrix. That is, source values at the different nodes 
%       are taken from the previous position of these nodes, the position 
%       in the previous timestep.
%INPUT:
%       previouscoordinates: Two matrices, representing the nodes on the 
%       grid, but with values describing the previous position of 
%       respective nodes.

%       source: matrix(yNodes, xNodes, [dim]), velocity or distance field
%       
%       NOTE: The built in Matlab function interp2 used for interpolation
%       within this function, allow use of different interpolation schemes, 
%       such as: 'spline', 'cubic', 'linear' or 'nearest'. use Matlab
%       command "help interp2" for details.
    
    mSize = size(source, 1);
    nSize = size(source, 2);
    
    % add padding, -1 because we don't want fluid to come out of nowhere
    source = padarray(source, [1 1], -1);% adds -1 around the source matrix
    
    %adds 1 to all elements. Interpret this as a need to shift position, as
    %element at index 1,1 lies outside the cells
    previouscoords = previouscoords + 1;
    mPaddedSize = size(source, 1);%number of rows in source
    nPaddedSize = size(source, 2);%number of columns in source
    %replicate all elements along the boundary of previouscoords, and
    %create a new boundary with these elements. The size of the matrix
    %increases by two in each dimension.
    previouscoords = padarray(previouscoords, [1 1], 'replicate');
    
    %keep coordinates in valid range, consequence of backwards advection
    %All coordinates less than one is set to 1. This takes care of
    %particles happen to be above or at the outside the left side of our
    %domain
    previouscoords(previouscoords < 1) = 1;
    xCoords = previouscoords(:, :, 1);% x-coordinates
    yCoords = previouscoords(:, :, 2);% y-coordinates
    %Address elements outside the domian at the right or bottom sides
    xCoords(xCoords > nPaddedSize) = nPaddedSize;
    yCoords(yCoords > mPaddedSize) = mPaddedSize;
    %Bring corrected values back to previouscoords
    previouscoords(:, :, 1) = xCoords;
    previouscoords(:, :, 2) = yCoords;
    
    % interpolate number of dimensions by use of built-in interp2 
    % function (1 dim for distance field, or 2 dims for velocity)
    % note that other interpolation algorithms can be chosen. 
    dimsToInterpolate = size(source, 3);
    advectedSource = zeros(size(source));
    for i = 1 : dimsToInterpolate
        
        advectedSource(:,:,i) = interp2(...
            source(:,:,i), previouscoords(:,:,1),...
            previouscoords(:,:,2),'linear');
    end
%bring matrix back to original dimension. This is the updated sourcefield.
    advectedSource = advectedSource(2 : mSize + 1, 2 : nSize + 1, :);
    
end

