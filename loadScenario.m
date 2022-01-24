function [ distanceField,...
           velocityField,...
           pressureField,...
           obstacle,...
           forceField] =  loadScenario(xCells, yCells)
%LOADSCENARIO Create domain, obstacles and fluid as defined
%   
%OUTPUT:
%       fluidNodes: A boolean matrix with indices representing nodes
%       (yNodes,xNodes). The values denote water with 1's and everthing
%       else (air) with 0's

%       obstacle: Boolean matrix, as above, obstacles are represented with
%       1's, these are domain borders and obstructions

%       distanceField: A matrix with indices representing nodes(yNodes, 
%       xNodes). The nodal values denote the distances to interfaces, see
%       further descriptions below.

%       velocityField: matrix(yNodes, xNodes,[dim]), velocity of fluid 
%       in m/s. Initial velocity given as 0.5 in the vertical.

%       pressureField: matrix(yNodes, xNodes), pressure at each node,
%       initially given as 0.

%       forceField: bodyforce matrix(yNodes, xNodes, [dim]), force applied to
%           fluid in m/s^2, due to gravitational acceleration




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       %% SCENARIO DESCRIPTION %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   from the initial work of Erler, the current scenario was known as the
%   singleblocker scenario

%   We will be describing the sccenario and prescribing values not to the 
%   cells, but at the nodes. Naturally, the number of nodes will be number 
%   of cells + 1, just visualize any square domain with square grid cells.
%   Each cell has four nodes which it also could be sharing with 
%   neighbouring cells.

    xNodes = xCells + 1;
    yNodes = yCells + 1;
    
%%   GENERATE WATER AND AIR WITHIN DOMAIN %%

%   Everything that is not water is air.
%   The water and air is represented by a boolean matrix, denoting water by
%   1's and not-water (air) by 0's. The boolean matrix describing the
%   position of air-water within the domain is called fluiNodes. 
%   
%   At first it might be challenging to observe that the domain is 
%   generated at the step: fluidNodes = false(yCells, xCells); despite 
%   being the boolean matrix keeping track of air-water whereabouts, it 
%   also represent the internal of the domain. In fact, if you are 
%   considering this simulation as new to the topic, you should be aware of
%   the fact that a single description of the overall domain dosent exist 
%   explicitly! It is described in multiple matrices. So, the typical
%   visual results you may have seen in lectures, or other work is a result
%   of multiple matrices. Below we describe the domain by use of two
%   matrices, (1) fluiNodes as already mentioned and (2) the obstacle
%   matrix, which contains the boundaries and obstacles of the domain.
%   Togehter they form what you would typically think of as the domain, 
%   some boundaries, physical limits, and some fluid, matter, within.
%
%   Initial position and size of water blob, in percentage:
    Xposition = 0.5;%0.5 is the middle for any direction
    Yposition = 0.3;
    radius = 0.25;
        
%   Create domain (internal fluids)
    fluidNodes = false(yNodes, xNodes); % creates boolean matrix with 0's
   
%   place initial water blob into the fluidNodes matrix 
    [x, y] = meshgrid(1 : xNodes, 1 : yNodes);
    x = x - Xposition * xNodes;
    y = y - Yposition * yNodes;
    distToCenter = sqrt(x .* x + y .* y);
    fluidNodes(distToCenter < radius * min(xNodes, yNodes)) = true;
%   The line above assigns the value 1 where we define water within the
%   fluidNodes matrix

        
%   Obstacles and borders
    obstacle = false(yNodes, xNodes); %Boolean matrix with 0's
   
%   Domain boundaries/borders represented with 1's. 
    obstacle(yNodes-1 : yNodes, :) = true; %bottom border
    obstacle(:, 1 : 1+1) = true; %left border
    obstacle(:, xNodes-1 : xNodes) = true;%right border
    
%   Internal obstruction, choose ball obstruction or box obstruction 
%   the latter being the original. 
%   Initial position and size of ball obstruction, in percentage:
        obsXposition = 0.5;
        obsYposition = 0.7;
        obsradius = 0.15;  
%   place obstruction  
    [x, y] = meshgrid(1 : xNodes, 1 : yNodes);
    x = x - obsXposition * xNodes;
    y = y - obsYposition * yNodes;
    distToCenter = sqrt(x .* x + y .* y);
    obstacle(distToCenter < obsradius * min(xNodes, yNodes)) = true;
    
%   Box obstruction (could replace the above ball obstruction)
%   uint32 allow fraction of cells
%   obstacle(uint32(yNodes * 0.65 : yNodes * 0.85), ...
%   uint32(xNodes * 0.4 : xNodes * 0.6)) = true;
    


% we have now generated what we typically consider to be the physical
% domain, with boundaries and fluid.


%%%%%%%%%%%%%%%%%%%%
%% INITIAL VALUES %%    
%%%%%%%%%%%%%%%%%%%%

% Other parameters are described separately in individual matrices.
    
%   PRESSURE
%   initial pressure shall be zero.
    pressureField = zeros(yNodes, xNodes);

%   VELOCITY
%   Generates to velocityfields, x- and y direction. The y direction,
%   representing the vertical, is given an initial velocity of 0.5, while
%   horizontal x direction is kept at 0. Note that the field is stated as
%   yNodes vs xNodes, this is because the x- direction is the column
%   direction and in matrice notation, the first indices refer to the rows
    velocityField(1 : yNodes, 1 : xNodes, 2) = 0.5;
    
%   BODYFORCES
%   Generate two initial forcefields, x- and y- direction,  no forces
    forceField = zeros(yNodes, xNodes, 2);
    
%   APPLY GRAVITY
%   apply force in the vertical direction, hence forcefield number 2
    forceField(:, :, 2) = 9.81;   

    
%   DISTANCEFIELD
%   Calculate distanceField based on the fluidNodes matrix
    distanceField = zeros(size(fluidNodes));%create matrix with zeros
    distancesToWater = (bwdist(fluidNodes) - 0.5) .* -1;
%   bwdist is used to compute the Euclidean distance from any node to
%   the nearest node of nonzero value (remember that within fluidNodes
%   water is represented by 1's and everything else by 0's). Further,
%   when using bwdist all initial nonzero indices become zero,that is,
%   positions initially representing water (by 1's) become zero. A distance
%   of 1 is when a neighbouring element is nonzero, it is one step away. 
%   The -0.5 part is placing the interface in the middle of the water-air 
%   elements. -1 flips the signs.

        
    distancesToAir = (bwdist(~fluidNodes) - 0.5);
%   Writing bwdist(~fluidNodes) as above, gives the distance from the
%   nonzero nodes to the nearest node with value zero.
%   tilde(~) is the "not" operator in Matlab. Subtracting 0.5
%   result in airnodes becoming -0.5 and waternodes = whatever value the
%   nodes had - 0.5
              
    distanceField(~fluidNodes) = distancesToWater(~fluidNodes);
    distanceField(fluidNodes) = distancesToAir(fluidNodes);

end

