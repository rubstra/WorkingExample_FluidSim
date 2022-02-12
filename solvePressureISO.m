function [ pressure ] = solvePressureISO( xCells,yCells, velocity, distanceField, obstacle )
%SOLVEPRESSURE Apply pressure solver with finite elements method on 
%   triangular elements. The solver allows higher order polynomials to be 
%   defined on the elements, up to forth order exact solution.
%   
%   For more information on the solver and finite elements methods a
%   recommended book isThe Finite Element Method: Theory, Implementation,
%   and Applications by Larson and Bengzon. The book contain several minor
%   errors in the examples, but explain the theory in a very understandable
%   manner.
%
%   OUTPUT: pressure matrix (yNodes, xNodes) containing the pressure at 
%           each node, solved by the finite element method
%
%    INPUT: DistanceField: A matrix with indices representing nodes(yNodes, 
%           xNodes). The nodal values denote the distances to interfaces, see
%           further descriptions in loadscenario or other.

%           Velocity: velocity field matrix(yNodes, xNodes,[dim]), velocity 
%           of fluid in m/s at the nodes.
%
%           Obstacle: Boolean matrix, obstacles are represented with
%           1's, these are domain borders and obstructions



%%%%%%%%%%%%%%%%%%%%   
%% initialization %%
%%%%%%%%%%%%%%%%%%%%


% Note that we use y,x only when accessing matrizes from caller, internal 
% is x,y


%renaming for clarity, easier to follow the procedure
    numCellsX = xCells; %number of cells in x-dir
    numCellsY = yCells; %number of cells in y-dir
%renaming for clarity, easier to follow the procedure   
    velocityFieldX = velocity(:,:,1); %renaming for 
    velocityFieldY = velocity(:,:,2);

% Allocating matrices
    % sparse matrix constructor parameters
    Kx = [];
    Ky = [];
    Kval = [];
    
    % pressure vector as sparse matrix
    fx = [];
    fval = [];
    
    % Get fluid cells and interface nodes
    % Convert a distance field, based on nodes, to a boolean matrix of 
    % cells. See getFluidCells for further description
    fluidCells = getFluidCells(distanceField);
    % See getInterfaceNodes for further description
    interfaceNodes = getInterfaceNodes(distanceField, obstacle, true);
    
    %Defining and allocating 
    nodeIndex = 0;
    nodeCoordToIndex = zeros(size(distanceField, 1), size(distanceField, 2));
    nodeIndexToCoord = zeros(size(distanceField, 1) * size(distanceField, 2), 2);
    nodesInElements = zeros(3,size(distanceField, 1) * size(distanceField, 2));
    
    % if we have a cell corresponding to 2,2 (row, column) then the node
    % offset would give the nodes surrounding that cell. Take the first row
    % vector 0 0 and add by elements to 2 2 it equals 2 2. Then node 2 2 is
    % the upper left node in the cell. Likewise 1,1 result in node 3,3
    % which is the bottom right node of the cell. The offsets refer to a
    % real coordinate system, not the matrix elements
    nodeOffsets = [0 0; 0 1; 1 1; 1 0];
    element = 1;
    

    
    %%%%%%%%%%%%%%%%%%%%%%
    %% pressure solving %%
    %%%%%%%%%%%%%%%%%%%%%%
    
    % for each cell in y-dir
    for cellY = 1 : numCellsY 
        % for each cell in x-dir
        for cellX = 1 : numCellsX
            % surrounded by 4 nodes
            % skip cells without water
            if ~fluidCells(cellY, cellX)
                continue;
            end
            
            %coordinates for the current iteration
            currentCoord = [cellX cellY];
            
            %Allocating/defining vectors
            wTildeEX = [0; 0; 0; 0];
            wTildeEY = [0; 0; 0; 0];
                     
                     
            % We need to gather all information for cell
            for corner = 1 : 4 % run onces per corner/node of the cell
                % gives corner coordinates, see nodeOffsets described above
                cornerCoord = currentCoord + nodeOffsets(corner, :);
                
                % Temporarily store nodal x and y components of velocity in
                % respective vectors
                wTildeEX(corner) = velocityFieldX(cornerCoord(2), cornerCoord(1));
                wTildeEY(corner) = velocityFieldY(cornerCoord(2), cornerCoord(1));
                
                % add this node to indexing
                if nodeCoordToIndex(cornerCoord(2), cornerCoord(1)) == 0
                    nodeIndex = nodeIndex + 1;
                    
                    nodeCoordToIndex(cornerCoord(2), cornerCoord(1)) = nodeIndex;
                    nodeIndexToCoord(nodeIndex, :) = cornerCoord;
                end                            
            end
            
            %Defining the nodes belonging to each triangle element, saving
            %the node indeces (x;y;z) of an element v in column v of the
            %nodesInElements matrix. This cause matrix dimensions of (nodes
            %per element, number of elements) or in this case (3, number of
            %elements)
            nodeCoord = repmat(currentCoord, [4,1])+ nodeOffsets;%create 
            % vector containing the 4 node coordinates belonging to the 
            % current cell. These coordinates give rise to two triangles
            %first triangle 1/2
            nodesInElements(:,element) =...
                   [nodeCoordToIndex(nodeCoord(1,2),nodeCoord(1,1)),... 
                   nodeCoordToIndex(nodeCoord(2,2),nodeCoord(2,1)),...
                   nodeCoordToIndex(nodeCoord(4,2),nodeCoord(4,1))];
            %second triangle 2/2
            nodesInElements(:,element+1) =...
                   [nodeCoordToIndex(nodeCoord(2,2),nodeCoord(2,1)),... 
                   nodeCoordToIndex(nodeCoord(3,2),nodeCoord(3,1)),...
                   nodeCoordToIndex(nodeCoord(4,2),nodeCoord(4,1))];
               
              

%ASSEMBLING STIFFNESS MATRIX AND LOAD VECTOR

% The current procedure for assembling the stiffness assembler and load 
% vector is based on the code presented in the mentioned book of Larson and
% Bengzon, see the above introduction, chapter 8, especially pages 211-217 
% in the 2013 version.


for t=1:2 %run for both triangles  t = 1:2

    %Stiffnessassembler ISO
    
[rspts,qwgts]=Gausspoints(2); % quadrature routine. Retrieves references 
%triangle coordinates and quadrature weights

loc2glb = nodesInElements(1:3,(element+t-1)); % node numbers of element
x = nodeIndexToCoord(loc2glb,2); % get x coordinates
y = nodeIndexToCoord(loc2glb,1); % get y coordiantes
AK=zeros(3,3); % allocate local element stiffness matrix

for q=1:length(qwgts) % loop for the number of quadrature weights
r=rspts(q,1); % refrerence triangle r coordinate
s=rspts(q,2); % refrerence triangle s coordinate
[dSdx,dSdy,detJ]=Isopmap(x,y,r,s); %The isoparametric map returns the
%determinant of the Jacobian as well as the partial derivatives of the
%shape function
wxarea=qwgts(q)*detJ/2; % weight times area (detJ = 2*area of triangle K)
AK=AK+(dSdx*dSdx'...
+dSdy*dSdy')*wxarea; % local element stiffness matrix AK
end

%Add local element matrix to global matrix. The global matrix is structured
%as columnvector with associated vectors containing the node indices. The
%global element column vector is Kval and the associated node vectors are
%Kx and Ky.

%Allocate vectors
     newKx = zeros(9, 1);
     newKy = zeros(9, 1);
     newKval = zeros(9, 1);
     newFx = zeros(3, 1);     
     
     vectorIndex = 1;
            
            for fromCorner = 1 : 3
                fromCornerNodeIndex = loc2glb(fromCorner);
                
                for toCorner = 1 : 3
                    toCornerNodeIndex = loc2glb(toCorner);
                
                    newKx(vectorIndex) = fromCornerNodeIndex;
                    newKy(vectorIndex) = toCornerNodeIndex;
                    newKval(vectorIndex) = AK(fromCorner, toCorner);
                    vectorIndex = vectorIndex + 1;
                end
                %This is for the load vector, to be calculated below
                newFx(fromCorner) = fromCornerNodeIndex;
            end
            Kx = [Kx; newKx];%adding newKx to the groing columnvector Kx
            Ky = [Ky; newKy];
            Kval = [Kval; newKval];
            fx = [fx; newFx];



%Load assembler

bK=0;
wEX = wTildeEX([1+t-1,2+t-1,4]);
wEY = wTildeEY([1+t-1,2+t-1,4]);
f(1,:) = dSdy'*wEX+dSdx'*wEY;
f(2,:) = dSdy'*wEX+dSdx'*wEY;
f(3,:) = dSdy'*wEX+dSdx'*wEY;
bK = [f(1,1);f(2,1);f(3,1)]*wxarea; % element load vector

fval = [fval; bK];



end
element = element+2;
        end
    end
 

% Remove boundary nodes from matrix and vector, pressure must be 0
% keep them if they are in obstacles
% see function further below
    [Kx, Ky, Kval, fx, fval] = removeBoundaryNodes(...
          Kx, Ky, Kval, fx, fval, nodeIndexToCoord, interfaceNodes);
      
    % map node indieces to matrix collumn indieces
    activeNodeIndieces = unique(fx);
    numActiveNodes = size(activeNodeIndieces, 1);
    nodeIndexToCollumn(activeNodeIndieces) = 1:numActiveNodes;
    Kx = nodeIndexToCollumn(Kx);
    Ky = nodeIndexToCollumn(Ky);
    fx = nodeIndexToCollumn(fx);
    
    % get pressurefield
    K = sparse(Kx, Ky, Kval);
    f = sparse(fx, ones(size(fx)), fval);
    
    pressureAtNodes = K \ f;    % solve
    
    % should never be necessary in a working system
    %pressureAtNodes(isnan(pressureAtNodes)) = 0; 
    %pressureAtNodes(isinf(pressureAtNodes)) = 0;
    
    % save resulting pressure values from node vector to matrix
    pressure = zeros(size(distanceField));
    currNodeIndex = 1 : numActiveNodes;
    currNodeY = nodeIndexToCoord(activeNodeIndieces(currNodeIndex), 2);
    currNodeX = nodeIndexToCoord(activeNodeIndieces(currNodeIndex), 1);
    matrixIndices = sub2ind(size(pressure), currNodeY, currNodeX);
    pressure(matrixIndices) = pressureAtNodes(currNodeIndex);
    
    %maxabspressure = max(max(abs(pressure)))
    %minpressure = min(min(pressure))
end

% The isoparametric map, computing the Jacobian J at a point (r,s) in the
% reference triangle given the node coordiantes x,y of the physical
% triangle. Note that the current Isopmap deviates from the book, as the
% book contains an error.
function [dSdx,dSdy,detJ] = Isopmap(x,y,r,s)
S=[1-r-s; r; s];
dSdr=[-1; 1; 0];
dSds=[-1; 0; 1];

j11=dot(dSdr,x); j12=dot(dSdr,y);
j21=dot(dSds,x); j22=dot(dSds,y);
detJ=j11*j22-j12*j21;
dSdx=( j22*dSdr-j12*dSds)/detJ;
dSdy=(-j21*dSdr+j11*dSds)/detJ;
end

%A routine tabulating the Gauss quadrature weigths and points on the
%reference triangle, exact solution up to fourth order polynomials. Taken
%from Larson and Bengzon page 214.
function [rspts,qwgts] = Gausspoints(precision)
switch precision
case 1
qwgts=[1];
rspts=[1/3 1/3];
case 2
qwgts=[1/3 1/3 1/3];
rspts=[1/6 1/6;
2/3 1/6;
1/6 2/3];
case 3
qwgts=[-27/48 25/48 25/48 25/48];
rspts=[1/3 1/3;
0.2 0.2;
0.6 0.2;
0.2 0.6];
case 4
qwgts=[0.223381589678011
0.223381589678011
0.223381589678011
0.109951743655322
0.109951743655322
0.109951743655322];
rspts=[0.445948490915965 0.445948490915965;
0.445948490915965 0.108103018168070;
0.108103018168070 0.445948490915965;
0.091576213509771 0.091576213509771;
0.091576213509771 0.816847572980459;
0.816847572980459 0.091576213509771];
otherwise
error('Quadrature precision too high on triangle')
end
end


% Remove boundary nodes from matrix and vector, pressure must be 0
% keep them if they are in obstacles
function [newKx, newKy, newKval, newfx, newfval] = removeBoundaryNodes(...
          Kx, Ky, Kval, fx, fval, nodeIndexToCoord, interfaceNodes)

    fromNodeCoords = nodeIndexToCoord(Kx, :);
    toNodeCoords = nodeIndexToCoord(Ky, :);

    fromNodeInBoundary = interfaceNodes(sub2ind(size(interfaceNodes), ...
        fromNodeCoords(:,2), fromNodeCoords(:,1)));
    toNodeInBoundary = interfaceNodes(sub2ind(size(interfaceNodes), ...
        toNodeCoords(:,2), toNodeCoords(:,1)));
    
    sparseIndiecesToRemove = fromNodeInBoundary | toNodeInBoundary;
        
    Kval(sparseIndiecesToRemove) = [];
    Kx(sparseIndiecesToRemove) = [];
    Ky(sparseIndiecesToRemove) = [];
    
    nodeCoords = nodeIndexToCoord(fx, :);
    inBoundary = interfaceNodes(sub2ind(size(interfaceNodes), ...
        nodeCoords(:,2), nodeCoords(:,1)));
    sparseIndiecesToRemove = inBoundary;
    
    fx(sparseIndiecesToRemove) = [];
    fval(sparseIndiecesToRemove) = [];
    
    newKx = Kx;
    newKy = Ky;
    newKval = Kval;
    newfx = fx;
    newfval = fval;
    
end
