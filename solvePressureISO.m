function [ pressure ] = solvePressureTestA( xCells,yCells, velocity, distanceField, obstacle )
%SOLVEPRESSURE Apply pressure solver with finite elements method
%   return pressure: matrix(yNodes, xNodes)
%   simulationSettings: struct with these fields:
%       cellSize: grid cell size in m
%   velocity: matrix(yNodes, xNodes, dim)
%   distanceField: matrix(yNodes, xNodes)
%   obstacle: boolean matrix(yNodes, xNodes)

    %% initialization
    % using y,x only when accessing matrizes from caller, internal is x,y
    numCellsX = xCells;
    numCellsY = yCells;
    
    velocityFieldX = velocity(:,:,1);
    velocityFieldY = velocity(:,:,2);

    % sparse matrix constructor parameters
    Kx = [];
    Ky = [];
    Kval = [];
    
    % pressure vector as sparse matrix
    fx = [];
    fval = [];
    
    % get fluid cells and interface nodes
    fluidCells = getFluidCells(distanceField);
    interfaceNodes = getInterfaceNodes(distanceField, obstacle, true);
    
    nodeIndex = 0;
    nodeCoordToIndex = zeros(size(distanceField, 1), size(distanceField, 2));
    nodeIndexToCoord = zeros(size(distanceField, 1) * size(distanceField, 2), 2);
    nodesInElements = zeros(3,size(distanceField, 1) * size(distanceField, 2));
    nodeOffsets = [0 0; 0 1; 1 1; 1 0];
    element = 1;
    

    
    
    %% pressure solving
    % for each cell
    for cellY = 1 : numCellsY 
        for cellX = 1 : numCellsX
            % surrounded by 4 nodes
            % skip cells without water
            if ~fluidCells(cellY, cellX)
                continue;
            end
            
            currentCoord = [cellX cellY];
            
            wTildeEX = [0; 0; 0; 0];
            wTildeEY = [0; 0; 0; 0];
          
           
          
           
            % gather all information for cell
            for corner = 1 : 4
                cornerCoord = currentCoord + nodeOffsets(corner, :);
                
                wTildeEX(corner) = velocityFieldX(cornerCoord(2), cornerCoord(1));
                wTildeEY(corner) = velocityFieldY(cornerCoord(2), cornerCoord(1));
                
                % add this node to indexing
                if nodeCoordToIndex(cornerCoord(2), cornerCoord(1)) == 0
                    nodeIndex = nodeIndex + 1;
                    
                    nodeCoordToIndex(cornerCoord(2), cornerCoord(1)) = nodeIndex;
                    nodeIndexToCoord(nodeIndex, :) = cornerCoord;
                end
            end
            nodeCoord = repmat(currentCoord, [4,1])+ nodeOffsets;
                nodesInElements(:,element) =...
                   [nodeCoordToIndex(nodeCoord(1,2),nodeCoord(1,1)),... 
                   nodeCoordToIndex(nodeCoord(2,2),nodeCoord(2,1)),...
                   nodeCoordToIndex(nodeCoord(4,2),nodeCoord(4,1))];
                nodesInElements(:,element+1) =...
                   [nodeCoordToIndex(nodeCoord(2,2),nodeCoord(2,1)),... 
                   nodeCoordToIndex(nodeCoord(3,2),nodeCoord(3,1)),...
                   nodeCoordToIndex(nodeCoord(4,2),nodeCoord(4,1))];
               
              



    

for t=1:2
    %Stiffnessassembler ISO
[rspts,qwgts]=Gausspoints(2); % quadrature rule

loc2glb = nodesInElements(1:3,(element+t-1)); % node numbers
x = nodeIndexToCoord(loc2glb,2);
y = nodeIndexToCoord(loc2glb,1);
AK=zeros(3,3); % elements stiffness
for q=1:length(qwgts) % quadrature loop
r=rspts(q,1); % quadrature r-coordinate
s=rspts(q,2); % s-
[S,dSdx,dSdy,detJ,dSdr,dSds]=Isopmap(x,y,r,s);
wxarea=qwgts(q)*detJ/2; % weight times area
AK=AK+(dSdx*dSdx'...
+dSdy*dSdy')*wxarea; % element stiffness
end
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
                
                newFx(fromCorner) = fromCornerNodeIndex;
            end
            Kx = [Kx; newKx];%adding newKx to the groing columnvector Kx
            Ky = [Ky; newKy];
            Kval = [Kval; newKval];
            fx = [fx; newFx];



    %Load assembler
    bK=0;
%loc2glb = nodesInElements(1:3,(element+t-1));
%x = nodeIndexToCoord(loc2glb,2);
%y = nodeIndexToCoord(loc2glb,1);
wTildeEX;
test=dSdx;
test2=dSdr;
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


function [S,dSdx,dSdy,detJ,dSdr,dSds] = Isopmap(x,y,r,s)
S=[1-r-s; r; s];
dSdr=[-1; 1; 0];
dSds=[-1; 0; 1];

j11=dot(dSdr,x); j12=dot(dSdr,y);
j21=dot(dSds,x); j22=dot(dSds,y);
detJ=j11*j22-j12*j21;
dSdx=( j22*dSdr-j12*dSds)/detJ;
dSdy=(-j21*dSdr+j11*dSds)/detJ;
end

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

function [newKx, newKy, newKval, newfx, newfval] = removeBoundaryNodes(...
          Kx, Ky, Kval, fx, fval, nodeIndexToCoord, interfaceNodes)

    % remove boundary nodes from matrix and vector, pressure must be 0
    % keep them if they are in obstacles
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
