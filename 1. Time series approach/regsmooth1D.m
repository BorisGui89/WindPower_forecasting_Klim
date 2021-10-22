function [x,yTilde] = regsmooth1D(data, points, order ,h);
% Locally weighted polynomial regression
% Tri-cube weight-function
% 
% [x,yTilde] = regsmooth1D(data, points, order ,h);
% 
% Input:
% Data:   [x ; y]
% ponits: Grid resolution on x axis
% order:  order of fitted polynomium
% h: bandwidth 0<h<1
%
% Stig B. Mortensen, 2004.




% number of nearest neighbors
NN = round(h*length(data));

% create grid of x,y for interpolation
% scale x and y axes to [-1,1]
minX = min(data(:,1)); maxX = max(data(:,1)); scaleX = maxX-minX;
dataSc(:,1) = ((data(:,1)-minX)/scaleX)*2-1; %IMPORTANT: scale data to [-1,1]

%create grid in scaled space
xVec = linspace(-1,1,points);
x = xVec;


for p = 1:points
    % find distance from current point x,y to data
    dist = abs(dataSc(:,1)-x(p));
    [sortDist sortIdx] = sort(dist);
    
    %extract nearest neighbors points
    sortIdx = sortIdx(1:NN);
    distNN = dist(sortIdx); 
    xNN = dataSc(sortIdx,1);
    zNN = data(sortIdx,2);
    
    % designmatrice and evaluation point
    U = ones(NN,1);
    ep = 1;
    for i = 1:order
        U = [U xNN.^i];
        ep = [ep x(p)^i];
    end
    
    % weights tri-cube
    distNN = distNN/max(distNN); %truncate distNN to [0,1] for proper wheigths
    W = diag((1-distNN.^3).^3); 
    
    % solve normal equation
    % theta = (U'*W*U)^(-1)*U'*W*zNN;
    theta = (U'*W*U)\(U'*W*zNN);
    
    % evaluate estimated locally wheighted polynomial
    yTilde(p) = ep*theta; %theta(1); %ep*theta;
    
    
end

% replace x by unscaled axes
x = linspace(minX,maxX,points);


