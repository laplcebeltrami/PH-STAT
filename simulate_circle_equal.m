function g=simulate_cirlce_equal(sigma, npoints, nSubjects)
% g=simulate_cirlce_equal(sigma, npoints)
%simulate topologically equivalent circular patterns with noise level sigma
%
% (C) 2021- Moo K. Chung
% Univeristy of Wisconsin-Madison

nGroup1 = nSubjects;
nGroup2 = nSubjects;
nGroup3 = nSubjects;
nGroup4 = nSubjects;
g1 = cell(nGroup1, 1);
g2 = cell(nGroup2, 1);
g3 = cell(nGroup3, 1);
g4 = cell(nGroup4, 1);


for i=1:nGroup1
    circle1 = graph_circle([1.5 0],1, npoints, sigma);
    circle2 = graph_circle([-1 0],1.5, npoints, sigma);

    obs =[circle1; circle2];
    %figure; imagesc(C); colorbar
    g1{i}=obs;
end

for i=1:nGroup2
    circle1 = graph_circle([-1.5 0],1, npoints, sigma);
    circle2 = graph_circle([1 0],1.5, npoints, sigma);
    obs =[circle1; circle2];
    g2{i}=obs;
end

for i=1:nGroup3
    circle1 = graph_circle([0 2],1, 60, sigma);
    circle2 = graph_circle([0 -1],1.5, 60, sigma);
    obs =[circle1; circle2];
    g3{i}=obs;
end

for i=1:nGroup4
    circle1 = graph_circle([0 1],1.5, 60, sigma);
    circle2 = graph_circle([0 -2],1, 60, sigma);
    obs =[circle1; circle2];
    g4{i}=obs;
end

g=[g1 g2 g3 g4];


% helper function
function coord = graph_circle(center,radius, n, sigma)
%function coord = graph_circle(center,radius)
%
%Draw a circle at center with radius

theta = linspace(0,2*pi,n); %sample n points along circle

x=radius*(cos(theta)+center(1));
y=radius*(sin(theta)+center(2));
coord=[x ; y]'; %2 x n matrix

coord=coord + normrnd(0,sigma,n,2);

function coord = graph_arc(center,radius, arc, n, sigma)
%function coord = graph_arc(center,radius, arc, n, sigma)
%
%Draw a circle at center with radius between arc


theta = linspace(arc(1),arc(2),n); %sample n points along circle

x=radius*(cos(theta)+center(1));
y=radius*(sin(theta)+center(2));
coord=[x ; y]'; %2 x n matrix

coord=coord + normrnd(0,sigma,n,2);