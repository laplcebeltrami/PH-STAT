function g=simulate_cirlce_difference(sigma,npoints, nSubjects)
% function g=simulate_cirlce_difference(sigma,npoints, nSubjects)
%
% Simulate topologically different circular patterns with noise level sigma
% Reproduce simulation setting given in 
%
% [2] Chung, M.K., Huang, S.G., Carroll, I.C., Calhoun, V.D. Goldsmith, H.H.
% 2024. Topological state-space estimation of functional human brain networks. 
% PLOS Computational Biology, 20(5), p.e1011869.
%
% (C) 2024 Moo K. Chung
% University of Wisonsin-Madison


nGroup1 = nSubjects;
nGroup2 = nSubjects;
nGroup3 = nSubjects;
nGroup4 = nSubjects;
g1 = cell(nGroup1, 1);
g2 = cell(nGroup2, 1);
g3 = cell(nGroup3, 1);
g4 = cell(nGroup4, 1);

%GROUP 1
%two circles horizontally
for i=1:nGroup1
    circle1 = graph_circle([-1 0],1.5, npoints, sigma);
    circle2 = graph_circle([2 0],1, npoints, sigma);

    %circle1 = graph_arc([-0.5 0],1.5, [0 2*pi], npoints, sigma);
    %circle2 = graph_arc([1 0],1, [0 2*pi], npoints, sigma);

    g1{i}=[circle1; circle2];
end

%GROUP 2
%two circles translated to right
for i=1:nGroup2
    circle1 = graph_circle([-1 0],1.5, npoints, sigma);
    circle2 = graph_circle([1 0],1.5, npoints, sigma);

      %circle1 = graph_arc([-0.5 0],1.5, [0 2*pi], npoints, sigma);
    %circle2 = graph_arc([1 0],1, [pi/15   2*pi-pi/15], npoints, sigma);
    g2{i}=[circle1; circle2];
end

%GROUP 3
%two circles vertically
for i=1:nGroup3
    circle1 = graph_circle([-1 0],1.5, npoints, sigma);
    circle2 = graph_circle([-1 0],1.5, npoints, sigma);

    %circle1 = graph_arc([-0.5 0],1.5, [pi+pi/20 pi-pi/20+2*pi], npoints, sigma);
    %circle2 = graph_arc([1 0],1, [pi/15   2*pi-pi/15], npoints, sigma);
    g3{i}=[circle1; circle2];
end

%GROUP 4
for i=1:nGroup4
    circle1 = graph_circle([-0.5 0],1.5, npoints, sigma);
    circle2 = graph_circle([0.5 0],1.5, npoints, sigma);

    %circle1 = graph_arc([-0.5 0],1.5, [pi/20 2*pi-pi/20], npoints, sigma);
    %circle2 = graph_arc([1 0],1, [pi+pi/20   pi-pi/20+ 2*pi], npoints, sigma);
    g4{i}=[circle1; circle2];
end


% for i=1:nGroup4
%     circle1 = graph_circle([0.5 0],1.5, 60, sigma);
%     circle2 = graph_circle([-0.5 0],1.5, 60, sigma);
%     g4{i}=[circle1; circle2];
% end

% for i=1:nGroup1
%     circle1 = graph_arc([-0.5 0],1.5, [0 2*pi], npoints, sigma)
%     circle2 = graph_arc([1 0],1, [0 2*pi], npoints, sigma)
%     obs =[circle1; circle2];
%     g1{i}=obs;
% end
% 
% for i=1:nGroup2
%     circle1 = graph_arc([-0.5 0],1.5, [0 2*pi], npoints, sigma)
%     circle2 = graph_arc([1 0],1, [pi/15   2*pi-pi/15], npoints, sigma)
%     obs =[circle1; circle2];
%     g2{i}=obs;
% end
% 
% for i=1:nGroup3
%     circle1 = graph_arc([-0.5 0],1.5, [pi+pi/20 pi-pi/20+2*pi], npoints, sigma)
%     circle2 = graph_arc([1 0],1, [pi/15   2*pi-pi/15], npoints, sigma)
%     obs =[circle1; circle2];
%     g3{i}=obs;
% end
% 
% for i=1:nGroup4
%     circle1 = graph_arc([-0.5 0],1.5, [pi/20 2*pi-pi/20], npoints, sigma)
%     circle2 = graph_arc([1 0],1, [pi+pi/20   pi-pi/20+ 2*pi], npoints, sigma)
%     obs =[circle1; circle2];
%     hold on; plot(obs(:,1), obs(:,2),'.k');
%     coord=[obs(:,1)  obs(:,2)];
%     g4{i}=obs;
% end

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