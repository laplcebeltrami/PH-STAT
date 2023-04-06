% Statistical Inference on Persistent Homology (PH-STAT)
% 
% The package can be downloaded from
% https://github.com/laplcebeltrami/PH-STAT
%
% (C) 2023 Moo K. Chung, Zijian Chen
% Universtiy of Wisconsin-Madison 
%
% email: mkchung@wisc.edu


%------------------
%% Morse filtration

%The construction for the persistence diagrm for 1D function is based on the 
%iterative pairing and deleation algorithm in
% 
% Chung, M.K., Singh, V., Kim, P.T., Dalton, K.M., Davidson, R.J. 2009. 
% Topological characterization of signal in brain images using the min-max diagram. 
% 12th International Conference on Medical Image Computing and Computer Assisted 
% Intervention (MICCAI). Lecture Notes in Computer Science (LNCS). 5762:158-166.

%Toy example
% In the interval x=[0, 1], we construct a signal s and add noise e to obtain simulated signal Y.
x=[0:0.002:1]';  %interval
s= x + 7*(x - 0.5).^2 + cos(7*pi*x)/2; %signal
e=normrnd(0,0.2,length(x),1); %noise
Y=s+e; %observed data = signal + noise

figure;
plot(x,Y,'ko','MarkerEdgeColor',[0.5 0.5 0.5], 'MarkerFaceColor',[0.7 0.7 0.7], 'MarkerSize',4)

%This produces a scatter points so we can't obtain critical values directly. 
% We smooth out the scatter plot using heat kernel smoothing with bandwidth 0.0001 
% and cosine series expansion upto degree 100. The resulting smoothed signal is stored in the variable wfs.

k=100; sigma=0.0001;
[wfs, beta]=WFS_COS(Y,x,k,sigma);

hold on;
plot(x,wfs,'k','LineWidth',2);

%Then using the iterative pairing and deleation algorithm pairing_1D, we determine 
% the pairing of minimums and maximums using the Elder's rule.  
% pairing_1D.m plots critical values as white notes.

pairs=PH_morse1D(x,wfs);
PH_PD_display(pairs)


% --------------------
%% Rips complex
%genrate random p number of scatter points in [0,1]^d hyper-cube.
%X is the matrix of size p x d of coordinates
rng(2023); %fixed random seeds
p=50; d=3;
X = rand(p, d);

%Compute the Rips complex up to 3-simplices.
radius = 0;
S= PH_rips(X, 2, radius)
%Display Ris complex
PH_rips_display(X,S);
%labels = cellstr(num2str((1:p)', '% d'));
%text(X(:,1)+0.01, X(:,2)+0.01, X(:,3)+0.01, labels, 'Color', 'r', 'FontSize',16) 

% Boundary matrices
B = PH_boundary(S);
betti = PH_boundary_betti(B);
title('\beta_0=1, \beta_1=4, \beta_2=0')


% ---------------
%% Sampling persistent diagram

% Sampling from a single Dirichlet distribution 
alpha=[3,2,3]'; %parameter
weight=1;   %mixing propotion
n=500;

x1 = Dirichlet_sample(alpha,weight,n);
f1 = Dirichlet_density(alpha, weight, x1);
figure; subplot(1,2,1); PH_PD_display(x1, f1); caxis([0 5])

% Sampling from a 3-components mixture of Dirichlet distributions
alpha=[3,8,2;8,3,2;7,3,8;3,7,8]'; %parameter
weight = [0.25,0.25,0.25,0.25]; %mixing propotion
n = 500;

x2 = Dirichlet_sample(alpha,weight,n);
f2 = Dirichlet_density(alpha, weight, x2);
subplot(1,2,2); PH_PD_display(x2, f2); caxis([0 5])


%-------------------
%% Graph filtrations
% The grph filtration is explained in 
%
% Chung, M.K., Lee, H. Ombao. H., Solo, V. 2019 Exact topological inference 
% of the resting-state brain networks in twins, Network Neuroscience 3:674-694
% http://www.stat.wisc.edu/~mchung/papers/chung.2019.NN.pdf
%
% Chung, M.K., Huang, S.-G., Gritsenko, A., Shen, L., Lee, H. 2o19
% Statistical inference on the number of cycles in brain networks. 
% IEEE International Symposium on Biomedical Imaging (ISBI) 113-116
%
% The graph filtration over filtration values epsilon is matheamticall equivaelnt 
% to the Rips filtration with filtration values maximum distance - epsilon up to 1-simplices only.  

rng(2023); %fixed random seeds
p=50; d=3;
X = rand(p, d);

w = pdist2(X,X); %edge weights
maxw = max(w(:)); %maximum edge weight

%Betti curves over graph filtration
thresholds=[0:0.05:maxw]; %filtration values
beta_gf = PH_betti(w, thresholds);

%display betti curves (beta-0 and beta-1)
PH_betti_display(beta_gf,thresholds)

%Equvalently, we can compute the graph filtration as Rips filtraions
%The graph filtration at filtraion value e is equivalent to
e = 0.7;
S= PH_rips(X, 1, maxw-e);
PH_rips_display(X,S);


