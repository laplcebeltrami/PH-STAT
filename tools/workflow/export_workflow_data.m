function export_workflow_data(repo, output)
% Export deterministic PH-STAT demonstration data for the README animation.
% Run with the pinned upstream toolbox on the MATLAB path.
addpath(repo, '-begin');
previous_folder = pwd;
restore_folder = onCleanup(@() cd(previous_folder));
cd(repo); % Current folder takes precedence over MATLAB's search path.
set(groot, 'defaultFigureVisible', 'off');
rng(20261016, 'twister');
data.source_commit = '835d5eecc760181ecef06883afbbb42c24bcfb0e';
data.matlab_version = version;
data.kind = 'Synthetic examples, adapted from SCRIPT.mlx; not clinical results';

% Same signal, smoothing settings, and pairing routine as SCRIPT.mlx.
x = (0:0.002:1)';
s = x + 7*(x - 0.5).^2 + cos(7*pi*x)/2;
Y = s + normrnd(0,0.2,length(x),1);
[wfs, ~] = WFS_COS(Y,x,100,0.0001);
pairs = PH_morse1D(x,wfs);
data.morse = struct('x',x,'observed',Y,'smooth',wfs,'pairs',pairs);
fprintf('Morse: %d finite pairs\n',size(pairs,1));
close all;

% Small ring example makes loop birth and filling readable in a GIF.
rng(20261016, 'twister');
theta = (0:15)'*2*pi/16;
X = [cos(theta), sin(theta)] + 0.04*randn(16,2);
radii = linspace(0.05,1.9,61);
betti = zeros(length(radii),2);
for i = 1:length(radii)
    S = PH_rips(X,2,radii(i));
    B = PH_boundary(S);
    betti(i,:) = PH_boundary_betti(B);
end
S = PH_rips(X,2,max(radii));
data.rips = struct('points',X,'radii',radii,'betti',betti, ...
    'edges',S{2,1},'edge_birth',S{2,2}, ...
    'faces',S{3,1},'face_birth',S{3,2});
fprintf('Rips: beta1 range %d to %d\n',min(betti(:,2)),max(betti(:,2)));

% PH_betti keeps weights strictly greater than threshold: edges disappear.
coords = [-1.0 0.6;-0.7 1.1;-0.3 0.7;-0.7 0.1; ...
           0.5 0.8;0.9 1.1;1.2 0.6;0.8 0.2; ...
          -0.2 -0.6;0.3 -0.6;0.5 -1.0;-0.3 -1.1];
C = exp(-pdist2(coords,coords).^2/1.2);
C = C + 0.015*rand(size(C)); C = (C+C')/2;
C(1:size(C,1)+1:end) = 0;
thresholds = linspace(0,1,61);
beta = PH_betti(C,thresholds);
data.graph = struct('coords',coords,'weights',C,'thresholds',thresholds, ...
    'beta0',beta.zero,'beta1',beta.one);
assert(all(diff(beta.zero)>=0) && all(diff(beta.one)<=0));

% The 4-versus-6, 50-node simulation in the current SCRIPT.mlx.
p=50; nA=4; nB=6;
rng(1); signal=randcorr(p);
conA=zeros(p,p,nA); conB=zeros(p,p,nB);
for i=1:nA
    conA(:,:,i)=signal+randcorr(p);
end
rng(2);
for j=1:nB
    conB(:,:,j)=signal+2*randcorr(p);
end
D=WS_pdist2(conA,conB);
observed=WS_ratio(D.D01,nA,nB);
% Enumerate every allocation for a reproducible exact small-sample display.
allocations=nchoosek(1:nA+nB,nA);
null=zeros(size(allocations,1),1);
for i=1:size(allocations,1)
    order=[allocations(i,:) setdiff(1:nA+nB,allocations(i,:))];
    null(i)=WS_ratio(D.D01(order,order),nA,nB);
end
p_exact=mean(null>=observed-1e-12);
rng(20261016);
[transStat,~]=WS_transpositions(D.D01,nA,nB,10000,100);
data.inference = struct('meanA',mean(conA,3),'meanB',mean(conB,3), ...
    'D01',D.D01,'observed',observed,'null',null,'p_exact',p_exact, ...
    'nA',nA,'nB',nB,'n_labelings',size(allocations,1), ...
    'p_transposition',mean(transStat>=observed));
assert(norm(D.D01-D.D01','fro')<1e-10);
assert(norm(D.D01-(D.D0+D.D1),'fro')<1e-8);
fprintf('Inference: R=%.6f, exact p=%.6f (%d allocations)\n',observed,p_exact,length(null));

% Four groups from the current live-script simulation, reduced for display.
rng(20261016);
clouds = simulate_circle_difference(0.18,24,8);
G=cell(size(clouds));
allC=zeros(48,48,numel(G));
features=[];
for j=1:size(clouds,2)
    for i=1:size(clouds,1)
        q=i+(j-1)*size(clouds,1);
        G{i,j}=pdist2(clouds{i,j},clouds{i,j});
        allC(:,:,q)=G{i,j};
        [b,d]=WS_decompose(G{i,j});
        features(q,:)=[b(:,3);d(:,3)]';
    end
end
Z=WS_embed(allC);
acc=WS_cluster(G);
rng(20261016);
labels=kmeans(features,4,'Replicates',20);
data.clustering = struct('xy',[Z.x'+Z.center(1),Z.y'+Z.center(2)], ...
    'groups',repelem((1:4)',8),'labels',labels,'clouds',{clouds(1,:)}, ...
    'accuracy_100_starts',acc,'n_per_group',8);
fprintf('Clustering: %d networks; 100-start mean accuracy %.4f\n',numel(G),acc.mean);

% The five-node, filled-triangle example in SCRIPT.mlx.
Y=[0 1 1 0 0; -1 0 2 1.5 0; -1 -2 0 0 1.5; ...
    0 -1.5 0 0 0.5; 0 0 -1.5 -0.5 0];
S=Hodge_2Skeleton(Y);
B=PH_boundary(S);
X=Hodge_vec(Y);
% Verify edge indexing before associating values with drawn arrows.
expected=Y(sub2ind(size(Y),S{2}(:,1),S{2}(:,2)));
assert(norm(X-expected)<1e-12);
[Xg,Xc,Xh]=Hodge_decompose(X,B);
Lgraph=Hodge_laplacian(B);
L=full(Lgraph+B{2}*B{2}');
[Phi,Lambda]=eig(L); lambda=diag(Lambda);
t=linspace(0,5,81);
flows=Phi*(exp(-lambda*t).*(Phi'*X));
err=[norm(X-Xg-Xc-Xh),norm(B{1}*Xh),norm(B{2}'*Xh), ...
    abs(Xg'*Xc),abs(Xg'*Xh),abs(Xc'*Xh)];
assert(max(err)<1e-8);
data.hodge=struct('coords',[0 1.2;-0.8 0.45;0.8 0.45;-0.8 -0.8;0.8 -0.8], ...
    'edges',S{2},'faces',S{3},'input',X,'gradient',Xg,'curl',Xc, ...
    'harmonic',Xh,'time',t,'diffusion',flows,'validation_errors',err);
fprintf('Hodge: maximum reconstruction/orthogonality error %.3g\n',max(err));

fid=fopen(output,'w');
fwrite(fid,jsonencode(data),'char');
fclose(fid);
close all;
fprintf('Exported %s\n',output);
end
