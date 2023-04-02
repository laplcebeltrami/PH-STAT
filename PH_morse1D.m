function  pairs = PH_morse1D(x, f)
% FUNCTION  pairs = PH_morse1D(x, f)
% computes the persitent diagram of 1D Morse function f(x) based on 
% the iterative pairing and deletion rule given in
%
%     Chung, M.K., Singh, V., Kim, P.T., Dalton, K.M., Davidson, R.J. 2009.
%     Topological characterization of signal in brain images using the min-max diagram.
%     12th International Conference on Medical Image Computing and Computer Assisted 
%     Intervention (MICCAI). Lecture Notes in Computer Science (LNCS). 5762:158-166.
%     http://www.stat.wisc.edu/~mchung/papers/miccai.2009.pdf 
%
%
% The function is modified from FUNCTION pairs = pairing_1D(x, f)
%
% (C) Moo K. Chung 2009
%     mkchung@wisc.edu
%     University of Wisconsin-Madison
%
% The code is downloaded from
% https://github.com/laplcebeltrami/PH-STAT

%---------------------
% construct a linear graph establishing connectivity between data points

n_vertex=length(x);
n_tri=n_vertex-1;

nbr=zeros(n_vertex,2);
for i=2:n_vertex-1
    nbr(i,:)=[i-1, i+1];
end;

nbr(1,1)=2;
nbr(n_vertex,1)=n_vertex-1;

degree=2*ones(n_vertex,1);
degree(1)=1;
degree(n_vertex)=1;

%--------------------
% fina all local minimum and maximum.

maxind=[];
minind=[];
for i=2:n_vertex-1
    flag=sum(f(i) >= f(nbr(i,1:degree(i))));
    % max
    if flag== degree(i)
        maxind=[maxind i];
        % min
    elseif flag ==0
        minind=[minind i];
    end;
end;
   

%----------------------------------------
% reconstruct graph using min and max only


xmin=x(minind);
xmax=x(maxind);
ff=f([minind maxind]);
xx=[xmin ; xmax];

n_vertex=length(xx);

xsort=sort(xx);
fsort=[];
for i=1:n_vertex
    fsort(i)=ff(find(xx==xsort(i)));
end;

hold on; plot(xsort,fsort,'o','MarkerEdgeColor','r', 'MarkerFaceColor','w', 'MarkerSize',7)
% xsort, fsort    contains min and max only

%---------------------
n_vertex=length(xsort);
n_tri=n_vertex-1;

nbr=zeros(n_vertex,2);
for i=2:n_vertex-1
    nbr(i,:)=[i-1, i+1];
end;
nbr(1,1)=2;
nbr(n_vertex,1)=n_vertex-1;

degree=2*ones(n_vertex,1);
degree(1)=1;
degree(n_vertex)=1;

%--------------------------------
% g: minimum
% h: maximum

%coordmin=x(minind);
%coordmax=x(maxind);


g =f(minind);
h= f(maxind);

gsort=sort(g);
hsort=sort(h);

m=length(g);
n=length(h);

hset=[];
for i=1:n
    hset=[hset find(fsort==hsort(i))];
end;
pairs=[];
for i=m:-1:1
    ind=find(fsort==gsort(i));
    nbr_ind=nbr(ind,1:degree(ind));
    nbr_ind=intersect(nbr_ind, hset);

    if ~isempty(nbr_ind)

        nbr_value =fsort(nbr_ind);
        argmin=find(nbr_value>gsort(i));
        if ~isempty(argmin)
            hstar = min(nbr_value(argmin));
            hstar_ind=find(fsort==hstar);
            hset=setdiff(hset,hstar_ind);
            pairs=[pairs; [gsort(i) hstar]];
        end;
    end;
end;

