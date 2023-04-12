function [transStat, elapsedTimes] = WS_transpositions(lossMtx,nGroup1, nGroup2, nTrans, permNo)
% function [transStat, elapsedTimes] = transposition_test(lossMtx,nGroup1, nGroup2, nTrans, permNo)
%
% This function performs the transposition procedure by transposing group labels and
% summing pair-wise distances from a pre-computed loss matrix. The transposition test 
% for ratio statistic is given in 
%
% [1] Songdechakraiwut, T., Shen, L., Chung, M.K. 2021 Topological learning and 
%its application to multimodal brain network integration, Medical Image 
%Computing and Computer Assisted Intervention (MICCAI), LNCS 12902:166-176 
%
% [2] Songdechakraiwut, T. Chung, M.K. 2023 Topological learning for brain networks, 
% Annals of Applied Statistics 17:403-433, arXiv:2012.00675.
% https://arxiv.org/pdf/2012.00675.pdf
% 
%
% INPUT
% lossMtx  : loss matrix whose entries are all possible pair-wise distances
% nGroup1  : sample size of group 1
% nGroup2  : sample size of group 2
% nTrans   : the number of transpositions being performed
% permNo   : intermix random permutation every 'permNo' transpositions
%
%
% OUTPUT
% transStat    : array of ratio statistic, 1 x nTrans
% elapsedTimes : elapsed time, 1 x nTrans
%
% %If you are using any part of the code, please reference one of the above paper.
%
% (C) 2020 Tananun Songdechakraiwut, Moo K. Chung
%          University of Wisconsin-Madison
%
%  Contact tananun@cs.wisc.edu or mkchung@wisc.edu
%  for support/permission with the codes 
%
% Update history
%     2020 November 11, created - Tananun Songdechakraiwut
%     2020 December 10, edited Moo Chung
%    

tStart = tic;
totalNo = nGroup1 + nGroup2;
% denomGroup1 = (nGroup1*(nGroup1-1) + nGroup2*(nGroup2-1))/2;
% denomGroup2 = nGroup1*nGroup2;
% constant = denomGroup1 / denomGroup2;

transStat = zeros(1, nTrans);
elapsedTimes = zeros(1, nTrans);

for t = 1:round(nTrans/permNo)
    startInd = (t - 1) * permNo + 1;

    %% Initialize sum of pair-wise distances for both within and between groups
    permutation = randperm(totalNo);

    % dividing group labels into two groups
    permutedG1 = permutation(1:nGroup1);
    permutedG2 = permutation(nGroup1+1:totalNo);

    % within groups
    within = 0;
    % sum of pair-wise distances within groups
    for i = 1:nGroup1 % group 1
        for j = i + 1:nGroup1
            within = within + lossMtx(permutedG1(i), permutedG1(j));
        end
    end
    for i = 1:nGroup2 % group 2
        for j = i + 1:nGroup2
            within = within + lossMtx(permutedG2(i), permutedG2(j));
        end
    end

    % between groups
    % sum of pair-wise distances between groups
    between = 0;
    for i = 1:nGroup1
        for j = 1:nGroup2
            between = between + lossMtx(permutedG1(i), permutedG2(j));
        end
    end

    % saving the ratio statistics
    transStat(startInd) = between/within; % * constant;
    elapsedTimes(startInd) = toc(tStart);

    %% transposition procedure

    prevWithin = within;
    prevBetween = between;
    for n = startInd + 1:startInd + permNo - 1
        % random transposition indices
        ind1 = randi(nGroup1);
        ind2 = randi(nGroup2);

        % within groups
        within = 0;
        % sum of offset distances being removed from within-group distances
        for i = 1:nGroup1 % group 1
            within = within + lossMtx(permutedG1(i), permutedG1(ind1));
        end
        for i = 1:nGroup2 % group 2
            within = within + lossMtx(permutedG2(i), permutedG2(ind2));
        end
        % remove excess distances to itself
        within = within - lossMtx(permutedG1(ind1), permutedG1(ind1));
        within = within - lossMtx(permutedG2(ind2), permutedG2(ind2));

        % between groups
        % sum of offset distances being removed from between-group distances
        between = 0;
        for i = 1:nGroup2
            between = between + lossMtx(permutedG1(ind1), permutedG2(i));
        end
        for i = 1:nGroup1
            between = between + lossMtx(permutedG1(i), permutedG2(ind2));
        end
        between = between - 2 * lossMtx(permutedG1(ind1), permutedG2(ind2));

        % update iteratively
        delta = between - within;
        prevWithin = prevWithin + delta;
        prevBetween = prevBetween - delta;

        % saving the ratio statistics
        transStat(n) = prevBetween/prevWithin; % * constant;
        elapsedTimes(n) = toc(tStart);

        % swapping networks corresponding to indices
        temp = permutedG1(ind1);
        permutedG1(ind1) = permutedG2(ind2);
        permutedG2(ind2) = temp;
    end
end
