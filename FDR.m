function pID = FDR(p,q)
% FORMAT pID = FDR(p,q)
%
% INPUT
%   p : vector of uncorrected p-values.
%   q : desired False Discovery Rate (FDR) level.
%
% OUTPUT
%   pID : Benjamini-Hochberg FDR p-value threshold.
%
% (C) Moo K. Chung
% University of Wisconsin-Madison

p = p(isfinite(p));
p = sort(p(:));

V = length(p);
I = (1:V)';

idxID = find(p <= I/V*q);

if isempty(idxID)
    pID = 0;
else
    pID = p(max(idxID));
end

fprintf('Number of tests = %d\n', V)
fprintf('Requested FDR level q = %.4f\n', q)
fprintf('BH FDR threshold = %.4f\n', pID)
fprintf('Number significant tests = %d\n', sum(p <= pID))

m = V;
threshold_line005 = (1:m) * 0.05 / m;
threshold_line001 = (1:m) * 0.01 / m;

figure;
plot(p,'ok')
hold on
plot(threshold_line005,'k','LineWidth',2)
plot(threshold_line001,'k-.','LineWidth',2)
text(round(0.85*m),0.06,'q=0.05','FontSize',16)
text(round(0.85*m),0.012,'q=0.01','FontSize',16)
ylabel('p-values')
figure_bigger(16)

end


