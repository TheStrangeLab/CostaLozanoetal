function [f_fdr, p_fdr] = do_fdr(pvec, bonf)
initialCorr = length(pvec);
[psort, indices] = sort(abs(pvec), 'ascend');
f_fdr = [];
p_fdr = 0;
for n = 1:length(psort),
    p0 = psort(n);
    if bonf == 0,
        critp = n * 0.05 / initialCorr;
    else,
        critp = 0.05 / initialCorr;
    end;
    if p0 >= critp,
        break;
    end;
    p_fdr = psort(n);
    f_fdr = [f_fdr; indices(n)];
end;