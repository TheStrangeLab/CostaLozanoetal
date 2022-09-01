function [f_fdr, p_fdr] = teg_do_fdr(pvec, fdr_crit)

f_fdr = [];
p_fdr = 0;
[p_sort, indices] = sort(abs(pvec), 'ascend');
N = length(p_sort);

for iP = 1:length(p_sort),
    p_try = p_sort(iP);
    [n_sig, n_exp] = inner_teg_do_fdr(p_sort, p_try, N);
    if n_sig * fdr_crit > n_exp,
        f_fdr = find(pvec <= p_try);
        p_fdr = p_try;
    end;
end;

function [n_sig, n_exp] = inner_teg_do_fdr(p_sort, p_try, N)
n_sig = length(find(p_sort <= p_try));
n_exp = p_try * N;
