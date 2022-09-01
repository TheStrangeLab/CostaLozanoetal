function teg_vv(M)

for n = 1:size(M, 1),
    fprintf([num2str([n M(n, :)]) '\n']);
end;
