function S = above_zero_alt(S,name)

pass_zero = find(S(:,4)<0);
if size(pass_zero,1)>0
    land = (pass_zero(1)-1);
    S = S(1:land,:);
    landing_time = S(end,1);
    str = strcat(name, ' landing time: %i\n');
    fprintf(str,landing_time);
end