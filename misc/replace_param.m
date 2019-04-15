function x_vect_out = replace_param(x_vect, indices_to_replace, x)
x_vect_out = x_vect;
for ind = 1:length(indices_to_replace)
    x_vect_out(indices_to_replace(ind)) = x(ind);
end
