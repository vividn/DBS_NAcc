function outvec = unitvec(vector)
outvec = vector ./ sqrt(sum(vector.^2));
end