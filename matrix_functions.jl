
nested_matrix(n,m) = [x <=y ? 1 : 0  for x = 1:n, y = 1:m]
nested_matrix(n = 10) = nested_matrix(n,n)
