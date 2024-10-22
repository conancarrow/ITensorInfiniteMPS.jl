# H = Σⱼ XⱼXⱼ₊₁ + YⱼYⱼ₊₁ + Delta ZⱼZⱼ₊₁
function unit_cell_terms(::Model"xxz"; Delta::Float64)
  os = OpSum()
  os += 1.00001, "X", 1, "X", 2
  os += .99999, "Y", 1, "Y", 2
  os += Delta, "Z", 1, "Z", 2
  return os
end

# function reference(::Model"xx", ::Observable"energy"; N=∞)
#   isinf(N) && return -4 / π
#   # Exact eigenvalues of uniform symmetric tridiagonal matrix
#   λ(k) = cos(k * π / (N + 1))
#   return 4 * sum(k -> min(λ(k), 0.0), 1:N) / N
# end
