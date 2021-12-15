using LambertW
using Ranges

a = 0.1 # lower bound
b = 1.3 # upper bound
c = -0.5  # exponential decay

# probability density for τ₁
Z = (1/c) * (exp(c*(b-a)) - 1) - (b - a)
P(τ) = (1/Z) * (exp(c*(b-τ)) - 1)

# cumulative probability
F(τ) = (1/Z) * ((1/c) * (exp(c*(b-a)) - exp(c*(b-τ))) - (τ - a))


function infer_τ(x)
    y = Z*x - (1/c)*exp(c*(b-a)) - a
    A = -(1/c)*exp(c*b)

    arg = max(-1/ℯ, A*c*exp(c*y))
    k = c < 0 ? 0 : -1
    τ = (1/c)*lambertw(arg, k) - y
    
    # Checks
    atol = 1e-10
    @assert isapprox(y, A*exp(-c*τ) - τ; atol=atol)
    @assert isapprox(F(τ), x; atol=atol)
    @assert a-atol <= τ <= b+atol

    return τ
end


# Check with some specific x values

@assert infer_τ(0) ≈ a
@assert a <= infer_τ(0.5) <= b
@assert infer_τ(1) ≈ b


# Check probability distribution
using Random: seed!
seed!(0)

samples = Float64[]
# for i = 1:100_000
for i = 1:10
    x = rand()
    push!(samples, infer_τ(x))
end

using Plots
p = histogram(samples, normalize=true)
taus = range(a, b, 100)
plot!(p, taus, P.(taus))
