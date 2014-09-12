#=

    Test case from Peifan's solution to Spring 2014 HW

    Log-linearized model should be written in Sims' form

    Γ0 * y(t) = Γ1 * y(t - 1) + C + Ψ * z(t) + Π * eta(t)

    y(t) = [Ex(t+1) xt Epi(t+1) pit it gt ut]

    z(t) = [epsilon_i epsilon_g epsilon_u]

    eta(t) = [err_x err_pi]
=#

using QuantEcon.GenSys
srand(42)  # reproducible!

σ = 1.0
κ = 0.2
β = 0.99
ϕ_π = 1.5
ϕ_x = 0.5
ρ_i = 0.9
σ_i = 1.0
ρ_g = 0.5
ρ_u = 0.3
σ_g = 1.0
σ_u = 1.0
Γ0 = [-1       1       -σ       0      σ  -1  0
      0        -κ      -β       1      0  0  -1
      0  -(1-ρ_i)*ϕ_x  0  -(1-ρ_i)*ϕ_π 1  0  0
      0        0       0        0      0  1  0
      0        0       0        0      0  0  1
      0        1       0        0      0  0  0
      0        0       0        1      0  0  0]
Γ1 = [0  0  0  0   0    0     0
      0  0  0  0   0    0     0
      0  0  0  0  ρ_i   0     0
      0  0  0  0   0   ρ_g    0
      0  0  0  0   0    0    ρ_u
      1  0  0  0   0    0     0
      0  0  1  0   0    0     0]
c = zeros(size(Γ0, 1), 1)

ψ = [0    0   0
     0    0   0
     σ_i  0   0
     0   σ_g  0
     0    0  σ_u
     0    0   0
     0    0   0]
π = [0  0
     0  0
     0  0
     0  0
     0  0
     1  0
     0  1]

G1, C, impact, fmat, fwt, ywt, gev, eu, loose  = gensys(Γ0, Γ1, c, ψ, π)

A, B = G1, impact

# Generate choice function
C = [0 1 0 0 0 0 0
     0 0 0 1 0 0 0
     0 0 0 0 1 0 0]

## -------- ##
#- Simulate -#
## -------- ##
st = zeros(7)
burn = 1000
T = 10000 + burn
obs = zeros(3, T)
for i=1:T
    shock = randn(3)
    st = A * st + B * shock
    obs[:, i] = C * st
end
obs = obs[:, burn+1:end]

# eu = (1, 1) => existence and uniqueness
