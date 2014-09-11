module GenSys


export gensys

function new_div(g0, g1)
    small = 1e-6

    F = schurfact(g0+0.0im, g1+0.0im)
    a, b = F[:S], F[:T]

    # TODO: unpack the rest of F like matlab does

    nunstab = 0
    zxz = 0

    for i=1:n
        # ------------------div calc------------
        if !fixdiv
            if abs(a[i, i]) > 0
                divhat = abs(b[i, i])/abs(a[i, i])
                if 1+realsmall < divhat & divhat <= div
                    div = .5 * (1 + divhat)
                end
            end
        end
        # ----------------------------------------
        nunstab=nunstab+(abs(b(i,i))>div*abs(a(i,i)));
        if abs(a(i,i))<realsmall & abs(b(i,i))<realsmall
            zxz=1;
        end
    end

end

function gensys(g0, g1, c, psi, pi, div=1.01)
    eu = [0; 0]
    realsmall = 1e-6
    n = size(g0, 1)

    # I want floating point equality here. This is a stand in check to
    # see if we just have div at its default value
    fixdiv = div == 1.01

    # TODO: this needs to be checked
    a, b, q, z, v = schur(g0, g1)

end  # function

end # module

#=

Test case from Peifan's solution

%Log-linearized model should be written in Sims' form
%Gamma0 * y(t) = Gamma1 * y(t - 1) + C + Psi * z(t) + Pi * eta(t)
%y(t) = [Ex(t+1) xt Epi(t+1) pit it gt ut]
%z(t) = [epsilon_i epsilon_g epsilon_u]
%eta(t) = [err_x err_pi]

sigma = 1.0
kappa = 0.2
beta = 0.99
phi_pi = 1.5
phi_x = 0.5
rho_i = 0.9
sigma_i = 1.0
rho_g = 0.5
rho_u = 0.3
sigma_g = 1.0
sigma_u = 1.0
Gamma0 = [-1  1  -sigma  0  sigma  -1  0;
           0  -kappa  -beta  1  0  0  -1;
           0  -(1-rho_i)*phi_x  0  -(1-rho_i)*phi_pi  1  0  0;
           0  0  0  0  0  1  0;
           0  0  0  0  0  0  1;
           0  1  0  0  0  0  0;
           0  0  0  1  0  0  0];
Gamma1 = [0  0  0  0  0  0  0;
          0  0  0  0  0  0  0;
          0  0  0  0  rho_i  0  0;
          0  0  0  0  0  rho_g  0;
          0  0  0  0  0  0  rho_u;
          1  0  0  0  0  0  0;
          0  0  1  0  0  0  0];
C = [0  0  0  0  0  0  0]';
Psi = [0  0  0;
       0  0  0;
       sigma_i  0  0;
       0  sigma_g  0;
       0  0  sigma_u;
       0  0  0;
       0  0  0];
Pi = [0  0;
      0  0;
      0  0;
      0  0;
      0  0;
      1  0;
      0  1];

%Thanks Sims!
[A, ~, B, ~, ~, ~, ~, eu, ~] = gensys(Gamma0, Gamma1, C, Psi, Pi); %#ok<ASGLU>
%eu = (1, 1) => existence and uniqueness


g0 = matrix(c(-1,  1,  -sigma,  0,  sigma,  -1,  0,
              0,  -kappa,  -beta,  1,  0,  0,  -1,
              0,  -(1-rho_i)*phi_x,  0,  -(1-rho_i)*phi_pi,  1,  0,  0,
              0,  0,  0,  0,  0,  1,  0,
              0,  0,  0,  0,  0,  0,  1,
              0,  1,  0,  0,  0,  0,  0,
              0,  0,  0,  1,  0,  0,  0), nrow=7, ncol=7, byrow=TRUE)
g1 = matrix(c(0,  0,  0,  0,  0,  0,  0,
          0,  0,  0,  0,  0,  0,  0,
          0,  0,  0,  0,  rho_i,  0,  0,
          0,  0,  0,  0,  0,  rho_g,  0,
          0,  0,  0,  0,  0,  0,  rho_u,
          1,  0,  0,  0,  0,  0,  0,
          0,  0,  1,  0,  0,  0,  0), nrow=7, ncol=7, byrow=TRUE)

=#


