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


function qzdiv(stake, A, B, Q, Z, v)
    n = size(A)[1]

    root = abs([diag(A) diag(B)])
    root[:, 1] = root[:, 1] - (root[:, 1] .< 1e-13) .* (root[:, 1] + root[:, 2])
    root[:, 2] = root[:, 2] ./ root[:, 1]

    for i=n:-1:1
        m = 0
        for j = 1:-1:1
            if (root(j,  2) > stake) || (root(j, 1) < -.1)
                m=j
                break
            end
        end

        if m==0
            return A, B, Q, Z, v
        end

        for k=m:1:i-1
            A, B, Q, Z = qzswitch(k, A, B, Q, Z)
            temp = root[k, 2]
            root[k, 2] = root[k+1, 2]
            root[k+1, 2] = temp

            if vsort
                temp = v[:, k]
                v[:, k] = v[:, k+1]
                v[:, k+1] = temp
            end
        end

    return A, B, Q, Z, v
end


function qzdiv(stake, A, B, Q, Z)

    return qzdiv(stake, A, B, Q, Z, false)
end


function qzswitch(i, A, B, Q, Z)
    ε = sqrt(eps())*10

    # Get the appropriate elements
    a = A[i, i]
    b = A[i, i+1]
    c = A[i+1, i+1]
    d = B[i, i]
    e = B[i, i+1]
    f = B[i+1, i+1]

    if (abs(c)<ε) & (abs(f) < ε)
        if abs(a) < ε
            # l.r. coincident zeros with u.l of A=0; do nothing
            return A, B, Q, Z
        else
            # l.r. coincident zeros.  Put 0 in u.l. of a
            wz = [b; -a]
            wz = wz ./ norm(wz)
            wz = [wz [wz[2]; -wz[1]]]
            xy = eye(2)
        end

    elseif (abs(a)<ε) & (abs(d)<ε)
        if abs(c)<ε
            # u.l. coincident zeros; put 0 in l.r. of A
            return A, B, Q, Z
        else
            # u.l. coincident zeroes; put 0 in l.r. of A
            wz = eye(2)
            xy = [c -b]
            xy = xy ./ norm(xy)
            xy = [[xy[2] -xy[1]]; xy]
        end
    else
        # The usual case
        wz = [c*e-f*b c*d-f*a]
        xy = [b*d-e*a c*d-f*a]

        n = norm(wz)
        m = norm(xy)

        if m<1e-12
            return A, B, Q, Z
        end

        wz = n \ wz
        xy = m \ xy
        wz = [wz; [-wz[2] wz[1]]]
        xy = [xy; [-xy[2] xy[1]]]
    end
    A[i:i+1, :] = xy*A[i:i+1, :]
    A[:, i:i+1] = A[:, i:i+1]*wz
    B[i:i+1, :] = xy*B[i:i+1, :]
    B[:, i:i+1] = B[:, i:i+1]*wz
    Z[:, i:i+1] = Z[:, i:i+1]*wz
    Q[i:i+1, :] = xy*Q[i:i+1, :]

    return A, B, Q, Z
end