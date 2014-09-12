module GenSys


export gensys

function new_div(F::Base.LinAlg.GeneralizedSchur)
    ε = 1e-6  # small number to check convergence
    n = size(F[:T], 1)

    a, b, q, z = real(F[:S]), real(F[:T]), F[:Q]', F[:Z]

    nunstab = 0.0
    zxz = 0

    div = 1.01

    for i=1:n
        if abs(a[i, i]) > 0
            divhat = abs(b[i, i]) / abs(a[i, i])
            if 1 + ε < divhat && divhat <= div
                div = .5 * (1 + divhat)
            end
        end
    end

    return div
end


## -------------- ##
#- gensys methods -#
## -------------- ##

# method if no div is given
function gensys(g0, g1, c, ψ, π)
    # TODO: right now a, b, q, and z *exactly* match MATLAB -- hence
    #       the somewhat ugly code. We should check to see if we actually
    #       need to do the converting to complex before computing schurfact
    zg0, zg1 = g0 + 0.0im, g1 + 0.0im
    F = schurfact(zg0, zg1)
    div = new_div(F)
    gensys(F, c, ψ, π, div)
end


# method if all arguments are given
function gensys(g0, g1, c, ψ, π, div)
    zg0, zg1 = g0 + 0.0im, g1 + 0.0im
    F = schurfact(zg0, zg1)
    gensys(F, c, ψ, π, div)
end


# Method that does the real work. Work directly on the decomposition
function gensys(F::Base.LinAlg.GeneralizedSchur, c, ψ, π, div)
    eu = [0, 0]
    ε = 1e-6  # small number to check convergence
    nunstab = 0.0
    zxz = 0
    a, b, q, z = map(real, {F[:S], F[:T], F[:Q]', F[:Z]})
    n = size(a, 1)

    for i=1:n
        nunstab += (abs(b[i, i]) > div * abs(a[i,i]))
        if abs(a[i, i]) < ε && abs(b[i, i]) < ε
            zxz = 1
        end
    end

    if zxz == 1
        msg = "Coincident zeros.  Indeterminacy and/or nonexistence."
        throw(InexactError(msg))
    end

    a, b, q, z = qzdiv(div, a, b, q, z)
    gev = [diag(a) diag(b)]

    q1 = q[1:n-nunstab, :]
    q2 = q[n-nunstab+1:n, ]
    z1 = z[:, 1:n-nunstab]'
    z2 = z[:, n-nunstab+1:n]'
    a2 = a[n-nunstab+1:n, n-nunstab+1:n]
    b2 = b[n-nunstab+1:n, n-nunstab+1:n]
    etawt = q2 * pi
    neta = size(pi, 2)

    # branch below is to handle case of no stable roots, which previously
    # (5/9/09) quit with an error in that case.
    if isapprox(nunstab, 0.0)
        etawt == zeros(0, neta)
        ueta = zeros(0, 0)
        deta = zeros(0, 0)
        veta = zeros(neta, 0)
        bigev = 0
    else
        ueta, deta, veta = svd(etawt)
        deta = diagm(deta)  # TODO: do we need to do this
        md = min(size(deta))
        bigev = find(diag(deta[1:md,1:md]) .> ε)
        ueta = ueta[:, bigev]
        veta = veta[:, bigev]
        deta = deta[bigev, bigev]
    end

    eu[1] = length(bigev) >= nunstab  # PERF: not type stable

    # ----------------------------------------------------
    # Note that existence and uniqueness are not just matters of comparing
    # numbers of roots and numbers of endogenous errors.  These counts are
    # reported below because usually they point to the source of the problem.
    # ------------------------------------------------------

    # branch below to handle case of no stable roots
    if nunstab == n
        etawt1 = zeros(0, neta)
        bigev = 0
        ueta1 = zeros(0, 0)
        veta1 = zeros(neta, 0)
        deta1 = zeros(0, 0)
    else
        etawt1 = q1 * pi
        ndeta1 = min(n-nunstab, neta)
        ueta1, deta1, veta1 = svd(etawt1)
        deta1 = diagm(deta1)  # TODO: do we need to do this
        md = min(size(deta1))
        bigev = find(diag(deta1[1:md, 1:md]) .> realsmall)
        ueta1 = ueta1[:, bigev]
        veta1 = veta1[:, bigev]
        deta1 = deta1[bigev, bigev]
    end

    if isempty(veta1)
        unique = 1
    else
        loose = veta1-veta*veta'*veta1
        ul, dl, vl = svd(loose)
        dl = diagm(dl)  # TODO: do we need to do this
        nloose = sum(abs(diag(dl)) .> realsmall*n)
        unique = (nloose == 0)
    end

    if unique
        info("gensys: Unique solution!")
        eu[2] = 1
    else
        println("Indeterminacy. $nloose loose endog errors.")
    end

    tmat = [eye(n-nunstab) -(ueta*(deta\veta')*veta1*deta1*ueta1')']
    G0 = [tmat*a; zeros(nunstab,n-nunstab) eye(nunstab)]
    G1 = [tmat*b; zeros(nunstab,n)]
    # ----------------------
    # G0 is always non-singular because by construction there are no zeros on
    # the diagonal of a(1:n-nunstab,1:n-nunstab), which forms G0's ul corner.
    # -----------------------
    G0I = inv(G0)
    G1 = G0I*G1
    usix = n-nunstab+1:n
    C = G0I * [tmat*q*c; (a[usix, usix] - b[usix,usix])\q2*c]
    impact = G0I * [tmat*q*psi; zeros(nunstab, size(psi, 2))]
    fmat = b(usix, usix)\a[usix,usix]
    fwt = -b[usix, usix]\q2*psi
    ywt = G0I[:, usix]

    # Correction 5/07/2009:  formerly had forgotten to premultiply by G0I
    loose = G0I * [etawt1 * (eye(neta) - veta * veta'); zeros(nunstab, neta)]

    # -------------------- above are output for system in terms of z'y -------
    G1 = real(z*G1*z')
    C = real(z*C)
    impact = real(z * impact)
    loose = real(z * loose)

    # Correction 10/28/96:  formerly line below had real(z*ywt) on rhs, an error.
    ywt=z*ywt

    return G1, C, impact, fmat, fwt, ywt, gev, eu, loose
end

end # module

#=

Test case from Peifan's solution

%Log-linearized model should be written in Sims' form
%Γ0 * y(t) = Γ1 * y(t - 1) + C + Ψ * z(t) + Π * eta(t)
%y(t) = [Ex(t+1) xt Epi(t+1) pit it gt ut]
%z(t) = [epsilon_i epsilon_g epsilon_u]
%eta(t) = [err_x err_pi]

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
Γ0 = [-1  1  -σ  0  σ  -1  0
           0  -κ  -β  1  0  0  -1
           0  -(1-ρ_i)*ϕ_x  0  -(1-ρ_i)*ϕ_π  1  0  0
           0  0  0  0  0  1  0
           0  0  0  0  0  0  1
           0  1  0  0  0  0  0
           0  0  0  1  0  0  0]
Γ1 = [0  0  0  0  0  0  0
          0  0  0  0  0  0  0
          0  0  0  0  ρ_i  0  0
          0  0  0  0  0  ρ_g  0
          0  0  0  0  0  0  ρ_u
          1  0  0  0  0  0  0
          0  0  1  0  0  0  0]
C = [0  0  0  0  0  0  0]'

Ψ = [0  0  0
       0  0  0
       σ_i  0  0
       0  σ_g  0
       0  0  σ_u
       0  0  0
       0  0  0]
Π = [0  0
      0  0
      0  0
      0  0
      0  0
      1  0
      0  1]

%Thanks Sims!
[A, ~, B, ~, ~, ~, ~, eu, ~] = gensys(Γ0, Γ1, C, Ψ, Π) %#ok<ASGLU>
%eu = (1, 1) => existence and uniqueness




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