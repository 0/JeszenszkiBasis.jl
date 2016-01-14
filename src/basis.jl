abstract AbstractSzbasis


"""
Basis of occupation vectors for K sites and N particles.
"""
immutable Szbasis <: AbstractSzbasis
    # Number of sites.
    K::Int
    # Number of particles.
    N::Int

    # Number of basis vectors.
    D::Int
    # Occupation vectors (K by D).
    vectors::Array{Int, 2}
end

function Szbasis(K::Int, N::Int)
    # At least 1 site.
    K >= 1 || throw(DomainError())
    # At least 0 particles.
    N >= 0 || throw(DomainError())

    # Basis size.
    D = num_vectors(N, K)

    v = zeros(Int, K)
    v[1] = N
    vectors = Array(Int, K, D)
    vectors[:, 1] = v

    for i=2:D
        if v[1] > 0
            v[1] -= 1
            v[2] += 1
        else
            j = findfirst(v)

            v[1] = v[j] - 1
            v[j] = 0
            v[j+1] += 1
        end
        for j=1:K
            vectors[j, i] = v[j]
        end
    end

    Szbasis(K, N, D, vectors)
end


"""
Basis of occupation vectors for K sites and N particles, with no more than M
particles per site.
"""
immutable RestrictedSzbasis <: AbstractSzbasis
    # Number of sites.
    K::Int
    # Number of particles.
    N::Int
    # Site capacity.
    M::Int

    # Number of basis vectors.
    D::Int
    # Occupation vectors (K by D).
    vectors::Array{Int, 2}
end

function RestrictedSzbasis(K::Int, N::Int, M::Int)
    # At least 1 site.
    K >= 1 || throw(DomainError())
    # At least 0 particles.
    N >= 0 || throw(DomainError())
    # No more than can fit.
    N <= K * M || throw(DomainError())

    # Basis size.
    D = num_vectors(N, K, M)
    dNM = M > 0 ? div(N, M) : 1

    v = zeros(Int, K)
    for j=1:dNM
        v[j] = M
    end
    if 1 <= dNM + 1 <= K
        v[dNM+1] = N - M * dNM
    end
    vectors = Array(Int, K, D)
    vectors[:, 1] = v

    for i=2:D
        if v[1] > 0
            if v[1] < M
                delta = M - v[1]
                v[1] = M
            else
                delta = 0
            end

            j = findfirst(v .< M)

            v[j] += 1
            v[j-1] -= 1 + delta
        else
            j = findfirst(v)
            k = j + findfirst(v[j+1:end] .< M)

            v[k-j] = v[j] - 1
            v[k] += 1
            for l=1:k-j-1
                v[l] = M
            end
            # The indices after the first one differ from those in the paper.
            for l=k-j+1:k-1
                v[l] = 0
            end
        end
        for j=1:K
            vectors[j, i] = v[j]
        end
    end

    RestrictedSzbasis(K, N, M, D, vectors)
end
