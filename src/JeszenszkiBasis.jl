module JeszenszkiBasis

export
    Szbasis,
    num_vectors,
    serial_num,
    sub_serial_num

num_vectors(N, K) = binomial(N+K-1, K-1)

"""
Basis of occupation vectors for K sites and N particles.
"""
immutable Szbasis
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
    # At least 1 particle.
    N >= 1 || throw(DomainError())

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
Deterimine the serial number of an occupation vector.
"""
function serial_num(K::Int, N::Int, v::AbstractArray{Int, 1})
    I = 1

    for mu=1:K
        s = 0
        for nu=mu+1:K
            s += v[nu]
        end
        for i=0:v[mu]-1
            I += num_vectors(N-s-i, mu-1)
        end
    end

    I
end

serial_num(basis::Szbasis, v::AbstractArray{Int, 1}) = serial_num(basis.K, basis.N, v)

"""
Determine the serial number of a reduced occupation vector (containing a subset
of the sites).
"""
function sub_serial_num(::Szbasis, v::AbstractArray{Int, 1})
    K = length(v)
    N = sum(v)

    # Only one way to have no sites.
    K >= 1 || return 1
    # Only one way to have no particles.
    N >= 1 || return 1

    # Count the zero-particle case.
    I = 1

    for n=1:N-1
        I += num_vectors(n, K)
    end

    I + serial_num(K, N, v)
end

end
