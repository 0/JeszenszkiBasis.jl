"""
Number of vectors in a basis.
"""
num_vectors(N, K) = binomial(N+K-1, K-1)

# Global cache of generalized Pascal's triangles for restricted num_vectors.
const triangles = Dict{Int, Array{Int, 2}}()

function num_vectors(N, K, M)
    0 <= N <= M * K || return 0
    N == K == 0 && return 1

    # Create a new triangle.
    if !haskey(triangles, M)
        triangles[M] = zeros(Int, 1, 1)
        triangles[M][1, 1] = 1
    end

    # Extend an existing triangle.
    if size(triangles[M], 1) < K + 1
        t_old = triangles[M]
        K_old = size(t_old, 1) - 1
        t = zeros(Int, K+1, M*K+1)
        for k=0:K_old
            for n=0:M*k
                t[k+1, n+1] = t_old[k+1, n+1]
            end
        end
        for k=K_old+1:K
            for n=0:M*k
                for m=0:min(M, n)
                    t[k+1, n+1] += t[k, n+1-m]
                end
            end
        end
        triangles[M] = t
    end

    triangles[M][K+1, N+1]
end

num_vectors(::Szbasis, N, K) = num_vectors(N, K)
num_vectors(basis::RestrictedSzbasis, N, K) = num_vectors(N, K, basis.M)


"""
Maximum number of particles in a site.
"""
site_max(basis::Szbasis) = basis.N
site_max(basis::RestrictedSzbasis) = basis.M


"""
Serial number of an occupation vector.
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
serial_num(::Szbasis, K::Int, N::Int, v::AbstractArray{Int, 1}) = serial_num(K, N, v)

function serial_num(K::Int, N::Int, M::Int, v::AbstractArray{Int, 1})
    I = 1

    for mu=1:K
        s = 0
        for nu=mu+1:K
            s += v[nu]
        end
        for i=0:v[mu]-1
            I += num_vectors(N-s-i, mu-1, M)
        end
    end

    I
end

serial_num(basis::RestrictedSzbasis, v::AbstractArray{Int, 1}) = serial_num(basis.K, basis.N, basis.M, v)
serial_num(basis::RestrictedSzbasis, K::Int, N::Int, v::AbstractArray{Int, 1}) = serial_num(K, N, basis.M, v)


"""
Serial number of a reduced occupation vector (with a subset of the sites).
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

function sub_serial_num(basis::RestrictedSzbasis, v::AbstractArray{Int, 1})
    K = length(v)
    N = sum(v)

    # Only one way to have no sites.
    K >= 1 || return 1
    # Only one way to have no particles.
    N >= 1 || return 1

    # Count the zero-particle case.
    I = 1

    for n=1:N-1
        I += num_vectors(n, K, basis.M)
    end

    I + serial_num(K, N, basis.M, v)
end
