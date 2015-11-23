module JeszenszkiBasis

export
    AbstractSzbasis,
    Szbasis,
    RestrictedSzbasis,

    num_vectors,
    serial_num,
    site_max,
    sub_serial_num


num_vectors(N, K) = binomial(N+K-1, K-1)

function num_vectors(N, K, M)
    0 <= N <= M * K || return 0
    N == K == 0 && return 1

    result = 0
    for L=N-M:N
        result += num_vectors(L, K-1, M)
    end
    result
end


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
    # At least 1 particle.
    N >= 1 || throw(DomainError())
    # No more than can fit.
    N <= K * M || throw(DomainError())

    # Basis size.
    D = num_vectors(N, K, M)

    v = zeros(Int, K)
    for j=1:div(N, M)
        v[j] = M
    end
    if 1 <= div(N, M) + 1 <= K
        v[div(N, M)+1] = N - M * div(N, M)
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


type SzbasisIterState
    i::Int
    v::Array{Int, 1}
end

function Base.start(basis::AbstractSzbasis)
    SzbasisIterState(0, Array(Int, basis.K))
end

function Base.next(basis::AbstractSzbasis, state::SzbasisIterState)
    state.i += 1

    for j=1:basis.K
        state.v[j] = basis.vectors[j, state.i]
    end

    state.v, state
end

function Base.done(basis::AbstractSzbasis, state::SzbasisIterState)
    state.i == basis.D
end

Base.eltype(::Type{AbstractSzbasis}) = Array{Int, 1}
Base.length(basis::AbstractSzbasis) = basis.D

function Base.in(v::AbstractArray{Int, 1}, basis::Szbasis)
    length(v) == basis.K && sum(v) == basis.N
end

function Base.in(v::AbstractArray{Int, 1}, basis::RestrictedSzbasis)
    length(v) == basis.K && sum(v) == basis.N && maximum(v) <= basis.M
end


"""
Maximum number of particles in a site.
"""
site_max(basis::Szbasis) = basis.N
site_max(basis::RestrictedSzbasis) = basis.M


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

end
