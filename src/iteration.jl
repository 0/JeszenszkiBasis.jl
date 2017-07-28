## Occupation vector iteration and array properties.

type SzbasisIterState
    i::Int
    v::Vector{Int}
end

function Base.start(basis::AbstractSzbasis)
    SzbasisIterState(0, Vector{Int}(basis.K))
end

function Base.next(basis::AbstractSzbasis, state::SzbasisIterState)
    state.i += 1

    for j in 1:basis.K
        state.v[j] = basis.vectors[j, state.i]
    end

    state.v, state
end

function Base.done(basis::AbstractSzbasis, state::SzbasisIterState)
    state.i == basis.D
end

Base.eltype(::Type{AbstractSzbasis}) = Vector{Int}
Base.length(basis::AbstractSzbasis) = basis.D

function Base.in(v::AbstractVector{Int}, basis::Szbasis)
    length(v) == basis.K && sum(v) == basis.N
end

function Base.in(v::AbstractVector{Int}, basis::RestrictedSzbasis)
    length(v) == basis.K && sum(v) == basis.N && maximum(v) <= basis.M
end
