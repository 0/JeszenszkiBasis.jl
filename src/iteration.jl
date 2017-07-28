## Occupation vector iteration and array properties.

mutable struct SzbasisIterState
    i::Int
end

function Base.start(basis::AbstractSzbasis)
    SzbasisIterState(0)
end

function Base.next(basis::AbstractSzbasis, state::SzbasisIterState)
    state.i += 1

    @view(basis.vectors[:, state.i]), state
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
