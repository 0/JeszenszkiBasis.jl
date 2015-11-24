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
