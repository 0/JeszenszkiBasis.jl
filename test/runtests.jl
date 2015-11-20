using JeszenszkiBasis
using Base.Test

### Szbasis, num_vectors, serial_num

for (K, N, D) in [(1, 1, 1), (2, 3, 4), (3, 2, 6), (4, 4, 35)]
    sb = Szbasis(K, N)

    # Numbers are correct.
    @test sb.K == K
    @test sb.N == N
    @test sb.D == D
    @test num_vectors(N, K) == D
    @test size(sb.vectors) == (sb.K, sb.D)

    # Occupations are correct.
    @test all(sum(sb.vectors, 1) .== N)

    # Each vector is unique.
    @test sb.vectors == unique(sb.vectors, 2)

    # Serial numbers match up.
    for i=1:sb.D
        @test i == serial_num(sb, sb.vectors[:, i])
    end
end

@test_throws DomainError Szbasis(0, 5)
@test_throws DomainError Szbasis(5, 0)


### sub_serial_num

let sb = Szbasis(3, 2)
    # Split each vector into 2 sites on the left and 1 on the right.
    left = []
    right = []

    for i=1:sb.D
        # Try range indexing.
        push!(left, sub_serial_num(sb, sb.vectors[1:2, i]))
        # Try a SubArray.
        push!(right, sub_serial_num(sb, sub(sb.vectors, 3:3, i)))
    end

    @test sort(left) == [1, 2, 3, 4, 5, 6]
    @test sort(right) == [1, 1, 1, 2, 2, 3]
end
