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
    @test length(sb) == D

    # Occupations are correct.
    @test all(sum(sb.vectors, 1) .== N)

    # Each vector is unique.
    @test sb.vectors == unique(sb.vectors, 2)

    # Iteration works, serial numbers match up.
    for (i, v) in enumerate(sb)
        @test v in sb
        @test v == sb.vectors[:, i]
        @test i == serial_num(sb, v)
    end

    # Invalid vectors.
    v = zeros(Int, K)
    @test !(v in sb)
    v[1] = N+1
    @test !(v in sb)
end

@test_throws DomainError Szbasis(0, 5)
@test_throws DomainError Szbasis(5, 0)


### sub_serial_num

let sb = Szbasis(3, 2)
    # Split each vector into 2 sites on the left and 1 on the right.
    left = []
    right = []

    for v in sb
        # Try range indexing.
        push!(left, sub_serial_num(sb, v[1:2]))
        # Try a SubArray.
        push!(right, sub_serial_num(sb, sub(v, 3:3)))
    end

    @test sort(left) == [1, 2, 3, 4, 5, 6]
    @test sort(right) == [1, 1, 1, 2, 2, 3]
end
