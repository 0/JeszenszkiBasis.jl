using JeszenszkiBasis
using Base.Test

### Szbasis, num_vectors, site_max, serial_num

for (K, N, D) in [(50, 0, 1), (1, 1, 1), (2, 3, 4), (3, 2, 6), (4, 4, 35)]
    sb = Szbasis(K, N)

    # Numbers are correct.
    @test sb.K == K
    @test sb.N == N
    @test sb.D == D
    @test num_vectors(N, K) == D
    @test size(sb.vectors) == (sb.K, sb.D)
    @test length(sb) == D
    @test site_max(sb) == N

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
    if N > 0
        @test !(v in sb)
    else
        @test v in sb
    end
    for n in [-1, N+1]
        v[1] = n
        @test !(v in sb)
    end
end

@test_throws DomainError Szbasis(0, 5)
@test_throws DomainError Szbasis(5, -1)


### RestrictedSzbasis

for (K, N, M, D) in [(50, 0, 0, 1), (1, 1, 1, 1), (2, 3, 2, 2), (3, 2, 1, 3), (4, 4, 2, 19)]
    sb = RestrictedSzbasis(K, N, M)

    # Numbers are correct.
    @test sb.K == K
    @test sb.N == N
    @test sb.M == M
    @test sb.D == D
    @test num_vectors(N, K, M) == D
    @test size(sb.vectors) == (sb.K, sb.D)
    @test length(sb) == D
    @test site_max(sb) == M

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
    if N > 0
        @test !(v in sb)
    else
        @test v in sb
    end
    for n in [-1, N+1, M+1]
        v[1] = n
        @test !(v in sb)
    end
end

@test_throws DomainError RestrictedSzbasis(2, 5, 2)


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

let sb = RestrictedSzbasis(5, 2, 2)
    # Split each vector into 2 sites on the left and 3 on the right.
    left = []
    right = []

    for v in sb
        # Try range indexing.
        push!(left, sub_serial_num(sb, v[1:2]))
        # Try a SubArray.
        push!(right, sub_serial_num(sb, sub(v, 3:5)))
    end

    @test sort(left) == [1, 1, 1, 1, 1, 1, 2, 2, 2, 3, 3, 3, 4, 5, 6]
    @test sort(right) == [1, 1, 1, 2, 2, 3, 3, 4, 4, 5, 6, 7, 8, 9, 10]
end
