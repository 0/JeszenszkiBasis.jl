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
