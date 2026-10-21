"""
    get_estimation_of_tau(A; target_ratio=10.0, max_evaluations=40)

Estimate an ILU drop tolerance on a leading principal submatrix. The discrete
fill ratio need not attain the requested target. Search is bounded and returns
the closest sampled ratio not exceeding the target (or the sparsest sample).
"""
function get_estimation_of_tau(A; target_ratio=10.0, max_evaluations::Integer=40)
    size(A, 1) == size(A, 2) && size(A, 1) > 0 ||
        throw(ArgumentError("expected a nonempty square matrix"))
    isfinite(target_ratio) && target_ratio > 0 ||
        throw(ArgumentError("target_ratio must be finite and positive"))
    max_evaluations >= 2 || throw(ArgumentError("max_evaluations must be at least 2"))
    n_small = max(1, div(size(A, 1), 5))
    A_small = A[1:n_small, 1:n_small]
    count_A = nnz(A_small)
    count_A > 0 || throw(ArgumentError("cannot estimate ILU fill for a zero matrix"))
    samples = Tuple{Float64, Float64}[]
    function evaluate(tau)
        ratio = nnz(ilu(A_small, τ=tau)) / count_A
        isfinite(ratio) || error("nonfinite ILU fill ratio")
        push!(samples, (Float64(tau), ratio))
        return ratio
    end
    # Exact no-drop factorization is safe only for tiny samples. On large
    # domains it can require far more memory than the requested fill target.
    tiny_sample = n_small <= 64
    if tiny_sample
        ratio_zero = evaluate(0.0)
        ratio_zero <= target_ratio && return 0.0
    end
    low = 0.0
    high = max(1.0, maximum(abs, nonzeros(A_small)))
    ratio_high = evaluate(high)
    while ratio_high > target_ratio && length(samples) < max_evaluations
        low = high
        high >= floatmax(Float64) / 10 && break
        high *= 10
        ratio_high = evaluate(high)
    end
    bracketed = tiny_sample || low > 0
    if !bracketed && ratio_high <= target_ratio
        low = high / 10
        while low > 0 && length(samples) < max_evaluations
            ratio_low = evaluate(low)
            if ratio_low > target_ratio
                bracketed = true
                break
            end
            high = low
            ratio_high = ratio_low
            low /= 10
        end
    end
    # Bisection needs no differentiability or division by a change in fill.
    while bracketed && ratio_high <= target_ratio && length(samples) < max_evaluations
        mid = low + (high - low) / 2
        (mid == low || mid == high) && break
        ratio = evaluate(mid)
        abs(ratio - target_ratio) <= 0.5 && return mid
        if ratio > target_ratio
            low = mid
        else
            high = mid
            ratio_high = ratio
        end
    end
    feasible = filter(sample -> sample[2] <= target_ratio, samples)
    candidates = isempty(feasible) ? samples : feasible
    best = candidates[argmin([abs(sample[2] - target_ratio) for sample in candidates])]
    @warn "ILU fill search reached its limit; using best sampled tolerance" tau=best[1] fill_ratio=best[2] target_ratio
    return best[1]
end
