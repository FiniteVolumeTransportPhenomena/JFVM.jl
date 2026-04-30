using JFVM

# Simple median benchmark helper without extra dependencies.
function median_metrics(f; samples::Int=5, warmup::Int=1)
    for _ in 1:warmup
        f()
    end

    times = Vector{Float64}(undef, samples)
    bytes = Vector{Int}(undef, samples)

    for i in 1:samples
        GC.gc()
        m = @timed f()
        times[i] = m.time
        bytes[i] = m.bytes
    end

    sort!(times)
    sort!(bytes)
    mid = cld(samples, 2)
    return times[mid], bytes[mid]
end

function assert_leq(name::AbstractString, value::Real, limit::Real)
    println(name, ": ", value, " (limit ", limit, ")")
    if value > limit
        error("Performance gate failed for $(name): $(value) > $(limit)")
    end
end

# Representative 2D workload.
m = createMesh2D(64, 64, 1.0, 1.0)
phi = createCellVariable(m, 1.0)
phi_face = createFaceVariable(m, 0.0)
D = harmonicMean(phi)
u = createFaceVariable(m, 0.1)

linear_t, linear_b = median_metrics(() -> linearMean!(phi, phi_face))
diff_t, diff_b = median_metrics(() -> diffusionTerm(D))
conv_t, conv_b = median_metrics(() -> convectionTerm(u))
trans_t, trans_b = median_metrics(() -> transientTerm(phi, 0.1))

println("--- JFVM performance gate report ---")
println("linearMean! time (s): ", linear_t)
assert_leq("linearMean! alloc (bytes)", linear_b, 800_000)

println("diffusionTerm time (s): ", diff_t)
assert_leq("diffusionTerm alloc (bytes)", diff_b, 10_000_000)

println("convectionTerm time (s): ", conv_t)
assert_leq("convectionTerm alloc (bytes)", conv_b, 10_000_000)

println("transientTerm time (s): ", trans_t)
assert_leq("transientTerm alloc (bytes)", trans_b, 2_000_000)

println("Performance gates passed.")
