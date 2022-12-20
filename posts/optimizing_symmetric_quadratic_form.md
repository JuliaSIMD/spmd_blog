+++
title = "Optimizing Symmetric Quadratic Form"
hascode = true
date = Date(2022, 12, 19)
rss = ""
+++

@def tags = ["julia", "optimization", "quadratic-form", "loop-transformation", "linear-algebra"]


This post is motivated by [this](https://discourse.julialang.org/t/faster-quadratic-expression-for-symmetric-matrices/91843) discourse thread, providing the following code:
```julia
using LinearAlgebra, MKL, BenchmarkTools, LoopVectorization
product_1(x, M) = dot(x, M, x)

function product_2(x, M) 
	Mx = M * x
	return dot(x, Mx)
end
function product_3(x, M)
    n = length(x)
    result = zero(eltype(x))

    @inbounds for i = 1:n
        result += x[i]^2 * M[i, i]

        @simd for j = 1 : i - 1
            result += 2* x[i] * x[j] * M[j, i]
        end 
    end

    return result
end
```
The author expressed surprise that making matrix `M` symmetric wasn't always faster, and also that `product_3`, explicitly taking advantage of that symmetry, was slowest.

Setting a baseline:
```julia
BLAS.set_num_threads(1) # single threaded `product_2`
x = rand(200);
A = rand(length(x), 2+length(x)) |> y -> y*y';
B = Symmetric(A, :U);
@btime product_1($x, $A)
@btime product_2($x, $A)
@btime product_3($x, $A)
@btime product_1($x, $B)
@btime product_2($x, $B)
@btime product_3($x, $B)
```
`product_1` and `product_3` are implemented in "pure Julia", while `product_2` calls out to the BLAS library for both the matrix-vector multiply and the vector-vector dot product. I get
```julia
julia> @btime product_1($x, $A)
  3.662 μs (0 allocations: 0 bytes)
508214.49429146934

julia> @btime product_2($x, $A)
  2.366 μs (1 allocation: 1.77 KiB)
508214.49429146945

julia> @btime product_3($x, $A)
  5.368 μs (0 allocations: 0 bytes)
508214.49429146847

julia> @btime product_1($x, $B)
  7.090 μs (0 allocations: 0 bytes)
508214.49429146847

julia> @btime product_2($x, $B)
  1.734 μs (1 allocation: 1.77 KiB)
508214.49429146945

julia> @btime product_3($x, $B)
  8.916 μs (0 allocations: 0 bytes)
508214.49429146847

julia> versioninfo()
Julia Version 1.10.0-DEV.175
Commit 0d469bd0f8* (2022-12-18 20:11 UTC)
Platform Info:
  OS: Linux (x86_64-redhat-linux)
  CPU: 36 × Intel(R) Core(TM) i9-10980XE CPU @ 3.00GHz
  WORD_SIZE: 64
  LIBM: libopenlibm
  LLVM: libLLVM-14.0.6 (ORCJIT, cascadelake)
  Threads: 36 on 36 virtual cores
```
My CPU is cascadelake, an Intel Skylake-X clone feature AVX512 and 1MiB L2 cache/core.
MKL is winning, especially for symmetric matrices, despite the fact we've limited it to only a single thread, and that it needs to allocate and write memory. `using MKL` is encouraged for best performance when possible, the default library OpenBLAS is likely to be slower.

Coming up second best will not do, especially when there's such an obvious deficit in `product_2` -- the need to allocate and write memory.

Asymptotically, the performance is limited by the rate at which we can loading memory from the matrix. Thus we should certainly see `Symmetric` winning, at least so long as we take advantage of the structure to avoid loading any elements more than once.

While LoopModels will be able to handle symmetry and arbitrary affine loop nests, LoopVectorization.jl is unfortunately limited to only rectangular loop nests, thus it cannot handle the symmetric case. I also plan on adding cache tiling to LoopModels, but this is also something LoopVectorization.jl does not support, so performance likely will not scale to larger sizes, nor do we consider things like alignment. But a 200x200 matrix fits neatly in cache, and has column sizes that are an integer multiple of 64 bytes, so none of this matters here.

```julia
function dotturbo(x, A)
    s = zero(promote_type(eltype(x), eltype(A)))
    @turbo for n ∈ axes(A, 2)
        t = zero(s)
        for m ∈ axes(A, 1)
            t += x[m] * A[m, n]
        end
        s += t * x[n]
    end
    s
end
```
Note that we manually implemented an optimization also used in `LinearAlgebra.dot`: we hoisted the multiplication by `x[n]` out of the loop.
It's also my intention to automate this, but `LoopVectorization.jl`'s transforms are fairly limited, and my performance development time has of course moved elsewhere.

I get
```julia
julia> @btime dotturbo($x, $A)
  1.875 μs (0 allocations: 0 bytes)
508214.49429146945
```
At least we're now winning the dense case. Thanks to the occasional GC activity from memory allocations, we're not too far behind the `Symmetric` case in average performance:
```julia
julia> @benchmark dotturbo($x, $A)
BenchmarkTools.Trial: 10000 samples with 10 evaluations.
 Range (min … max):  1.869 μs …  6.652 μs  ┊ GC (min … max): 0.00% … 0.00%
 Time  (median):     1.887 μs              ┊ GC (median):    0.00%
 Time  (mean ± σ):   1.892 μs ± 56.790 ns  ┊ GC (mean ± σ):  0.00% ± 0.00%

      ▄█▇▅▃
  ▂▂▃▆██████▅▃▃▂▂▂▂▂▂▂▂▂▂▂▂▁▂▁▂▂▁▂▂▁▁▁▁▂▁▁▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂ ▃
  1.87 μs        Histogram: frequency by time        2.03 μs <

 Memory estimate: 0 bytes, allocs estimate: 0.

julia> @benchmark product_2($x, $B)
BenchmarkTools.Trial: 10000 samples with 10 evaluations.
 Range (min … max):  1.733 μs … 107.171 μs  ┊ GC (min … max): 0.00% … 93.21%
 Time  (median):     1.768 μs               ┊ GC (median):    0.00%
 Time  (mean ± σ):   1.843 μs ±   1.972 μs  ┊ GC (mean ± σ):  2.04% ±  1.88%

    ▁█▇▁
  ▁▃████▅▃▃▂▂▂▂▂▂▃▃▃▃▃▃▂▂▂▂▂▂▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁▁ ▂
  1.73 μs         Histogram: frequency by time        2.12 μs <

 Memory estimate: 1.77 KiB, allocs estimate: 1.
 ```
 1.892 vs 1.843 microseconds.

Smelling blood in the water, lets take advantage of symmetry to win this benchmark.

Lets start with `product_3`, and apply one transform at a time to discuss it's problems.
First, we can hoist a lot of arithmetic.
```julia
function product_4(x, _M)
    n = length(x)
    result_diag = zero(eltype(x))
    result_supd = zero(eltype(x))
	M = _M isa Symmetric ? parent(_M) : _M
	Base.require_one_based_indexing(x)
	Base.require_one_based_indexing(M)
    @inbounds for i = 1:n
        result_diag += x[i]^2 * M[i, i]
	    result_innr = zero(eltype(x))
        @simd for j = 1 : i - 1
            result_innr += x[j] * M[j, i]
        end
        result_supd += x[i]*result_innr
    end
    return 2result_supd + result_diag
end
```
yielding
```julia
julia> @btime product_4($x, $B)
  3.977 μs (0 allocations: 0 bytes)
508214.49429146934
```
Better, now we're more in line with `product_1`.
If our math is slow, maybe we just need to ask the compiler to make it fast?
```julia
@fastmath function product_5(x, _M) # only change is `@fastmath`
    n = length(x)
    result_diag = zero(eltype(x))
    result_supd = zero(eltype(x))
	M = _M isa Symmetric ? parent(_M) : _M
	Base.require_one_based_indexing(x)
	Base.require_one_based_indexing(M)
    @inbounds for i = 1:n
        result_diag += x[i]^2 * M[i, i]
	    result_innr = zero(eltype(x))
        @simd for j = 1 : i - 1
            result_innr += x[j] * M[j, i]
        end
        result_supd += x[i]*result_innr
    end
    return 2result_supd + result_diag
end
```
so we're now beating `product_1` with dense arrays:
```julia
julia> @btime product_5($x, $B)
  3.503 μs (0 allocations: 0 bytes)
508214.49429146934
```
But we're still a far cry from `@turbo`. Unfortunately, because `LoopVectorization.jl` only supports rectangular loops, we're limited to applying it to the inner-most loop, hamstringing it:
```julia
@fastmath function product_6(x, _M) # only change is `@fastmath`
    n = length(x)
	n == 0 && return zero(eltype(x))
	M = _M isa Symmetric ? parent(_M) : _M
	Base.require_one_based_indexing(x)
	Base.require_one_based_indexing(M)
    result_diag = x[1]^2 * M[1, 1]
    result_supd = zero(eltype(x))
	# `@turbo` doesn't support length-0 iteration
	# so we need to make sure length(1:i-1) > 0
    @inbounds for i = 2:n 
        result_diag += x[i]^2 * M[i, i]
	    result_innr = zero(eltype(x))
        @turbo for j = 1 : i - 1
            result_innr += x[j] * M[j, i]
        end
        result_supd += x[i]*result_innr
    end
    return 2result_supd + result_diag
end
```
yet
```julia
julia> @btime product_6($x, $B)
  1.795 μs (0 allocations: 0 bytes)
508214.4942914694

julia> @benchmark product_6($x, $B)
BenchmarkTools.Trial: 10000 samples with 10 evaluations.
 Range (min … max):  1.663 μs …  5.431 μs  ┊ GC (min … max): 0.00% … 0.00%
 Time  (median):     1.680 μs              ┊ GC (median):    0.00%
 Time  (mean ± σ):   1.717 μs ± 73.691 ns  ┊ GC (mean ± σ):  0.00% ± 0.00%

   ▂▇█▅▁                               ▁
  ▃█████▅▅▃▂▂▂▂▂▂▁▁▁▂▁▁▂▁▁▂▁▂▂▂▁▁▁▁▂▂▂▄█▇▆▅▃▃▂▂▂▂▂▂▂▂▂▂▂▂▂▂▂ ▃
  1.66 μs        Histogram: frequency by time        1.88 μs <

 Memory estimate: 0 bytes, allocs estimate: 0.
```
I'm not sure why the minimum time from `@btime` was so much higher than the mean time from `@benchmark`.

But we're now at least in clear contention with MKL's `symv!`+`dot`.
So, if we could only apply `@turbo` to the inner loop, why did it help so much? LLVM is really bad at producing code for architectures with big vectors, and AVX512 in particular. On a Zen3 CPU with AVX2, I got
```julia
julia> @btime product_5($x, $B)
  2.537 μs (0 allocations: 0 bytes)
502848.6107780409

julia> @btime product_6($x, $B)
  2.063 μs (0 allocations: 0 bytes)
502848.61077804095
```
while on an Apple M1, I got
```julia
julia> @btime product_5($x, $B)
  3.604 μs (0 allocations: 0 bytes)
580770.6876433745

julia> @btime product_6($x, $B)
  3.521 μs (0 allocations: 0 bytes)
580770.6876433745
```
so the performance advantage decreases as the vector width decreases. This is mostly because LLVM doesn't handle remainders efficiently. The larger the vectors, the larger the remainders on average, the worse LLVM's performance.

So, what's left, what's missing?
What would `@turbo` on an outerloop do differently?
For one thing, reducing vectors to scalars is expensive, so `@turbo` will only do so once at the end, while our code is doing so on every iteration of the outer loop.
How can we avoid doing so?

We can tell `@turbo` we want to accumulate into a SIMD vector instead of a scalar by defining our accumulator to be a `VectorizationBase.Vec` instead of a scalar.
Then, we just have to accumulate the vector into a scalar ourselves before returning the value:
```julia
using VectorizationBase
@fastmath function product_7(x, _M) # only change is `@fastmath`
    n = length(x)
	T = eltype(x)
	n == 0 && return zero(T)
	M = _M isa Symmetric ? parent(_M) : _M
	Base.require_one_based_indexing(x)
	Base.require_one_based_indexing(M)
    result_diag = Vec(zero(T))
	@turbo for i = 1:n
		result_diag += x[i]^2 * M[i, i]
	end
    result_supd = Vec(zero(T))
	# `@turbo` doesn't support length-0 iteration
	# so we need to make sure length(1:i-1) > 0
    @inbounds for i = 2:n 
	    result_innr = Vec(zero(T))
        @turbo for j = 1 : i - 1
            result_innr += x[j] * M[j, i]
        end
        result_supd += x[i]*result_innr
    end
    return VectorizationBase.vsum(2result_supd + result_diag)
end
```
Now I get:
```julia
julia> @btime product_7($x, $B)
  1.354 μs (0 allocations: 0 bytes)
508214.4942914695
```
Well ahead of any other implementation.

The other major optimization `@turbo` does when applied to the outer loop is to cut down on the `x[j]` reloads in the inner `j` loop by unrolling the outer `i` loop. That is, if we unroll `M[j,i]` by 8, we can reuse an `x[j]` load 8 times before discarding it and moving on to the next inner loop iteration. In this way, we need $\frac{1}{8}$ as many loads from $x$ in the inner loop.

```julia
function product_8(x, _M)
    N = length(x)
	T = eltype(x)
	N == 0 && return zero(T)
	M = _M isa Symmetric ? parent(_M) : _M
	Base.require_one_based_indexing(x)
	Base.require_one_based_indexing(M)
	# we unroll the outer loop by 8
	remainder = N & 7
	n = remainder
	if n != N # guard
	    # 8 accumulators (for the outer loop)
		Base.Cartesian.@nexprs 4 i -> acc_i = Vec(zero(T))
		while true
			Base.Cartesian.@nexprs 8 i -> a_i = Vec(zero(T))
			# iterations i = n+1:n+8
			if n != 0
			    # we'll loop normally up through n
			    @turbo for j = 1:n
				    Base.Cartesian.@nexprs 8 i -> begin
						a_i += x[j] * M[j,i+n]
					end
				end
			end
			# 8x8 diagonal block
			@turbo for j = 1:8
			    Base.Cartesian.@nexprs 8 i -> begin
			        mask_i = (j<i) + 0.5*(j==i)
			        a_i += x[j+n] * (M[j+n,i+n] * mask_i)
				end
			end
			Base.Cartesian.@nexprs 4 i -> acc_i += x[i+n]*a_i
			Base.Cartesian.@nexprs 4 i -> acc_i += x[4+i+n]*a_{4+i}
			n += 8
			n == N && break
		end # while true
		acc_1 += acc_3
		acc_2 += acc_4
		acc_1 += acc_2
		ret = 2*VectorizationBase.vsum(acc_1)
	    remainder == 0 && return ret
	else
	    ret = zero(T)
	end
	# we know remainder != 0, because N == 0 returned early
	r = Base.oneto(remainder)
    # product_5 seemed fast for small inputs
	return ret + @views product_5(x[r], M[r,r])
end
```
Now performance is starting to look pretty good:
```julia
julia> @btime product_8($x, $B)
  987.308 ns (0 allocations: 0 bytes)
508214.49429146945
```
We're now under 1 microsecond, leaving all the competition well behind, so I think I'll call it quits for today.

Note however that this implementation optimistically assumes that `stride(A,2)%8==0`.

If we want to try and address alignment, we could consider taking `product_7`, and then initially pealing off enough initial iterations of the inner loop to align it. Not crossing 64 byte boundaries with loads from `A`  will make them cheaper, but at the cost of reloading `x` more often.
Given aligned loads and a vector width of `8`, we could estimate the cost of this as roughly $\left(\frac{N}{8}\right)\left(\frac{N}{8}\right) + \frac{N(N+1)}{2\cdot 8} + N = \frac{N^2}{64} + \frac{(N(N+1))}{16} + N$.
The $\left(\frac{N}{8}\right)\left(\frac{N}{8}\right)$ reflects the $\frac{N}{8}$ outer loop iterations times the $\frac{N}{W}=\frac{N}{8}$ loads per inner loop iteration.
The $\frac{N(N+1)}{2}$ are the number of elements we're loading from `A`, then the $\frac{1}{8}$ is again dividing by $W$. 
The final $+ N$ are the loads `x[i]`, insignificant thanks to the smaller exponent.

The base pointer of a large array, e.g. `A`, will normally be 64 byte aligned.
```julia
julia> Int(pointer(A))%64
0
```
Assuming this, the worst case scenario for alignment would be that `isodd(stride(A,2)%8)`. Then, with AVX512, only 1 out of 8 columns will be aligned.

With that, the cost of loading from `A` becomes 
\begin{align}
\left(\frac{N(N+1)}{16}\right)
\left(\frac{2\cdot 7}{8}\right)
+
\left(\frac{N(N+1)}{16}\right)
\left(\frac{1}{8}\right)
=
\frac{15N(N+1)}{128}
\end{align}
That is, loading from $A$ becomes $\frac{15}{8}$ times more expensive. We have an increase in cost of $\frac{7N(N+1)}{128}$.

If we try to address this by pealing loop iterations at the start of each column to align the loads from $A$, this causes one of two complications. Either
1. We can no longer unroll the outer loop, increasing the cost of reloading $x$ in the inner loop from $\frac{N^2}{64}$ to $\frac{N^2}{8}$, an increase in cost of $\frac{16N^2}{128} - \frac{2N^2}{128} = \frac{14}{128}$, meaning this is potentially a smaller increase in cost than what the unrolled `product_8` suffers.
2. Still unroll, and use shuffle/permute instructions to try and line up the loads of `x` correctly with the loads from `A` that may now all be offset to align them. This approach would require generating different specializations for the different offset patterns, where you'd have fixed shuffle sequences in each one. As we're dominated by loads, the shuffles may be "free"/able to execute concurrently, making this likely the fastest, but most complicated, approach. This is especially the case for Zen4, which has fairly fast shuffles but slow loads, or for consumer AVX512 Intel CPUs, which can use port 5 for shuffle instructions but not for floating point arithmetic, again making shuffles more or less "free" to execute concurrently with loads and floating point.
