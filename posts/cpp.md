+++
title = "C++"
hascode = true
date = Date(2023, 4, 06)
rss = "c++ intro"
+++

@def tags = ["c++", "intro"]

### Intro

Note, I am only talking about C++20 here.

C++ can be an intimidating language, with a lot of scary negative sentiment expressed online, often (but not always) from people with little to no experience with the language.

Here, I'll give a brief introduction to the language from the perspective of a Julia programmer, describing differences from the language, and giving suggestions for a simple style that should make for relatively simple and maintainable code.

The first thing to recognize is that things like function boundaries, assignments, and leaving a scope are significant in C++, unlike Julia. In Julia, `a = b` simply binds `a` to the same thing `b` is bound to. In C++, this instead calls `operator=`.
Similarly, there are many ways we can pass an object to a function in C++, but it's probably best to restrict ourselves to three for most circumstances:
1. pass by value. This copies an object, e.g. passing a vector by value copies the entire vector and all of its contents. The type signature of the function you call may look like `foo(int);` or `bar(std::vector<double>);`, i.e. no decoration. If you no longer need the object you pass in, you could `std::move` it to a function expecting the value. If you do this, instead of copying, you give away ownership.
2. Const reference, e.g. `foo(const int&);` or `bar(const std::vector<double>&);`. This lets functions `foo` and `bar` read and use the variables, but not modify them. For things cheap/trivial to copy, like `int`, pass by value should be preferred, but you should generally pass things like `std::vector`s by const reference.
3. By pointer, e.g. `foo(int*)` or `bar(std::vector<double>*)`. Do this when you want to mutate/modify the variable. Note that unary `&` almost always returns the pointer to a variable. Only use this when you're actually mutating the variable!

There are all sorts of combinations, but sticking to these three covers most use cases, is simple to remember, and should help make code readable. For example, if you see the call `buz(a, b, &c, d)`, then you know that `c` is being mutated, while `a`, `b`, and `d` are not. Think of this convention as similar to the convention in Julia of ending mutating function names with a `!`, except here we see specifically which variables are being mutated.

The next thing to keep in mind is that a value's destructor is called as soon as it leaves a scope.
Most objects in C++ have RAII; this means that creating the object allocates its resources (so creating or pushing into a vector allocates the memory it needs), and its destruction frees that memory.

```c++
#include <vector>

bar(const std::vector<double>&);

{
  std::vector<double> x{1.0, 5.0, 8.0}; // initializes x
  x.push_back(2); // x now contains {1.0, 5.0, 8.0, 2.0}
  bar(x); // we do something with `x`, the ref leaves scope, but that doesn't trigger anything
  bar(x); // same again
} // scope ends, `x` is destroyed and the memory is freed.
```

Returning such a value from a function or adding it to another collection will copy it (note, you can move it to avoid a copy).
Thus, safe memory management is rather effortless; you don't need to use `new`/`malloc` or `delete`/`free`, nor do you need a garbage collector.

### Generic programming

Julia features/makes heavy use of multiple dispatch.
Multiple dispatch is semantically a dynamic/runtime dispatch; C++ only has built-in support for single dispatch, because this is easy to optimize with vtables.
However, performant Julia code is heavily reliant on an optimization called devirtualization. That is, when the types are known at compile time, the runtime multiple dispatch can be resolved/turned into a fast static dispatch, or even inlined.

While C++ only has single dispatch, it has function overloading. Function overloading lets you pick the method that gets called based on the combination of all types, but this must happen at compile time.
That is, as long as all types are known at compile time -- which is the vast majority of use cases in Julia -- we can treat overloading in C++ just like multiple dispatch in Julia:
```c++
#include <concepts>
#include <cstdio>

auto foo(const auto &x, const auto &y) { printf("generic fallback method\n"); }
auto foo(const std::integral auto &x, const auto &y) {
  printf("first arg integral\n");
}
auto foo(const auto &x, const std::integral auto &y) {
  printf("second arg integral\n");
}
auto foo(const std::integral auto &x, const std::integral auto &y) {
  printf("both args integral\n");
}

int main() {
  foo(1.0, 2.0); // generic fallback method
  foo(1, 2.0); // first arg integral
  foo(1.0, 2); // second arg integral
  foo(1, 2); // both args integral
  return 0;
}
```
That is, you can think of overloading as multiple dispatch, with the requirement that your code must be type stable.

### Concepts

C++20 also brings introduces concepts, which allow you to control dispatches based on what types are capable of.
```c++
#include <array>
#include <cstddef>
#include <cstdio>
#include <vector>

template <typename T>
concept AbstractVector = requires(T t, size_t i) {
  { t.size() } -> std::convertible_to<size_t>;
  { t[i] } -> std::convertible_to<typename T::value_type>;
};

static_assert(AbstractVector<std::vector<double>>);
static_assert(AbstractVector<std::vector<int>>);
static_assert(AbstractVector<std::vector<bool>>);
static_assert(AbstractVector<std::array<ptrdiff_t, 4>>);
static_assert(!AbstractVector<double>);

// if an AbstractVector, get first index
auto first(const AbstractVector auto &v) { return v[0]; }
auto last(const AbstractVector auto &v) { return v[v.size() - 1]; }

auto first(const auto &x) { return x; }
auto last(const auto &x) { return x; }

int main() {
  std::array<double, 4> a = {1, 2, 3, 4};
  std::vector<double> v = {5, 4, 3, 2};
  printf("first(a) = %f, last(a) = %f\nfirst(v) = %f, last(v) = %f\nfirst(3) = "
         "%d, last(3) = %d\n",
         first(a), last(a), first(v), last(v), first(3), last(3));
  return 0;
}
```
I get
```
first(a) = 1.000000, last(a) = 4.000000
first(v) = 5.000000, last(v) = 2.000000
first(3) = 3, last(3) = 3
```

Thus, rather than defining abstract type trees like in Julia, we can dispatch based on a type's capabilities.
Type trees are not always sufficient. Things like `Base.StridedArray` are messy
```julia
julia> Base.StridedArray
StridedArray (alias for Union{DenseArray{T, N}, Base.ReinterpretArray{T, N, S, A, IsReshaped} where {A<:Union{SubArray{T, N, A, I, true} where {T, N, A<:DenseArray, I<:Union{Tuple{Vararg{Real}}, Tuple{AbstractUnitRange, Vararg{Any}}}}, DenseArray}, IsReshaped, S}, Base.ReshapedArray{T, N, A} where A<:Union{Base.ReinterpretArray{T, N, S, A, IsReshaped} where {T, N, A<:Union{SubArray{T, N, A, I, true} where {T, N, A<:DenseArray, I<:Union{Tuple{Vararg{Real}}, Tuple{AbstractUnitRange, Vararg{Any}}}}, DenseArray}, IsReshaped, S}, SubArray{T, N, A, I, true} where {T, N, A<:DenseArray, I<:Union{Tuple{Vararg{Real}}, Tuple{AbstractUnitRange, Vararg{Any}}}}, DenseArray}, SubArray{T, N, A, I} where {A<:Union{Base.ReinterpretArray{T, N, S, A, IsReshaped} where {T, N, A<:Union{SubArray{T, N, A, I, true} where {T, N, A<:DenseArray, I<:Union{Tuple{Vararg{Real}}, Tuple{AbstractUnitRange, Vararg{Any}}}}, DenseArray}, IsReshaped, S}, Base.ReshapedArray{T, N, A} where {T, N, A<:Union{Base.ReinterpretArray{T, N, S, A, IsReshaped} where {T, N, A<:Union{SubArray{T, N, A, I, true} where {T, N, A<:DenseArray, I<:Union{Tuple{Vararg{Real}}, Tuple{AbstractUnitRange, Vararg{Any}}}}, DenseArray}, IsReshaped, S}, SubArray{T, N, A, I, true} where {T, N, A<:DenseArray, I<:Union{Tuple{Vararg{Real}}, Tuple{AbstractUnitRange, Vararg{Any}}}}, DenseArray}}, DenseArray}, I<:Tuple{Vararg{Union{Base.AbstractCartesianIndex, Base.ReshapedArray{T, N, A, Tuple{}} where {T, N, A<:AbstractUnitRange}, Union{AbstractRange{<:Union{Int128, Int16, Int32, Int64, Int8, UInt128, UInt16, UInt32, UInt64, UInt8}}, var"#s93"} where var"#s93"<:Union{Int128, Int16, Int32, Int64, Int8, UInt128, UInt16, UInt32, UInt64, UInt8}}}}}} where {T, N})
```
and excludes actually strided arrays like `StaticArrays.MArray` that could dispatch on `BLAS` and `LAPACK` routines, while including arrays like [GPUArrays](https://github.com/JuliaGPU/GPUArrays.jl/blob/cd237a4f77ddfaeea6d01c6ea84ca02e71da3108/src/device/abstractarray.jl#L15) that can't.
Not to mention the mess around what is a mutable array.

```c++
#include <array>
#include <cstddef>
#include <cstdio>
#include <vector>

template <typename T, typename S>
concept CanAssign = requires(T t, S s, size_t i) {
   { t[i] = s } -> std::convertible_to<S>;
};

static_assert(CanAssign<std::vector<double>, double>);
static_assert(CanAssign<std::vector<double>, int>);
static_assert(!CanAssign<std::vector<double>, std::vector<double>>);

template <typename T>
void setInd(CanAssign<T> auto *v, size_t i, T t) {
  (*v)[i] = t;
}
void setInd(const auto *v, size_t i, const auto &t) {
  printf("oops! Should implement an out of place operator?\n");
}

int main() {
  std::vector<double> v = {5, 4, 3, 2};
  setInd(&v, 1, 1.0);
  printf("v[1] = %f\n", v[1]);
  setInd(&v, 1, v);
  const std::vector<double> &vref = v;
  setInd(&v, 1, v);
  return 0;
}
```
I get
```
v[1] = 1.000000
oops! Should implement an out of place operator?
oops! Should implement an out of place operator?
```
So the first `setInd` worked, but setting an illegal type (we can't store a `vector<double>` into a `vector<double>`!), or setting into a constant array of course do not!

