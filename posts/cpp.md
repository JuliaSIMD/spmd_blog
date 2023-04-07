+++
title = "C++"
hascode = true
date = Date(2023, 4, 06)
rss = "c++ intro"
+++

@def tags = ["c++", "intro"]

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


