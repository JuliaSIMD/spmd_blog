+++
title = "Building an LLVM Toolchain"
hascode = true
date = Date(2023, 9, 24)
rss = "building llvm toolchain"
+++

@def tags = ["c++", "llvm"]

### Intro

If you want to use the latest Julia release, it is easy to install it via [juliaup](https://github.com/JuliaLang/juliaup), [download a binary](https://julialang.org/downloads/), or simply check out the git repo and `make` -- I've rarely had trouble building it from source.
Then, all you need to do is use this binary.
You might need to re `instantiate` a few projects, but that's about all it takes.
`juliaup` lets you choose which `julia` to launch, otherwise you may want to set aliases, paths, symlinks or whatever your preferred method is of making it your default `julia`.

Open source Julia projects only supporting the latest release is common, because each new Julia release is much better than what came before it, why wouldn't you be on the latest release?
Only if you place extreme value on stability, in which case staying on the latest package releases probably isn't for you either.

Each new release of Clang and GCC are also major improvements over old ones, especially when it comes to supporting the latest C++ standards, which offer tons of niceties. Only the insane would use SFINAE when C++20 concepts are available.

Yet, upgrading to the latest standards and toolchains seems comparatively slow in the world of C++.
I suspect this is because setting up and using a custom toolchain with a new standard library is comparatively difficult.

### Three approaches

1. Use [Fedora](https://fedoraproject.org/). Fedora releases twice a year, and these releases coincide fairly closely with the GCC & even-LLVM version releases in spring, and then the odd-LLVM release in fall. If you're impatient, you can even checkout the beta releases of Fedora, ahead of a release.
Most other popular distributions, e.g. Arch Linux and Ubuntu, are slower to upgrade, e.g. Arch Linux finally upgraded to LLVM 16 in September, about 5 months after you could've gotten your hands on it with Fedora.

The major downside is that this requires changing the operating system, and telling others that building your project requires changing operating systems will not go over well.

2. Build from source or [download a release](https://releases.llvm.org/) and set paths. You need to set paths for the binaries, includes, and libraries. Preferably in `.profile`, or at some point early enough that any tooling, e.g. `clangd`, is also using your fresh version and correctly finds the new headers. This is probably the best approach.

3. Build a compiler that knows how to find the rest of the toolchain you built it with. Coming from more modern languages, this seems like the obvious approach, or how things should just work. It'd be ridiculous if any Julia not installed through your operating system's blessed package manager didn't know how to find the standard libraries it shipped with. But, if you simply build clang, it won't be able to find any headers at all, even when `/usr/bin/clang++` and `/usr/bin/g++` just work.

### Configuring LLVM

With help from [this question](https://discourse.llvm.org/t/how-to-build-libc-with-pstl-support/69341), I landed on this cmake config
```sh
cmake -G "Ninja" \
      -DCMAKE_C_COMPILER=/usr/bin/clang \
      -DCMAKE_CXX_COMPILER=/usr/bin/clang++ \
      -DLLVM_USE_LINKER=lld \
      -DLLVM_ENABLE_PROJECTS="clang;clang-tools-extra;lld;lldb" \
      -DLLVM_TARGETS_TO_BUILD=X86 \
      -DCLANG_ENABLE_BOOTSTRAP=ON \
      -DCMAKE_BUILD_TYPE=Release \
      -DLLVM_ENABLE_ASSERTIONS=OFF \
      -DLLVM_BUILD_LLVM_DYLIB=ON \
      -DLLVM_LINK_LLVM_DYLIB=ON \
      -DLLVM_INSTALL_UTILS=ON \
      -DCMAKE_INSTALL_PREFIX=~/.local/stow/llvm \
      -DCLANG_DEFAULT_CXX_STDLIB=libc++ \
      -DCLANG_DEFAULT_RTLIB=compiler-rt \
      -DLLVM_ENABLE_RUNTIMES="compiler-rt;libcxx;libcxxabi;libunwind" \
      -DLIBCXX_USE_COMPILER_RT=YES \
      -DLIBCXXABI_USE_COMPILER_RT=YES \
      -DLIBCXXABI_USE_LLVM_UNWINDER=YES \
      -DLIBUNWIND_USE_COMPILER_RT=Yes \
      -DBOOTSTRAP_CMAKE_BUILD_TYPE=Release \
      -DBOOTSTRAP_LLVM_ENABLE_PROJECTS="clang;lld;lldb" \
      -DBOOTSTRAP_LLVM_ENABLE_RUNTIMES="compiler-rt;libc;libcxx;libcxxabi;libunwind" \
      -DBOOTSTRAP_CLANG_DEFAULT_CXX_STDLIB=libc++ \
      -DBOOTSTRAP_CLANG_DEFAULT_RTLIB=compiler-rt \
      -DBOOTSTRAP_LIBCXX_USE_COMPILER_RT=YES \
      -DBOOTSTRAP_LIBCXXABI_USE_COMPILER_RT=YES \
      -DBOOTSTRAP_LIBCXXABI_USE_LLVM_UNWINDER=YES \
      -DBOOTSTRAP_LIBUNWIND_USE_COMPILER_RT=Yes \
      -DBOOTSTRAP_LLVM_USE_LINKER=lld \
      ../llvm
```
This would produce an llvm binary that would install in `~/.local/stow/llvm`.

```diff
diff --git a/compiler-rt/cmake/config-ix.cmake b/compiler-rt/cmake/config-ix.cmake
index 8d3dc8d208b2..c6daa35f9622 100644
--- a/compiler-rt/cmake/config-ix.cmake
+++ b/compiler-rt/cmake/config-ix.cmake
@@ -642,6 +642,7 @@ if(APPLE)
     SANITIZER_COMMON_SUPPORTED_ARCH)
 
 else()
+  detect_target_arch()
   # Architectures supported by compiler-rt libraries.
   filter_available_targets(SANITIZER_COMMON_SUPPORTED_ARCH
     ${ALL_SANITIZER_COMMON_SUPPORTED_ARCH})
```




