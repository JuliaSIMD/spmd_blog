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

With help from [this question](https://discourse.llvm.org/t/how-to-build-libc-with-pstl-support/69341), I landed on this cmake config (NOTE: if you are not on Apple, do not run this before applying the patch presented below)
```sh
$ git clone https://github.com/llvm/llvm-project.git
$ mkdir llvm-project/builddir && llvm-project/builddir
$ cmake -G "Ninja" \
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
      -DLIBUNWIND_USE_COMPILER_RT=YES \
      -DBOOTSTRAP_CMAKE_BUILD_TYPE=Release \
      -DBOOTSTRAP_LLVM_ENABLE_PROJECTS="clang;lld;lldb" \
      -DBOOTSTRAP_LLVM_ENABLE_RUNTIMES="compiler-rt;libc;libcxx;libcxxabi;libunwind" \
      -DBOOTSTRAP_CLANG_DEFAULT_CXX_STDLIB=libc++ \
      -DBOOTSTRAP_CLANG_DEFAULT_RTLIB=compiler-rt \
      -DBOOTSTRAP_LIBCXX_USE_COMPILER_RT=YES \
      -DBOOTSTRAP_LIBCXXABI_USE_COMPILER_RT=YES \
      -DBOOTSTRAP_LIBCXXABI_USE_LLVM_UNWINDER=YES \
      -DBOOTSTRAP_LIBUNWIND_USE_COMPILER_RT=YES \
      -DBOOTSTRAP_LLVM_USE_LINKER=lld \
      ../llvm
$ cmake --build .
$ cmake --build . --target runtimes
$ cmake --build . --target install
```
This would produce an llvm binary that would install in `~/.local/stow/llvm`. You can then run `cd ~/.local/stow && sudo stow llvm -t /usr/local` to install it.
Alternatively, setting `~/.local` as the prefix and setting your paths may work just as well; I still had to run
```sh
sudo ldconfig /usr/local/lib/x86_64-unknown-linux-gnu/
```
The install target also performed compilation, so I really wanted to avoid using `sudo` for it.

With this, I had an installed `clang++` binary on my path that could find its headers and `libc++`, which it used by default. This could successfully compile my test C++20 and C++23 projects (that were restricted to [the features clang 17 supports](https://en.cppreference.com/w/cpp/23)) when not using sanitizers.

However, it failed to find the sanitizer runtime libraries, because they were not there!
They should be a part of `compiler-rt`, and are built as part of `compiler-rt` when including it among the `*LLVM_ENABLE_PROJECTS` instead of `*LLVM_ENABLE_RUNTIMES`.
However, I didn't try the minimal diff from the above cmake config to see if this would still produce a working clang. Starting from a more minimal diff, clang failed to find include files.

I didn't want to just make random changes without understanding what is going on and recompiling in hopes of fixing the problem, so I decided to dig into the `CMakeLists.txt` files to look at options.
I debugging the `CMakeLists.txt` via inserting tons of `message` commands of the form:
```cmake
message(STATUS, "some_variable = ${some_variable}")
message(FATAL_ERROR "some_other_variable = ${some_other_variable}")
```
and rerunning the big cmake configuration command (which I placed in a script). The messages let me see the values of variables at different points, and whether certain files were being included.
This print debugging can lead you through the call graph, eventually noticing that a file is being included that is meant to define the valid sanitizers, but ended up defining them all as invalid due to not detecting the host's CPU architecture. This bug is on the non-Apple path of the control flow, so as a work around I added a `detect_target_arch()` call at the top of the `else()`/non-Apple branch ahead of target filtering.
Workaround only, probably not the correct fix -- would need to double check if we reach this code when building compiler-rt as a project, and if so where the architecture was detected to determine if there should be a change there to support it when building the rt as a runtime.
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

With this, and a `ninja runtimes`, I finally built the sanitizers.

My project still wouldn't compile with sanitizers, because CMake said clang couldn't compile a simple test program, getting undefined references to the unwind library. I'm not sure what the best way to address this is, but I simply added 
```cmake
if(((USE_SANITIZER MATCHES "([Aa]ddress)") OR (USE_SANITIZER MATCHES "([Aa]ddress);([Uu]ndefined)"
                                               )
) AND (CMAKE_CXX_COMPILER_ID MATCHES "Clang"))
  set(CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} -lunwind -Wno-unused-command-line-argument")
endif()
```
to link `lunwind`. This had to be added to the top of the `CMakeLists.txt`, so that it could apply to all dependencies. This also required the `-Wno-unused-command-line-argument` argument to avoid errors from unused command line arguments; apparently we only need to link libunwind occasionally. It'd be great if Clang could do so automatically when needed and we didn't need to suppress a warning. It'd also be great if we could use better hygiene in our CMakeLists, but I fear we need to apply this to all targets as a workaround for this failure to bring in `libunwind` automatically (which I'd like to call a bug, but perhaps it is only behavior I dislike and is working as intended?).

With this, I could finally compile and run using asan + ubsan.

