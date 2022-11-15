+++
title = "Dependency Analysis of Julia Loops Example"
hascode = true
date = Date(2022, 11, 13)
rss = "Example of running dependence analysis on Julia loops"
+++

@def tags = ["dependency analysis", "llvm-ir"]

In order to try examples of LoopModels on Julia IR, the current workflow requires a custom build of Julia with vectorization disabled:
```diff
diff --git a/src/aotcompile.cpp b/src/aotcompile.cpp
index 98777ddd17..7da2f53bd5 100644
--- a/src/aotcompile.cpp
+++ b/src/aotcompile.cpp
@@ -867,14 +867,14 @@ void addOptimizationPasses(legacy::PassManagerBase *PM, int opt_level,
     PM->add(createAllocOptPass());
     PM->add(createLoopDeletionPass());
     PM->add(createInstructionCombiningPass());
-    PM->add(createLoopVectorizePass());
+    // PM->add(createLoopVectorizePass());
     PM->add(createLoopLoadEliminationPass());
     // Cleanup after LV pass
     PM->add(createInstructionCombiningPass());
     PM->add(createCFGSimplificationPass( // Aggressive CFG simplification
         aggressiveSimplifyCFGOptions
     ));
-    PM->add(createSLPVectorizerPass());
+    // PM->add(createSLPVectorizerPass());
     // might need this after LLVM 11:
     //PM->add(createVectorCombinePass());
```

Now, we can define a function to run Julia's optimization passes (minus the  vectorization passes, as we've disabled them), an then a series of additional passes (in this case, `loop-simplify` followed by loop-invariant code motion (`licm`)):
```julia
function write_code_llvm(f, args...; filename = nothing, passes = "function(simplifycfg,early-cse,loop-simplify)")
    path, io = mktemp()
    T = args isa Tuple{Vararg{DataType}} ? args : map(typeof,args)
    code_llvm(io, f, T;
        raw = true,
        dump_module = true,
        optimize = true,
        debuginfo = :none)
    close(io)
    newfilename = filename === nothing ? path * ".ll" : filename
    mv(path, newfilename, force=true)
    run(pipeline(`opt --passes="$passes" $newfilename`, `opt -licm -S -o $newfilename`))
    return newfilename
end
```
This function accepts either types or instances (in which case it calls `typeof`) as arguments.
Currently, I'm not running `loop-rotate` or `lcssa`. Finding the set of passes to actually run will take some work.

Now we can create a `.ll` file, e.g.:
```julia-repl
julia> function triangular_solve0!(A,B,U)
           M,N = size(A)
           @assert M == size(B,1)
           @assert N == size(B,2)
           @assert N == LinearAlgebra.checksquare(U)
           @inbounds for m = 1:M
               for n = 1:N
                   A[m,n] = B[m,n]
               end
               for n = 1:N
                   A[m,n] /= U[n,n]
                   for k = n+1:N
                       A[m,k] -= A[m,n]*U[n,k]
                   end
               end
           end
       end
triangular_solve0! (generic function with 1 method)

julia> B = rand(5,4); U = UpperTriangular(rand(4,4)); A = similar(B);

julia> triangular_solve0!(A, B, parent(U)); A ≈ B / U
true

julia> write_code_llvm(triangular_solve0!, (Matrix{Float64} for _ in 1:3)..., filename = "/home/chriselrod/Documents/progwork/cxx/LoopPlayground/LoopInductTests/test/triangular_solve.ll")
"/home/chriselrod/Documents/progwork/cxx/LoopPlayground/LoopInductTests/test/triangular_solve.ll"
```

For reference, here is the resulting function:
```llvm
define nonnull {} addrspace(10)* @"japi1_triangular_solve0!_699"({} addrspace(10)* %0, {} addrspace(10)** noalias nocapture noundef readonly %1, i32 %2) #0 !dbg !4 {
top:
  %3 = alloca [2 x {} addrspace(10)*], align 8
  %gcframe74 = alloca [3 x {} addrspace(10)*], align 16
  %gcframe74.sub = getelementptr inbounds [3 x {} addrspace(10)*], [3 x {} addrspace(10)*]* %gcframe74, i64 0, i64 0
  %.sub = getelementptr inbounds [2 x {} addrspace(10)*], [2 x {} addrspace(10)*]* %3, i64 0, i64 0
  %4 = bitcast [3 x {} addrspace(10)*]* %gcframe74 to i8*
  call void @llvm.memset.p0i8.i32(i8* noundef nonnull align 16 dereferenceable(24) %4, i8 0, i32 24, i1 false), !tbaa !7
  %5 = alloca {} addrspace(10)**, align 8
  store volatile {} addrspace(10)** %1, {} addrspace(10)*** %5, align 8
  %thread_ptr = call i8* asm "movq %fs:0, $0", "=r"() #8
  %ppgcstack_i8 = getelementptr i8, i8* %thread_ptr, i64 -8
  %ppgcstack = bitcast i8* %ppgcstack_i8 to {}****
  %pgcstack = load {}***, {}**** %ppgcstack, align 8
  %6 = bitcast [3 x {} addrspace(10)*]* %gcframe74 to i64*
  store i64 4, i64* %6, align 16, !tbaa !7
  %7 = getelementptr inbounds [3 x {} addrspace(10)*], [3 x {} addrspace(10)*]* %gcframe74, i64 0, i64 1
  %8 = bitcast {} addrspace(10)** %7 to {}***
  %9 = load {}**, {}*** %pgcstack, align 8
  store {}** %9, {}*** %8, align 8, !tbaa !7
  %10 = bitcast {}*** %pgcstack to {} addrspace(10)***
  store {} addrspace(10)** %gcframe74.sub, {} addrspace(10)*** %10, align 8
  %11 = load {} addrspace(10)*, {} addrspace(10)** %1, align 8, !tbaa !11, !nonnull !6, !dereferenceable !13, !align !14
  %12 = getelementptr inbounds {} addrspace(10)*, {} addrspace(10)** %1, i64 1
  %13 = load {} addrspace(10)*, {} addrspace(10)** %12, align 8, !tbaa !11, !nonnull !6, !dereferenceable !13, !align !14
  %14 = getelementptr inbounds {} addrspace(10)*, {} addrspace(10)** %1, i64 2
  %15 = load {} addrspace(10)*, {} addrspace(10)** %14, align 8, !tbaa !11, !nonnull !6, !dereferenceable !13, !align !14
  %16 = bitcast {} addrspace(10)* %11 to {} addrspace(10)* addrspace(10)*, !dbg !15
  %17 = addrspacecast {} addrspace(10)* addrspace(10)* %16 to {} addrspace(10)* addrspace(11)*, !dbg !15
  %18 = getelementptr inbounds {} addrspace(10)*, {} addrspace(10)* addrspace(11)* %17, i64 3, !dbg !15
  %19 = bitcast {} addrspace(10)* addrspace(11)* %18 to i64 addrspace(11)*, !dbg !15
  %20 = load i64, i64 addrspace(11)* %19, align 8, !dbg !15, !tbaa !11, !range !19
  %21 = getelementptr inbounds {} addrspace(10)*, {} addrspace(10)* addrspace(11)* %17, i64 4, !dbg !15
  %22 = bitcast {} addrspace(10)* addrspace(11)* %21 to i64 addrspace(11)*, !dbg !15
  %23 = load i64, i64 addrspace(11)* %22, align 8, !dbg !15, !tbaa !11, !range !19
  %24 = bitcast {} addrspace(10)* %13 to {} addrspace(10)* addrspace(10)*, !dbg !20
  %25 = addrspacecast {} addrspace(10)* addrspace(10)* %24 to {} addrspace(10)* addrspace(11)*, !dbg !20
  %26 = getelementptr inbounds {} addrspace(10)*, {} addrspace(10)* addrspace(11)* %25, i64 3, !dbg !20
  %27 = bitcast {} addrspace(10)* addrspace(11)* %26 to i64 addrspace(11)*, !dbg !20
  %28 = load i64, i64 addrspace(11)* %27, align 8, !dbg !20, !tbaa !11, !range !19
  %.not = icmp eq i64 %20, %28, !dbg !22
  br i1 %.not, label %L6, label %L163, !dbg !21

L6:                                               ; preds = %top
  %29 = getelementptr inbounds {} addrspace(10)*, {} addrspace(10)* addrspace(11)* %25, i64 4, !dbg !25
  %30 = bitcast {} addrspace(10)* addrspace(11)* %29 to i64 addrspace(11)*, !dbg !25
  %31 = load i64, i64 addrspace(11)* %30, align 8, !dbg !25, !tbaa !11, !range !19
  %.not47 = icmp eq i64 %23, %31, !dbg !27
  br i1 %.not47, label %L10, label %L160, !dbg !26

L10:                                              ; preds = %L6
  %32 = bitcast {} addrspace(10)* %15 to {} addrspace(10)* addrspace(10)*, !dbg !28
  %33 = addrspacecast {} addrspace(10)* addrspace(10)* %32 to {} addrspace(10)* addrspace(11)*, !dbg !28
  %34 = getelementptr inbounds {} addrspace(10)*, {} addrspace(10)* addrspace(11)* %33, i64 3, !dbg !28
  %35 = bitcast {} addrspace(10)* addrspace(11)* %34 to i64 addrspace(11)*, !dbg !28
  %36 = load i64, i64 addrspace(11)* %35, align 8, !dbg !28, !tbaa !11, !range !19
  %37 = getelementptr inbounds {} addrspace(10)*, {} addrspace(10)* addrspace(11)* %33, i64 4, !dbg !28
  %38 = bitcast {} addrspace(10)* addrspace(11)* %37 to i64 addrspace(11)*, !dbg !28
  %39 = load i64, i64 addrspace(11)* %38, align 8, !dbg !28, !tbaa !11, !range !19
  %.not48 = icmp eq i64 %36, %39, !dbg !33
  br i1 %.not48, label %L23, label %L16, !dbg !34

L16:                                              ; preds = %L10
  %ptls_field81 = getelementptr inbounds {}**, {}*** %pgcstack, i64 2, !dbg !35
  %40 = bitcast {}*** %ptls_field81 to i8**, !dbg !35
  %ptls_load8283 = load i8*, i8** %40, align 8, !dbg !35, !tbaa !7
  %41 = call noalias nonnull {} addrspace(10)* @ijl_gc_pool_alloc(i8* %ptls_load8283, i32 1440, i32 32) #3, !dbg !35
  %42 = bitcast {} addrspace(10)* %41 to i64 addrspace(10)*, !dbg !35
  %43 = getelementptr inbounds i64, i64 addrspace(10)* %42, i64 -1, !dbg !35
  store atomic i64 139879779354576, i64 addrspace(10)* %43 unordered, align 8, !dbg !35, !tbaa !38
  %44 = bitcast {} addrspace(10)* %41 to i8 addrspace(10)*, !dbg !35
  store i64 %36, i64 addrspace(10)* %42, align 8, !dbg !35, !tbaa !41
  %.sroa.2.0..sroa_idx = getelementptr inbounds i8, i8 addrspace(10)* %44, i64 8, !dbg !35
  %.sroa.2.0..sroa_cast = bitcast i8 addrspace(10)* %.sroa.2.0..sroa_idx to i64 addrspace(10)*, !dbg !35
  store i64 %39, i64 addrspace(10)* %.sroa.2.0..sroa_cast, align 8, !dbg !35, !tbaa !41
  %45 = getelementptr inbounds [3 x {} addrspace(10)*], [3 x {} addrspace(10)*]* %gcframe74, i64 0, i64 2
  store {} addrspace(10)* %41, {} addrspace(10)** %45, align 16
  store {} addrspace(10)* addrspacecast ({}* inttoptr (i64 139879848027248 to {}*) to {} addrspace(10)*), {} addrspace(10)** %.sub, align 8, !dbg !35
  %46 = getelementptr inbounds [2 x {} addrspace(10)*], [2 x {} addrspace(10)*]* %3, i64 0, i64 1, !dbg !35
  store {} addrspace(10)* %41, {} addrspace(10)** %46, align 8, !dbg !35
  %47 = call nonnull {} addrspace(10)* @j1_print_to_string_700({} addrspace(10)* addrspacecast ({}* inttoptr (i64 139879781117552 to {}*) to {} addrspace(10)*), {} addrspace(10)** nonnull %.sub, i32 2), !dbg !35
  store {} addrspace(10)* %47, {} addrspace(10)** %45, align 16
  store {} addrspace(10)* %47, {} addrspace(10)** %.sub, align 8, !dbg !34
  %48 = call nonnull {} addrspace(10)* @ijl_apply_generic({} addrspace(10)* addrspacecast ({}* inttoptr (i64 139879785329264 to {}*) to {} addrspace(10)*), {} addrspace(10)** nonnull %.sub, i32 1), !dbg !34
  %49 = addrspacecast {} addrspace(10)* %48 to {} addrspace(12)*, !dbg !34
  call void @ijl_throw({} addrspace(12)* %49), !dbg !34
  unreachable, !dbg !34

L23:                                              ; preds = %L10
  %.not49 = icmp eq i64 %23, %36, !dbg !42
  br i1 %.not49, label %L25, label %L157, !dbg !32

L25:                                              ; preds = %L23
  %.not51.not = icmp eq i64 %20, 0, !dbg !43
  br i1 %.not51.not, label %L156, label %L42.preheader, !dbg !54

L42.preheader:                                    ; preds = %L25
  %cond = icmp eq i64 %23, 0
  %50 = bitcast {} addrspace(10)* %13 to double addrspace(13)* addrspace(10)*
  %51 = addrspacecast double addrspace(13)* addrspace(10)* %50 to double addrspace(13)* addrspace(11)*
  %52 = load double addrspace(13)*, double addrspace(13)* addrspace(11)* %51, align 8
  %53 = bitcast {} addrspace(10)* %11 to double addrspace(13)* addrspace(10)*
  %54 = addrspacecast double addrspace(13)* addrspace(10)* %53 to double addrspace(13)* addrspace(11)*
  %55 = load double addrspace(13)*, double addrspace(13)* addrspace(11)* %54, align 8
  %56 = bitcast {} addrspace(10)* %15 to double addrspace(13)* addrspace(10)*
  %57 = addrspacecast double addrspace(13)* addrspace(10)* %56 to double addrspace(13)* addrspace(11)*
  %58 = load double addrspace(13)*, double addrspace(13)* addrspace(11)* %57, align 8
  br i1 %cond, label %L156, label %L60.preheader.preheader, !dbg !55

L60.preheader.preheader:                          ; preds = %L42.preheader
  br label %L60.preheader, !dbg !56

L60.preheader:                                    ; preds = %L60.preheader.preheader, %L145
  %value_phi4 = phi i64 [ %91, %L145 ], [ 1, %L60.preheader.preheader ]
  %59 = add nsw i64 %value_phi4, -1
  br label %L60, !dbg !56

L60:                                              ; preds = %L60, %L60.preheader
  %value_phi10 = phi i64 [ %66, %L60 ], [ 1, %L60.preheader ]
  %60 = add nsw i64 %value_phi10, -1, !dbg !57
  %61 = mul i64 %60, %20, !dbg !57
  %62 = add i64 %59, %61, !dbg !57
  %63 = getelementptr inbounds double, double addrspace(13)* %52, i64 %62, !dbg !57
  %64 = load double, double addrspace(13)* %63, align 8, !dbg !57, !tbaa !61
  %65 = getelementptr inbounds double, double addrspace(13)* %55, i64 %62, !dbg !63
  store double %64, double addrspace(13)* %65, align 8, !dbg !63, !tbaa !61
  %.not54.not = icmp eq i64 %value_phi10, %23, !dbg !65
  %66 = add nuw nsw i64 %value_phi10, 1, !dbg !66
  br i1 %.not54.not, label %L91.preheader, label %L60, !dbg !56

L91.preheader:                                    ; preds = %L60
  br label %L91, !dbg !67

L91:                                              ; preds = %L91.preheader, %L134
  %value_phi19 = phi i64 [ %77, %L134 ], [ 1, %L91.preheader ]
  %67 = add nsw i64 %value_phi19, -1, !dbg !68
  %68 = mul i64 %67, %20, !dbg !68
  %69 = add i64 %59, %68, !dbg !68
  %70 = getelementptr inbounds double, double addrspace(13)* %55, i64 %69, !dbg !68
  %71 = load double, double addrspace(13)* %70, align 8, !dbg !68, !tbaa !61
  %72 = mul i64 %23, %67, !dbg !68
  %73 = add i64 %67, %72, !dbg !68
  %74 = getelementptr inbounds double, double addrspace(13)* %58, i64 %73, !dbg !68
  %75 = load double, double addrspace(13)* %74, align 8, !dbg !68, !tbaa !61
  %76 = fdiv double %71, %75, !dbg !70
  store double %76, double addrspace(13)* %70, align 8, !dbg !73, !tbaa !61
  %77 = add nuw nsw i64 %value_phi19, 1, !dbg !74
  %.not57.not = icmp ult i64 %value_phi19, %23, !dbg !76
  %value_phi21 = select i1 %.not57.not, i64 %23, i64 %value_phi19, !dbg !80
  %.not58.not.not = icmp sgt i64 %value_phi21, %value_phi19, !dbg !86
  br i1 %.not58.not.not, label %L115.preheader, label %L134, !dbg !67

L115.preheader:                                   ; preds = %L91
  br label %L115, !dbg !90

L115:                                             ; preds = %L115.preheader, %L115.L115_crit_edge
  %78 = phi double [ %.pre, %L115.L115_crit_edge ], [ %76, %L115.preheader ], !dbg !91
  %value_phi25 = phi i64 [ %90, %L115.L115_crit_edge ], [ %77, %L115.preheader ]
  %79 = add i64 %value_phi25, -1, !dbg !91
  %80 = mul i64 %79, %20, !dbg !91
  %81 = add i64 %59, %80, !dbg !91
  %82 = getelementptr inbounds double, double addrspace(13)* %55, i64 %81, !dbg !91
  %83 = load double, double addrspace(13)* %82, align 8, !dbg !91, !tbaa !61
  %84 = mul i64 %23, %79, !dbg !91
  %85 = add i64 %67, %84, !dbg !91
  %86 = getelementptr inbounds double, double addrspace(13)* %58, i64 %85, !dbg !91
  %87 = load double, double addrspace(13)* %86, align 8, !dbg !91, !tbaa !61
  %88 = fmul double %78, %87, !dbg !93
  %89 = fsub double %83, %88, !dbg !95
  store double %89, double addrspace(13)* %82, align 8, !dbg !97, !tbaa !61
  %.not59.not = icmp eq i64 %value_phi25, %value_phi21, !dbg !98
  br i1 %.not59.not, label %L134.loopexit, label %L115.L115_crit_edge, !dbg !90

L115.L115_crit_edge:                              ; preds = %L115
  %90 = add i64 %value_phi25, 1, !dbg !99
  %.pre = load double, double addrspace(13)* %70, align 8, !dbg !91, !tbaa !61
  br label %L115, !dbg !90

L134.loopexit:                                    ; preds = %L115
  br label %L134, !dbg !100

L134:                                             ; preds = %L134.loopexit, %L91
  %.not60.not = icmp eq i64 %value_phi19, %23, !dbg !100
  br i1 %.not60.not, label %L145, label %L91, !dbg !102

L145:                                             ; preds = %L134
  %.not61 = icmp eq i64 %value_phi4, %20, !dbg !103
  %91 = add nuw nsw i64 %value_phi4, 1, !dbg !104
  br i1 %.not61, label %L156.loopexit, label %L60.preheader, !dbg !105

L156.loopexit:                                    ; preds = %L145
  br label %L156

L156:                                             ; preds = %L156.loopexit, %L42.preheader, %L25
  %92 = load {} addrspace(10)*, {} addrspace(10)** %7, align 8, !tbaa !7
  %93 = bitcast {}*** %pgcstack to {} addrspace(10)**
  store {} addrspace(10)* %92, {} addrspace(10)** %93, align 8, !tbaa !7
  ret {} addrspace(10)* addrspacecast ({}* inttoptr (i64 139880021618696 to {}*) to {} addrspace(10)*), !dbg !105

L157:                                             ; preds = %L23
  store {} addrspace(10)* addrspacecast ({}* inttoptr (i64 139879957513232 to {}*) to {} addrspace(10)*), {} addrspace(10)** %.sub, align 8, !dbg !32
  %94 = call nonnull {} addrspace(10)* @ijl_apply_generic({} addrspace(10)* addrspacecast ({}* inttoptr (i64 139879780146624 to {}*) to {} addrspace(10)*), {} addrspace(10)** nonnull %.sub, i32 1), !dbg !32
  %95 = addrspacecast {} addrspace(10)* %94 to {} addrspace(12)*, !dbg !32
  call void @ijl_throw({} addrspace(12)* %95), !dbg !32
  unreachable, !dbg !32

L160:                                             ; preds = %L6
  %ptls_field6878 = getelementptr inbounds {}**, {}*** %pgcstack, i64 2, !dbg !26
  %96 = bitcast {}*** %ptls_field6878 to i8**, !dbg !26
  %ptls_load697980 = load i8*, i8** %96, align 8, !dbg !26, !tbaa !7
  %97 = call noalias nonnull {} addrspace(10)* @ijl_gc_pool_alloc(i8* %ptls_load697980, i32 1392, i32 16) #3, !dbg !26
  %98 = bitcast {} addrspace(10)* %97 to i64 addrspace(10)*, !dbg !26
  %99 = getelementptr inbounds i64, i64 addrspace(10)* %98, i64 -1, !dbg !26
  store atomic i64 139879780146624, i64 addrspace(10)* %99 unordered, align 8, !dbg !26, !tbaa !38
  %100 = bitcast {} addrspace(10)* %97 to {} addrspace(10)* addrspace(10)*, !dbg !26
  store {} addrspace(10)* addrspacecast ({}* inttoptr (i64 139879957512400 to {}*) to {} addrspace(10)*), {} addrspace(10)* addrspace(10)* %100, align 8, !dbg !26, !tbaa !106
  %101 = addrspacecast {} addrspace(10)* %97 to {} addrspace(12)*, !dbg !26
  call void @ijl_throw({} addrspace(12)* %101), !dbg !26
  unreachable, !dbg !26

L163:                                             ; preds = %top
  %ptls_field7175 = getelementptr inbounds {}**, {}*** %pgcstack, i64 2, !dbg !21
  %102 = bitcast {}*** %ptls_field7175 to i8**, !dbg !21
  %ptls_load727677 = load i8*, i8** %102, align 8, !dbg !21, !tbaa !7
  %103 = call noalias nonnull {} addrspace(10)* @ijl_gc_pool_alloc(i8* %ptls_load727677, i32 1392, i32 16) #3, !dbg !21
  %104 = bitcast {} addrspace(10)* %103 to i64 addrspace(10)*, !dbg !21
  %105 = getelementptr inbounds i64, i64 addrspace(10)* %104, i64 -1, !dbg !21
  store atomic i64 139879780146624, i64 addrspace(10)* %105 unordered, align 8, !dbg !21, !tbaa !38
  %106 = bitcast {} addrspace(10)* %103 to {} addrspace(10)* addrspace(10)*, !dbg !21
  store {} addrspace(10)* addrspacecast ({}* inttoptr (i64 139879957462352 to {}*) to {} addrspace(10)*), {} addrspace(10)* addrspace(10)* %106, align 8, !dbg !21, !tbaa !106
  %107 = addrspacecast {} addrspace(10)* %103 to {} addrspace(12)*, !dbg !21
  call void @ijl_throw({} addrspace(12)* %107), !dbg !21
  unreachable, !dbg !21
}
```

Now you can clone LoopModels, and once you've installed all the dependencies:
```sh
CC_LD=lld CXX_LD=lld CXXFLAGS="" meson setup builddir -Db_santize=address -Db_coverage=true
cd builddir
meson compile TurboLoop
opt -mcpu=native --disable-output --load-pass-plugin=/home/chriselrod/Documents/progwork/cxx/LoopModels/builddir/libTurboLoop.so --passes="function(turbo-loop)" /home/chriselrod/Documents/progwork/cxx/LoopPlayground/LoopInductTests/test/triangular_solve.ll
```
You'll have to substitute the path to the `libTurboLoop.so` plugin and to the `triangular_solve.ll` as appropriate.
The output is currently extremely verbose, and it does not yet actually transform the code or perform cost modeling, unrolling, or vectorization analysis.
But it does perform some dependence analysis; sample output from the above:
```
LoopBlock graph (#nodes = 3):
v_0:
mem =
  %64 = load double, double addrspace(13)* %63, align 8, !dbg !57, !tbaa !61
  store double %64, double addrspace(13)* %65, align 8, !dbg !63, !tbaa !61
inNeighbors =
outNeighbors = v_1, v_2,

v_1:
mem =
  %71 = load double, double addrspace(13)* %70, align 8, !dbg !68, !tbaa !61
  %75 = load double, double addrspace(13)* %74, align 8, !dbg !68, !tbaa !61
  store double %76, double addrspace(13)* %70, align 8, !dbg !73, !tbaa !61
inNeighbors = v_0, v_1, v_2,
outNeighbors = v_1, v_2,

v_2:
mem =
  store double %76, double addrspace(13)* %70, align 8, !dbg !73, !tbaa !61
  %83 = load double, double addrspace(13)* %82, align 8, !dbg !91, !tbaa !61
  %87 = load double, double addrspace(13)* %86, align 8, !dbg !91, !tbaa !61
  store double %89, double addrspace(13)* %82, align 8, !dbg !97, !tbaa !61
inNeighbors = v_0, v_1, v_2,
outNeighbors = v_1, v_2,


LoopBlock Edges (#edges = 10):
        Edge = Dependence Poly y -> x:
v_2 <= -1 + %23
v_3 <= -1 + %20
v_2 >= 0
v_3 >= 0
-v_0 + v_2 == 0
-v_1 + v_3 == 0

A =
[ -1  1  0  0  0 -1  0
  -1  0  1  0  0  0 -1
   0  0  0  0  0  1  0
   0  0  0  0  0  0  1 ]
E =
[  0  0  0  1  0 -1  0
   0  0  0  0  1  0 -1 ]
Schedule Constraints:
Simplex; tableau =
[  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
   0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
   0  0  1 -1 -1  0  0  0  0  0  0  0  0  0  0 -1  1 -1
   0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0
   0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0
   0  0  0  0  0  0  0  1  0 -1  0 -1  0  0  0  0  0  0
   0  0  0  0  0  0  0  0  1  0 -1  0 -1  0  0  0  0  0
   0  0  0 -1  0  1  0 -1  0  1  0  0  0  1  0  0  0  0
   0  0  0  0 -1  0  1  0 -1  0  1  0  0  0  1  0  0  0 ]
Bounding Constraints:
Simplex; tableau =
[  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
   0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
   0  0  1 -1 -1  0  0  0  0  0  0  0  0  0  0  1 -1 -1  0  0
   0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0 -1  0
   0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0 -1
   0  0  0  0  0  0  0  1  0 -1  0  1  0  0  0  0  0  0  0  0
   0  0  0  0  0  0  0  0  1  0 -1  0  1  0  0  0  0  0  0  0
   0  0  0 -1  0  1  0 -1  0  1  0  0  0 -1  0  0  0  0  0  0
   0  0  0  0 -1  0  1  0 -1  0  1  0  0  0 -1  0  0  0  0  0 ]
        Input:
Store:   store double %64, double addrspace(13)* %65, align 8, !dbg !63, !tbaa !61
ar.indexMatrix() =
[  1  0
   0  1 ]
ArrayReference %55 (dim = 2, num loops: 2, element size: 8):
A.numRow() = 2; A.numCol() = 2
Sizes: [ unknown, %20 ]
Subscripts: [ i_1 , i_0 ]
Schedule Omega: [ 0, 0, 1 ]
AffineLoopNest: alnb.getNumLoops() = 2
aln.getNumLoops() = 2
alnb.A =
[ -1  1  0 -1  0
  -1  0  1  0 -1
   0  0  0  1  0
   0  0  0  0  1 ]
Loop 1 lower bounds:
i_1 >= 0
Loop 1 upper bounds:
i_1 <= -1 + %23
i = 0; getNumSymbols() = 3
Loop 0 lower bounds:
i_0 >= 0
Loop 0 upper bounds:
i_0 <= -1 + %20

        Output:
Load:   %71 = load double, double addrspace(13)* %70, align 8, !dbg !68, !tbaa !61
ar.indexMatrix() =
[  1  0
   0  1 ]
ArrayReference %55 (dim = 2, num loops: 2, element size: 8):
A.numRow() = 2; A.numCol() = 2
Sizes: [ unknown, %20 ]
Subscripts: [ i_1 , i_0 ]
Schedule Omega: [ 0, 1, 0 ]
AffineLoopNest: alnb.getNumLoops() = 2
aln.getNumLoops() = 2
alnb.A =
[ -1  1  0 -1  0
  -1  0  1  0 -1
   0  0  0  1  0
   0  0  0  0  1 ]
Loop 1 lower bounds:
i_1 >= 0
Loop 1 upper bounds:
i_1 <= -1 + %23
i = 0; getNumSymbols() = 3
Loop 0 lower bounds:
i_0 >= 0
Loop 0 upper bounds:
i_0 <= -1 + %20

Schedule In:
nodeIndex = BitSet[0]; ref = ar.indexMatrix() =
[  1  0
   0  1 ]
ArrayReference %55 (dim = 2, num loops: 2, element size: 8):
A.numRow() = 2; A.numCol() = 2
Sizes: [ unknown, %20 ]
Subscripts: [ i_1 , i_0 ]
s.getPhi()
[  1  0
   0  1 ]
s.getFusionOmega() = [ 0, 0, 0 ]
s.getOffsetOmega() = [ 1, 0 ]

Schedule Out:
nodeIndex = BitSet[1]; ref = ar.indexMatrix() =
[  1  0
   0  1 ]
ArrayReference %55 (dim = 2, num loops: 2, element size: 8):
A.numRow() = 2; A.numCol() = 2
Sizes: [ unknown, %20 ]
Subscripts: [ i_1 , i_0 ]
s.getPhi()
[  1  0
   0  1 ]
s.getFusionOmega() = [ 0, 0, 0 ]
s.getOffsetOmega() = [ 1, 0 ]

        Edge = Dependence Poly y -> x:
v_2 <= -1 + %23
v_3 <= -1 + %20
v_2 >= 0
v_3 >= 0
-v_0 + v_2 == 0
-v_1 + v_3 == 0

A =
[ -1  1  0  0  0 -1  0
  -1  0  1  0  0  0 -1
   0  0  0  0  0  1  0
   0  0  0  0  0  0  1 ]
E =
[  0  0  0  1  0 -1  0
   0  0  0  0  1  0 -1 ]
Schedule Constraints:
Simplex; tableau =
[  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
   0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
   0  0  1 -1 -1  0  0  0  0  0  0  0  0  0  0 -1  1 -1
   0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0
   0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0
   0  0  0  0  0  0  0  1  0 -1  0 -1  0  0  0  0  0  0
   0  0  0  0  0  0  0  0  1  0 -1  0 -1  0  0  0  0  0
   0  0  0 -1  0  1  0 -1  0  1  0  0  0  1  0  0  0  0
   0  0  0  0 -1  0  1  0 -1  0  1  0  0  0  1  0  0  0 ]
Bounding Constraints:
Simplex; tableau =
[  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
   0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
   0  0  1 -1 -1  0  0  0  0  0  0  0  0  0  0  1 -1 -1  0  0
   0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0 -1  0
   0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0 -1
   0  0  0  0  0  0  0  1  0 -1  0  1  0  0  0  0  0  0  0  0
   0  0  0  0  0  0  0  0  1  0 -1  0  1  0  0  0  0  0  0  0
   0  0  0 -1  0  1  0 -1  0  1  0  0  0 -1  0  0  0  0  0  0
   0  0  0  0 -1  0  1  0 -1  0  1  0  0  0 -1  0  0  0  0  0 ]
        Input:
Store:   store double %64, double addrspace(13)* %65, align 8, !dbg !63, !tbaa !61
ar.indexMatrix() =
[  1  0
   0  1 ]
ArrayReference %55 (dim = 2, num loops: 2, element size: 8):
A.numRow() = 2; A.numCol() = 2
Sizes: [ unknown, %20 ]
Subscripts: [ i_1 , i_0 ]
Schedule Omega: [ 0, 0, 1 ]
AffineLoopNest: alnb.getNumLoops() = 2
aln.getNumLoops() = 2
alnb.A =
[ -1  1  0 -1  0
  -1  0  1  0 -1
   0  0  0  1  0
   0  0  0  0  1 ]
Loop 1 lower bounds:
i_1 >= 0
Loop 1 upper bounds:
i_1 <= -1 + %23
i = 0; getNumSymbols() = 3
Loop 0 lower bounds:
i_0 >= 0
Loop 0 upper bounds:
i_0 <= -1 + %20

        Output:
Store:   store double %76, double addrspace(13)* %70, align 8, !dbg !73, !tbaa !61
ar.indexMatrix() =
[  1  0
   0  1 ]
ArrayReference %55 (dim = 2, num loops: 2, element size: 8):
A.numRow() = 2; A.numCol() = 2
Sizes: [ unknown, %20 ]
Subscripts: [ i_1 , i_0 ]
Schedule Omega: [ 0, 1, 2 ]
AffineLoopNest: alnb.getNumLoops() = 2
aln.getNumLoops() = 2
alnb.A =
[ -1  1  0 -1  0
  -1  0  1  0 -1
   0  0  0  1  0
   0  0  0  0  1 ]
Loop 1 lower bounds:
i_1 >= 0
Loop 1 upper bounds:
i_1 <= -1 + %23
i = 0; getNumSymbols() = 3
Loop 0 lower bounds:
i_0 >= 0
Loop 0 upper bounds:
i_0 <= -1 + %20

Schedule In:
nodeIndex = BitSet[0]; ref = ar.indexMatrix() =
[  1  0
   0  1 ]
ArrayReference %55 (dim = 2, num loops: 2, element size: 8):
A.numRow() = 2; A.numCol() = 2
Sizes: [ unknown, %20 ]
Subscripts: [ i_1 , i_0 ]
s.getPhi()
[  1  0
   0  1 ]
s.getFusionOmega() = [ 0, 0, 0 ]
s.getOffsetOmega() = [ 1, 0 ]

Schedule Out:
nodeIndex = BitSet[1, 2]; ref = ar.indexMatrix() =
[  1  0
   0  1 ]
ArrayReference %55 (dim = 2, num loops: 2, element size: 8):
A.numRow() = 2; A.numCol() = 2
Sizes: [ unknown, %20 ]
Subscripts: [ i_1 , i_0 ]
s.getPhi()
[  1  0
   0  1 ]
s.getFusionOmega() = [ 0, 0, 0 ]
s.getOffsetOmega() = [ 1, 0 ]

Schedule Out:
nodeIndex = BitSet[1, 2]; ref = ar.indexMatrix() =
[  1  0
   0  1 ]
ArrayReference %55 (dim = 2, num loops: 2, element size: 8):
A.numRow() = 2; A.numCol() = 2
Sizes: [ unknown, %20 ]
Subscripts: [ i_1 , i_0 ]
s.getPhi()
[  1  1  0
   0  0  1
   0  1  0 ]
s.getFusionOmega() = [ 0, 0, 0, 0 ]
s.getOffsetOmega() = [ 0, 0, 0 ]

        Edge = Dependence Poly y -> x:
v_2 <= -1 + %23
v_3 <= -1 + %20
v_2 >= 0
v_3 >= 0
-v_0 + v_2 == 0
-v_1 + v_3 == 0

A =
[ -1  1  0  0  0 -1  0
  -1  0  1  0  0  0 -1
   0  0  0  0  0  1  0
   0  0  0  0  0  0  1 ]
E =
[  0  0  0  1  0 -1  0
   0  0  0  0  1  0 -1 ]
Schedule Constraints:
Simplex; tableau =
[  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
   0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
   0  0  1 -1 -1  0  0  0  0  0  0  0  0  0  0 -1  1 -1
   0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0
   0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0
   0  0  0  0  0  0  0  1  0 -1  0 -1  0  0  0  0  0  0
   0  0  0  0  0  0  0  0  1  0 -1  0 -1  0  0  0  0  0
   0  0  0 -1  0  1  0 -1  0  1  0  0  0  1  0  0  0  0
   0  0  0  0 -1  0  1  0 -1  0  1  0  0  0  1  0  0  0 ]
Bounding Constraints:
Simplex; tableau =
[  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
   0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
   0  0  1 -1 -1  0  0  0  0  0  0  0  0  0  0  1 -1 -1  0  0
   0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0 -1  0
   0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0 -1
   0  0  0  0  0  0  0  1  0 -1  0  1  0  0  0  0  0  0  0  0
   0  0  0  0  0  0  0  0  1  0 -1  0  1  0  0  0  0  0  0  0
   0  0  0 -1  0  1  0 -1  0  1  0  0  0 -1  0  0  0  0  0  0
   0  0  0  0 -1  0  1  0 -1  0  1  0  0  0 -1  0  0  0  0  0 ]
        Input:
Load:   %71 = load double, double addrspace(13)* %70, align 8, !dbg !68, !tbaa !61
ar.indexMatrix() =
[  1  0
   0  1 ]
ArrayReference %55 (dim = 2, num loops: 2, element size: 8):
A.numRow() = 2; A.numCol() = 2
Sizes: [ unknown, %20 ]
Subscripts: [ i_1 , i_0 ]
Schedule Omega: [ 0, 1, 0 ]
AffineLoopNest: alnb.getNumLoops() = 2
aln.getNumLoops() = 2
alnb.A =
[ -1  1  0 -1  0
  -1  0  1  0 -1
   0  0  0  1  0
   0  0  0  0  1 ]
Loop 1 lower bounds:
i_1 >= 0
Loop 1 upper bounds:
i_1 <= -1 + %23
i = 0; getNumSymbols() = 3
Loop 0 lower bounds:
i_0 >= 0
Loop 0 upper bounds:
i_0 <= -1 + %20

        Output:
Store:   store double %76, double addrspace(13)* %70, align 8, !dbg !73, !tbaa !61
ar.indexMatrix() =
[  1  0
   0  1 ]
ArrayReference %55 (dim = 2, num loops: 2, element size: 8):
A.numRow() = 2; A.numCol() = 2
Sizes: [ unknown, %20 ]
Subscripts: [ i_1 , i_0 ]
Schedule Omega: [ 0, 1, 2 ]
AffineLoopNest: alnb.getNumLoops() = 2
aln.getNumLoops() = 2
alnb.A =
[ -1  1  0 -1  0
  -1  0  1  0 -1
   0  0  0  1  0
   0  0  0  0  1 ]
Loop 1 lower bounds:
i_1 >= 0
Loop 1 upper bounds:
i_1 <= -1 + %23
i = 0; getNumSymbols() = 3
Loop 0 lower bounds:
i_0 >= 0
Loop 0 upper bounds:
i_0 <= -1 + %20

Schedule In:
nodeIndex = BitSet[1]; ref = ar.indexMatrix() =
[  1  0
   0  1 ]
ArrayReference %55 (dim = 2, num loops: 2, element size: 8):
A.numRow() = 2; A.numCol() = 2
Sizes: [ unknown, %20 ]
Subscripts: [ i_1 , i_0 ]
s.getPhi()
[  1  0
   0  1 ]
s.getFusionOmega() = [ 0, 0, 0 ]
s.getOffsetOmega() = [ 1, 0 ]

Schedule Out:
nodeIndex = BitSet[1, 2]; ref = ar.indexMatrix() =
[  1  0
   0  1 ]
ArrayReference %55 (dim = 2, num loops: 2, element size: 8):
A.numRow() = 2; A.numCol() = 2
Sizes: [ unknown, %20 ]
Subscripts: [ i_1 , i_0 ]
s.getPhi()
[  1  0
   0  1 ]
s.getFusionOmega() = [ 0, 0, 0 ]
s.getOffsetOmega() = [ 1, 0 ]

Schedule Out:
nodeIndex = BitSet[1, 2]; ref = ar.indexMatrix() =
[  1  0
   0  1 ]
ArrayReference %55 (dim = 2, num loops: 2, element size: 8):
A.numRow() = 2; A.numCol() = 2
Sizes: [ unknown, %20 ]
Subscripts: [ i_1 , i_0 ]
s.getPhi()
[  1  1  0
   0  0  1
   0  1  0 ]
s.getFusionOmega() = [ 0, 0, 0, 0 ]
s.getOffsetOmega() = [ 0, 0, 0 ]

        Edge = Dependence Poly y -> x:
-v_0 - v_1 + v_3 + v_4 <= 0 + %20
v_4 >= 0
-v_0 - v_1 + 2v_3 <= 0 + %23
v_0 >= 0
v_1 >= 0
-v_0 - v_1 + v_3 == 1
-v_2 + v_4 == 0

A =
[  0  0  1  1  1  0 -1 -1
   0  0  0  0  0  0  0  1
   0  1  0  1  1  0 -2  0
   0  0  0  1  0  0  0  0
   0  0  0  0  1  0  0  0 ]
E =
[  1  0  0  1  1  0 -1  0
   0  0  0  0  0  1  0 -1 ]
Schedule Constraints:
Simplex; tableau =
[  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
   0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
   0  0  1  0  0  0  0  0  1  0 -1  0  0  0  0  0  0 -1  1 -1
   0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0
   0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
   0  0  0  1  0  1  1  0  1  0 -1  0 -1  0  0  0  0  0  0  0
   0  0  0  1  0  1  0  1  1  0 -1  0  0 -1  0  0  0  0  0  0
   0  0  0  0  0  0  0  0  0  1  0 -1  0  0 -1  0  0  0  0  0
   0  0  0 -1  0 -2  0  0 -1  0  1  0  0  0  0  1  0  0  0  0
   0  0  0 -1  1  0  0  0  0 -1  0  1  0  0  0  0  1  0  0  0 ]
Bounding Constraints:
Simplex; tableau =
[  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
   0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
   0  0  1  0  0  0  0  0  1  0 -1  0  0  0  0  0  0  1 -1 -1  0  0
   0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0 -1  0
   0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 -1
   0  0  0  1  0  1  1  0  1  0 -1  0  1  0  0  0  0  0  0  0  0  0
   0  0  0  1  0  1  0  1  1  0 -1  0  0  1  0  0  0  0  0  0  0  0
   0  0  0  0  0  0  0  0  0  1  0 -1  0  0  1  0  0  0  0  0  0  0
   0  0  0 -1  0 -2  0  0 -1  0  1  0  0  0  0 -1  0  0  0  0  0  0
   0  0  0 -1  1  0  0  0  0 -1  0  1  0  0  0  0 -1  0  0  0  0  0 ]
        Input:
Store:   store double %64, double addrspace(13)* %65, align 8, !dbg !63, !tbaa !61
ar.indexMatrix() =
[  1  0
   0  1 ]
ArrayReference %55 (dim = 2, num loops: 2, element size: 8):
A.numRow() = 2; A.numCol() = 2
Sizes: [ unknown, %20 ]
Subscripts: [ i_1 , i_0 ]
Schedule Omega: [ 0, 0, 1 ]
AffineLoopNest: alnb.getNumLoops() = 2
aln.getNumLoops() = 2
alnb.A =
[ -1  1  0 -1  0
  -1  0  1  0 -1
   0  0  0  1  0
   0  0  0  0  1 ]
Loop 1 lower bounds:
i_1 >= 0
Loop 1 upper bounds:
i_1 <= -1 + %23
i = 0; getNumSymbols() = 3
Loop 0 lower bounds:
i_0 >= 0
Loop 0 upper bounds:
i_0 <= -1 + %20

        Output:
Load:   %83 = load double, double addrspace(13)* %82, align 8, !dbg !91, !tbaa !61
ar.indexMatrix() =
[  1  0
   1  0
   0  1 ]
ArrayReference %55 (dim = 2, num loops: 3, element size: 8):
A.numRow() = 3; A.numCol() = 2
Sizes: [ unknown, %20 ]
Subscripts: [ i_1  + i_2  + 1, i_0 ]
Schedule Omega: [ 0, 1, 3, 0 ]
AffineLoopNest: alnb.getNumLoops() = 3
aln.getNumLoops() = 3
alnb.A =
[ -2  1  0 -1 -1  0
   0  0  0  0  0  1
  -1  0  1  0  0 -1
   0  0  0  1  0  0
   0  0  0  0  1  0 ]
Loop 2 lower bounds:
i_2 >= 0
Loop 2 upper bounds:
i_2 <= -2 + %23 - i_1
i = 0; getNumSymbols() = 3
Loop 1 lower bounds:
i_1 >= 0
Loop 1 upper bounds:
i_1 <= -2 + %23
i = 1; getNumSymbols() = 3
Loop 0 lower bounds:
i_0 >= 0
Loop 0 upper bounds:
i_0 <= -1 + %20

Schedule In:
nodeIndex = BitSet[0]; ref = ar.indexMatrix() =
[  1  0
   0  1 ]
ArrayReference %55 (dim = 2, num loops: 2, element size: 8):
A.numRow() = 2; A.numCol() = 2
Sizes: [ unknown, %20 ]
Subscripts: [ i_1 , i_0 ]
s.getPhi()
[  1  0
   0  1 ]
s.getFusionOmega() = [ 0, 0, 0 ]
s.getOffsetOmega() = [ 1, 0 ]

Schedule Out:
nodeIndex = BitSet[2]; ref = ar.indexMatrix() =
[  1  0
   1  0
   0  1 ]
ArrayReference %55 (dim = 2, num loops: 3, element size: 8):
A.numRow() = 3; A.numCol() = 2
Sizes: [ unknown, %20 ]
Subscripts: [ i_1  + i_2  + 1, i_0 ]
s.getPhi()
[  1  1  0
   0  0  1
   0  1  0 ]
s.getFusionOmega() = [ 0, 0, 0, 0 ]
s.getOffsetOmega() = [ 0, 0, 0 ]

        Edge = Dependence Poly x -> y:
-v_0 - v_1 + v_3 + v_4 <= 0 + %20
v_4 >= 0
-v_0 - v_1 + 2v_3 <= 0 + %23
v_0 >= 0
v_1 >= 0
-v_0 - v_1 + v_3 == 1
-v_2 + v_4 == 0

A =
[  0  0  1  1  1  0 -1 -1
   0  0  0  0  0  0  0  1
   0  1  0  1  1  0 -2  0
   0  0  0  1  0  0  0  0
   0  0  0  0  1  0  0  0 ]
E =
[  1  0  0  1  1  0 -1  0
   0  0  0  0  0  1  0 -1 ]
Schedule Constraints:
Simplex; tableau =
[  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
   0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
   0  0  1  0  0  0  0  0  1  0 -1  0  0  0  0  0  0  1 -1 -1
   0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0
   0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
   0  0  0  1  0  1  1  0  1  0 -1  0  1  0  0  0  0  0  0  0
   0  0  0  1  0  1  0  1  1  0 -1  0  0  1  0  0  0  0  0  0
   0  0  0  0  0  0  0  0  0  1  0 -1  0  0  1  0  0  0  0  0
   0  0  0 -1  0 -2  0  0 -1  0  1  0  0  0  0 -1  0  0  0  0
   0  0  0 -1  1  0  0  0  0 -1  0  1  0  0  0  0 -1  0  0  0 ]
Bounding Constraints:
Simplex; tableau =
[  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
   0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
   0  0  1  0  0  0  0  0  1  0 -1  0  0  0  0  0  0 -1  1 -1  0  0
   0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0 -1  0
   0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 -1
   0  0  0  1  0  1  1  0  1  0 -1  0 -1  0  0  0  0  0  0  0  0  0
   0  0  0  1  0  1  0  1  1  0 -1  0  0 -1  0  0  0  0  0  0  0  0
   0  0  0  0  0  0  0  0  0  1  0 -1  0  0 -1  0  0  0  0  0  0  0
   0  0  0 -1  0 -2  0  0 -1  0  1  0  0  0  0  1  0  0  0  0  0  0
   0  0  0 -1  1  0  0  0  0 -1  0  1  0  0  0  0  1  0  0  0  0  0 ]
        Input:
Load:   %83 = load double, double addrspace(13)* %82, align 8, !dbg !91, !tbaa !61
ar.indexMatrix() =
[  1  0
   1  0
   0  1 ]
ArrayReference %55 (dim = 2, num loops: 3, element size: 8):
A.numRow() = 3; A.numCol() = 2
Sizes: [ unknown, %20 ]
Subscripts: [ i_1  + i_2  + 1, i_0 ]
Schedule Omega: [ 0, 1, 3, 0 ]
AffineLoopNest: alnb.getNumLoops() = 3
aln.getNumLoops() = 3
alnb.A =
[ -2  1  0 -1 -1  0
   0  0  0  0  0  1
  -1  0  1  0  0 -1
   0  0  0  1  0  0
   0  0  0  0  1  0 ]
Loop 2 lower bounds:
i_2 >= 0
Loop 2 upper bounds:
i_2 <= -2 + %23 - i_1
i = 0; getNumSymbols() = 3
Loop 1 lower bounds:
i_1 >= 0
Loop 1 upper bounds:
i_1 <= -2 + %23
i = 1; getNumSymbols() = 3
Loop 0 lower bounds:
i_0 >= 0
Loop 0 upper bounds:
i_0 <= -1 + %20

        Output:
Store:   store double %76, double addrspace(13)* %70, align 8, !dbg !73, !tbaa !61
ar.indexMatrix() =
[  1  0
   0  1 ]
ArrayReference %55 (dim = 2, num loops: 2, element size: 8):
A.numRow() = 2; A.numCol() = 2
Sizes: [ unknown, %20 ]
Subscripts: [ i_1 , i_0 ]
Schedule Omega: [ 0, 1, 2 ]
AffineLoopNest: alnb.getNumLoops() = 2
aln.getNumLoops() = 2
alnb.A =
[ -1  1  0 -1  0
  -1  0  1  0 -1
   0  0  0  1  0
   0  0  0  0  1 ]
Loop 1 lower bounds:
i_1 >= 0
Loop 1 upper bounds:
i_1 <= -1 + %23
i = 0; getNumSymbols() = 3
Loop 0 lower bounds:
i_0 >= 0
Loop 0 upper bounds:
i_0 <= -1 + %20

Schedule In:
nodeIndex = BitSet[2]; ref = ar.indexMatrix() =
[  1  0
   1  0
   0  1 ]
ArrayReference %55 (dim = 2, num loops: 3, element size: 8):
A.numRow() = 3; A.numCol() = 2
Sizes: [ unknown, %20 ]
Subscripts: [ i_1  + i_2  + 1, i_0 ]
s.getPhi()
[  1  1  0
   0  0  1
   0  1  0 ]
s.getFusionOmega() = [ 0, 0, 0, 0 ]
s.getOffsetOmega() = [ 0, 0, 0 ]

Schedule Out:
nodeIndex = BitSet[1, 2]; ref = ar.indexMatrix() =
[  1  0
   0  1 ]
ArrayReference %55 (dim = 2, num loops: 2, element size: 8):
A.numRow() = 2; A.numCol() = 2
Sizes: [ unknown, %20 ]
Subscripts: [ i_1 , i_0 ]
s.getPhi()
[  1  0
   0  1 ]
s.getFusionOmega() = [ 0, 0, 0 ]
s.getOffsetOmega() = [ 1, 0 ]

Schedule Out:
nodeIndex = BitSet[1, 2]; ref = ar.indexMatrix() =
[  1  0
   0  1 ]
ArrayReference %55 (dim = 2, num loops: 2, element size: 8):
A.numRow() = 2; A.numCol() = 2
Sizes: [ unknown, %20 ]
Subscripts: [ i_1 , i_0 ]
s.getPhi()
[  1  1  0
   0  0  1
   0  1  0 ]
s.getFusionOmega() = [ 0, 0, 0, 0 ]
s.getOffsetOmega() = [ 0, 0, 0 ]

        Edge = Dependence Poly y -> x:
-v_0 - v_1 + v_3 + v_4 <= 0 + %20
v_4 >= 0
-v_0 - v_1 + 2v_3 <= 0 + %23
v_0 >= 0
v_1 >= 0
-v_0 - v_1 + v_3 == 1
-v_2 + v_4 == 0

A =
[  0  0  1  1  1  0 -1 -1
   0  0  0  0  0  0  0  1
   0  1  0  1  1  0 -2  0
   0  0  0  1  0  0  0  0
   0  0  0  0  1  0  0  0 ]
E =
[  1  0  0  1  1  0 -1  0
   0  0  0  0  0  1  0 -1 ]
Schedule Constraints:
Simplex; tableau =
[  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
   0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
   0  0  1  0  0  0  0  0  1  0 -1  0  0  0  0  0  0 -1  1 -1
   0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0
   0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
   0  0  0  1  0  1  1  0  1  0 -1  0 -1  0  0  0  0  0  0  0
   0  0  0  1  0  1  0  1  1  0 -1  0  0 -1  0  0  0  0  0  0
   0  0  0  0  0  0  0  0  0  1  0 -1  0  0 -1  0  0  0  0  0
   0  0  0 -1  0 -2  0  0 -1  0  1  0  0  0  0  1  0  0  0  0
   0  0  0 -1  1  0  0  0  0 -1  0  1  0  0  0  0  1  0  0  0 ]
Bounding Constraints:
Simplex; tableau =
[  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
   0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
   0  0  1  0  0  0  0  0  1  0 -1  0  0  0  0  0  0  1 -1 -1  0  0
   0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0 -1  0
   0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 -1
   0  0  0  1  0  1  1  0  1  0 -1  0  1  0  0  0  0  0  0  0  0  0
   0  0  0  1  0  1  0  1  1  0 -1  0  0  1  0  0  0  0  0  0  0  0
   0  0  0  0  0  0  0  0  0  1  0 -1  0  0  1  0  0  0  0  0  0  0
   0  0  0 -1  0 -2  0  0 -1  0  1  0  0  0  0 -1  0  0  0  0  0  0
   0  0  0 -1  1  0  0  0  0 -1  0  1  0  0  0  0 -1  0  0  0  0  0 ]
        Input:
Store:   store double %64, double addrspace(13)* %65, align 8, !dbg !63, !tbaa !61
ar.indexMatrix() =
[  1  0
   0  1 ]
ArrayReference %55 (dim = 2, num loops: 2, element size: 8):
A.numRow() = 2; A.numCol() = 2
Sizes: [ unknown, %20 ]
Subscripts: [ i_1 , i_0 ]
Schedule Omega: [ 0, 0, 1 ]
AffineLoopNest: alnb.getNumLoops() = 2
aln.getNumLoops() = 2
alnb.A =
[ -1  1  0 -1  0
  -1  0  1  0 -1
   0  0  0  1  0
   0  0  0  0  1 ]
Loop 1 lower bounds:
i_1 >= 0
Loop 1 upper bounds:
i_1 <= -1 + %23
i = 0; getNumSymbols() = 3
Loop 0 lower bounds:
i_0 >= 0
Loop 0 upper bounds:
i_0 <= -1 + %20

        Output:
Store:   store double %89, double addrspace(13)* %82, align 8, !dbg !97, !tbaa !61
ar.indexMatrix() =
[  1  0
   1  0
   0  1 ]
ArrayReference %55 (dim = 2, num loops: 3, element size: 8):
A.numRow() = 3; A.numCol() = 2
Sizes: [ unknown, %20 ]
Subscripts: [ i_1  + i_2  + 1, i_0 ]
Schedule Omega: [ 0, 1, 3, 2 ]
AffineLoopNest: alnb.getNumLoops() = 3
aln.getNumLoops() = 3
alnb.A =
[ -2  1  0 -1 -1  0
   0  0  0  0  0  1
  -1  0  1  0  0 -1
   0  0  0  1  0  0
   0  0  0  0  1  0 ]
Loop 2 lower bounds:
i_2 >= 0
Loop 2 upper bounds:
i_2 <= -2 + %23 - i_1
i = 0; getNumSymbols() = 3
Loop 1 lower bounds:
i_1 >= 0
Loop 1 upper bounds:
i_1 <= -2 + %23
i = 1; getNumSymbols() = 3
Loop 0 lower bounds:
i_0 >= 0
Loop 0 upper bounds:
i_0 <= -1 + %20

Schedule In:
nodeIndex = BitSet[0]; ref = ar.indexMatrix() =
[  1  0
   0  1 ]
ArrayReference %55 (dim = 2, num loops: 2, element size: 8):
A.numRow() = 2; A.numCol() = 2
Sizes: [ unknown, %20 ]
Subscripts: [ i_1 , i_0 ]
s.getPhi()
[  1  0
   0  1 ]
s.getFusionOmega() = [ 0, 0, 0 ]
s.getOffsetOmega() = [ 1, 0 ]

Schedule Out:
nodeIndex = BitSet[2]; ref = ar.indexMatrix() =
[  1  0
   1  0
   0  1 ]
ArrayReference %55 (dim = 2, num loops: 3, element size: 8):
A.numRow() = 3; A.numCol() = 2
Sizes: [ unknown, %20 ]
Subscripts: [ i_1  + i_2  + 1, i_0 ]
s.getPhi()
[  1  1  0
   0  0  1
   0  1  0 ]
s.getFusionOmega() = [ 0, 0, 0, 0 ]
s.getOffsetOmega() = [ 0, 0, 0 ]

        Edge = Dependence Poly x -> y:
-v_0 - v_1 + v_3 + v_4 <= 0 + %20
v_4 >= 0
-v_0 - v_1 + 2v_3 <= 0 + %23
v_0 >= 0
v_1 >= 0
-v_0 - v_1 + v_3 == 1
-v_2 + v_4 == 0

A =
[  0  0  1  1  1  0 -1 -1
   0  0  0  0  0  0  0  1
   0  1  0  1  1  0 -2  0
   0  0  0  1  0  0  0  0
   0  0  0  0  1  0  0  0 ]
E =
[  1  0  0  1  1  0 -1  0
   0  0  0  0  0  1  0 -1 ]
Schedule Constraints:
Simplex; tableau =
[  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
   0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
   0  0  1  0  0  0  0  0  1  0 -1  0  0  0  0  0  0  1 -1 -1
   0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0
   0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
   0  0  0  1  0  1  1  0  1  0 -1  0  1  0  0  0  0  0  0  0
   0  0  0  1  0  1  0  1  1  0 -1  0  0  1  0  0  0  0  0  0
   0  0  0  0  0  0  0  0  0  1  0 -1  0  0  1  0  0  0  0  0
   0  0  0 -1  0 -2  0  0 -1  0  1  0  0  0  0 -1  0  0  0  0
   0  0  0 -1  1  0  0  0  0 -1  0  1  0  0  0  0 -1  0  0  0 ]
Bounding Constraints:
Simplex; tableau =
[  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
   0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
   0  0  1  0  0  0  0  0  1  0 -1  0  0  0  0  0  0 -1  1 -1  0  0
   0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0 -1  0
   0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 -1
   0  0  0  1  0  1  1  0  1  0 -1  0 -1  0  0  0  0  0  0  0  0  0
   0  0  0  1  0  1  0  1  1  0 -1  0  0 -1  0  0  0  0  0  0  0  0
   0  0  0  0  0  0  0  0  0  1  0 -1  0  0 -1  0  0  0  0  0  0  0
   0  0  0 -1  0 -2  0  0 -1  0  1  0  0  0  0  1  0  0  0  0  0  0
   0  0  0 -1  1  0  0  0  0 -1  0  1  0  0  0  0  1  0  0  0  0  0 ]
        Input:
Store:   store double %89, double addrspace(13)* %82, align 8, !dbg !97, !tbaa !61
ar.indexMatrix() =
[  1  0
   1  0
   0  1 ]
ArrayReference %55 (dim = 2, num loops: 3, element size: 8):
A.numRow() = 3; A.numCol() = 2
Sizes: [ unknown, %20 ]
Subscripts: [ i_1  + i_2  + 1, i_0 ]
Schedule Omega: [ 0, 1, 3, 2 ]
AffineLoopNest: alnb.getNumLoops() = 3
aln.getNumLoops() = 3
alnb.A =
[ -2  1  0 -1 -1  0
   0  0  0  0  0  1
  -1  0  1  0  0 -1
   0  0  0  1  0  0
   0  0  0  0  1  0 ]
Loop 2 lower bounds:
i_2 >= 0
Loop 2 upper bounds:
i_2 <= -2 + %23 - i_1
i = 0; getNumSymbols() = 3
Loop 1 lower bounds:
i_1 >= 0
Loop 1 upper bounds:
i_1 <= -2 + %23
i = 1; getNumSymbols() = 3
Loop 0 lower bounds:
i_0 >= 0
Loop 0 upper bounds:
i_0 <= -1 + %20

        Output:
Load:   %71 = load double, double addrspace(13)* %70, align 8, !dbg !68, !tbaa !61
ar.indexMatrix() =
[  1  0
   0  1 ]
ArrayReference %55 (dim = 2, num loops: 2, element size: 8):
A.numRow() = 2; A.numCol() = 2
Sizes: [ unknown, %20 ]
Subscripts: [ i_1 , i_0 ]
Schedule Omega: [ 0, 1, 0 ]
AffineLoopNest: alnb.getNumLoops() = 2
aln.getNumLoops() = 2
alnb.A =
[ -1  1  0 -1  0
  -1  0  1  0 -1
   0  0  0  1  0
   0  0  0  0  1 ]
Loop 1 lower bounds:
i_1 >= 0
Loop 1 upper bounds:
i_1 <= -1 + %23
i = 0; getNumSymbols() = 3
Loop 0 lower bounds:
i_0 >= 0
Loop 0 upper bounds:
i_0 <= -1 + %20

Schedule In:
nodeIndex = BitSet[2]; ref = ar.indexMatrix() =
[  1  0
   1  0
   0  1 ]
ArrayReference %55 (dim = 2, num loops: 3, element size: 8):
A.numRow() = 3; A.numCol() = 2
Sizes: [ unknown, %20 ]
Subscripts: [ i_1  + i_2  + 1, i_0 ]
s.getPhi()
[  1  1  0
   0  0  1
   0  1  0 ]
s.getFusionOmega() = [ 0, 0, 0, 0 ]
s.getOffsetOmega() = [ 0, 0, 0 ]

Schedule Out:
nodeIndex = BitSet[1]; ref = ar.indexMatrix() =
[  1  0
   0  1 ]
ArrayReference %55 (dim = 2, num loops: 2, element size: 8):
A.numRow() = 2; A.numCol() = 2
Sizes: [ unknown, %20 ]
Subscripts: [ i_1 , i_0 ]
s.getPhi()
[  1  0
   0  1 ]
s.getFusionOmega() = [ 0, 0, 0 ]
s.getOffsetOmega() = [ 1, 0 ]

        Edge = Dependence Poly x -> y:
-v_0 - v_1 + v_3 + v_4 <= 0 + %20
v_4 >= 0
-v_0 - v_1 + 2v_3 <= 0 + %23
v_0 >= 0
v_1 >= 0
-v_0 - v_1 + v_3 == 1
-v_2 + v_4 == 0

A =
[  0  0  1  1  1  0 -1 -1
   0  0  0  0  0  0  0  1
   0  1  0  1  1  0 -2  0
   0  0  0  1  0  0  0  0
   0  0  0  0  1  0  0  0 ]
E =
[  1  0  0  1  1  0 -1  0
   0  0  0  0  0  1  0 -1 ]
Schedule Constraints:
Simplex; tableau =
[  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
   0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
   0  0  1  0  0  0  0  0  1  0 -1  0  0  0  0  0  0  1 -1 -1
   0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0
   0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
   0  0  0  1  0  1  1  0  1  0 -1  0  1  0  0  0  0  0  0  0
   0  0  0  1  0  1  0  1  1  0 -1  0  0  1  0  0  0  0  0  0
   0  0  0  0  0  0  0  0  0  1  0 -1  0  0  1  0  0  0  0  0
   0  0  0 -1  0 -2  0  0 -1  0  1  0  0  0  0 -1  0  0  0  0
   0  0  0 -1  1  0  0  0  0 -1  0  1  0  0  0  0 -1  0  0  0 ]
Bounding Constraints:
Simplex; tableau =
[  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
   0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
   0  0  1  0  0  0  0  0  1  0 -1  0  0  0  0  0  0 -1  1 -1  0  0
   0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0 -1  0
   0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 -1
   0  0  0  1  0  1  1  0  1  0 -1  0 -1  0  0  0  0  0  0  0  0  0
   0  0  0  1  0  1  0  1  1  0 -1  0  0 -1  0  0  0  0  0  0  0  0
   0  0  0  0  0  0  0  0  0  1  0 -1  0  0 -1  0  0  0  0  0  0  0
   0  0  0 -1  0 -2  0  0 -1  0  1  0  0  0  0  1  0  0  0  0  0  0
   0  0  0 -1  1  0  0  0  0 -1  0  1  0  0  0  0  1  0  0  0  0  0 ]
        Input:
Store:   store double %89, double addrspace(13)* %82, align 8, !dbg !97, !tbaa !61
ar.indexMatrix() =
[  1  0
   1  0
   0  1 ]
ArrayReference %55 (dim = 2, num loops: 3, element size: 8):
A.numRow() = 3; A.numCol() = 2
Sizes: [ unknown, %20 ]
Subscripts: [ i_1  + i_2  + 1, i_0 ]
Schedule Omega: [ 0, 1, 3, 2 ]
AffineLoopNest: alnb.getNumLoops() = 3
aln.getNumLoops() = 3
alnb.A =
[ -2  1  0 -1 -1  0
   0  0  0  0  0  1
  -1  0  1  0  0 -1
   0  0  0  1  0  0
   0  0  0  0  1  0 ]
Loop 2 lower bounds:
i_2 >= 0
Loop 2 upper bounds:
i_2 <= -2 + %23 - i_1
i = 0; getNumSymbols() = 3
Loop 1 lower bounds:
i_1 >= 0
Loop 1 upper bounds:
i_1 <= -2 + %23
i = 1; getNumSymbols() = 3
Loop 0 lower bounds:
i_0 >= 0
Loop 0 upper bounds:
i_0 <= -1 + %20

        Output:
Store:   store double %76, double addrspace(13)* %70, align 8, !dbg !73, !tbaa !61
ar.indexMatrix() =
[  1  0
   0  1 ]
ArrayReference %55 (dim = 2, num loops: 2, element size: 8):
A.numRow() = 2; A.numCol() = 2
Sizes: [ unknown, %20 ]
Subscripts: [ i_1 , i_0 ]
Schedule Omega: [ 0, 1, 2 ]
AffineLoopNest: alnb.getNumLoops() = 2
aln.getNumLoops() = 2
alnb.A =
[ -1  1  0 -1  0
  -1  0  1  0 -1
   0  0  0  1  0
   0  0  0  0  1 ]
Loop 1 lower bounds:
i_1 >= 0
Loop 1 upper bounds:
i_1 <= -1 + %23
i = 0; getNumSymbols() = 3
Loop 0 lower bounds:
i_0 >= 0
Loop 0 upper bounds:
i_0 <= -1 + %20

Schedule In:
nodeIndex = BitSet[2]; ref = ar.indexMatrix() =
[  1  0
   1  0
   0  1 ]
ArrayReference %55 (dim = 2, num loops: 3, element size: 8):
A.numRow() = 3; A.numCol() = 2
Sizes: [ unknown, %20 ]
Subscripts: [ i_1  + i_2  + 1, i_0 ]
s.getPhi()
[  1  1  0
   0  0  1
   0  1  0 ]
s.getFusionOmega() = [ 0, 0, 0, 0 ]
s.getOffsetOmega() = [ 0, 0, 0 ]

Schedule Out:
nodeIndex = BitSet[1, 2]; ref = ar.indexMatrix() =
[  1  0
   0  1 ]
ArrayReference %55 (dim = 2, num loops: 2, element size: 8):
A.numRow() = 2; A.numCol() = 2
Sizes: [ unknown, %20 ]
Subscripts: [ i_1 , i_0 ]
s.getPhi()
[  1  0
   0  1 ]
s.getFusionOmega() = [ 0, 0, 0 ]
s.getOffsetOmega() = [ 1, 0 ]

Schedule Out:
nodeIndex = BitSet[1, 2]; ref = ar.indexMatrix() =
[  1  0
   0  1 ]
ArrayReference %55 (dim = 2, num loops: 2, element size: 8):
A.numRow() = 2; A.numCol() = 2
Sizes: [ unknown, %20 ]
Subscripts: [ i_1 , i_0 ]
s.getPhi()
[  1  1  0
   0  0  1
   0  1  0 ]
s.getFusionOmega() = [ 0, 0, 0, 0 ]
s.getOffsetOmega() = [ 0, 0, 0 ]

        Edge = Dependence Poly y -> x:
v_5 <= -1 + %20
v_3 >= 0
v_4 >= 0
-v_1 - v_3 + v_4 <= 0
v_1 >= 0
v_3 + v_4 <= -2 + %23
v_5 >= 0
-v_0 - v_1 + v_3 + v_4 == 0
-2v_1 + 2v_4 == 0
-v_2 + v_5 == 0

A =
[ -1  0  1  0  0  0  0  0 -1
   0  0  0  0  0  0  1  0  0
   0  0  0  0  0  0  0  1  0
   0  0  0  0  1  0  1 -1  0
   0  0  0  0  1  0  0  0  0
  -2  1  0  0  0  0 -1 -1  0
   0  0  0  0  0  0  0  0  1 ]
E =
[  0  0  0  1  1  0 -1 -1  0
   0  0  0  0  2  0  0 -2  0
   0  0  0  0  0  1  0  0 -1 ]
Schedule Constraints:
Simplex; tableau =
[  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
   0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
   0  0  1 -1  0  0  0  0 -2  0  0  0  0  0  0  0  0  0  0  0  0  0 -1  1 -1
   0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
   0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
   0  0  0  0  0  0  0  0  0  0  1  0  0 -1  0  0 -1  0  0  0  0  0  0  0  0
   0  0  0  0  0  0  1  1  0  0  1  2  0 -1 -2  0  0 -1  0  0  0  0  0  0  0
   0  0  0  0  0  0  0  0  0  0  0  0  1  0  0 -1  0  0 -1  0  0  0  0  0  0
   0  0  0  0  1  0  1  0 -1  0 -1  0  0  1  0  0  0  0  0  1  0  0  0  0  0
   0  0  0  0  0  1 -1  0 -1  0 -1 -2  0  1  2  0  0  0  0  0  1  0  0  0  0
   0  0  0 -1  0  0  0  0  0  1  0  0 -1  0  0  1  0  0  0  0  0  1  0  0  0 ]
Bounding Constraints:
Simplex; tableau =
[  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
   0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
   0  0  1 -1  0  0  0  0 -2  0  0  0  0  0  0  0  0  0  0  0  0  0  1 -1 -1  0  0
   0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 -1  0
   0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 -1
   0  0  0  0  0  0  0  0  0  0  1  0  0 -1  0  0  1  0  0  0  0  0  0  0  0  0  0
   0  0  0  0  0  0  1  1  0  0  1  2  0 -1 -2  0  0  1  0  0  0  0  0  0  0  0  0
   0  0  0  0  0  0  0  0  0  0  0  0  1  0  0 -1  0  0  1  0  0  0  0  0  0  0  0
   0  0  0  0  1  0  1  0 -1  0 -1  0  0  1  0  0  0  0  0 -1  0  0  0  0  0  0  0
   0  0  0  0  0  1 -1  0 -1  0 -1 -2  0  1  2  0  0  0  0  0 -1  0  0  0  0  0  0
   0  0  0 -1  0  0  0  0  0  1  0  0 -1  0  0  1  0  0  0  0  0 -1  0  0  0  0  0 ]
        Input:
Load:   %83 = load double, double addrspace(13)* %82, align 8, !dbg !91, !tbaa !61
ar.indexMatrix() =
[  1  0
   1  0
   0  1 ]
ArrayReference %55 (dim = 2, num loops: 3, element size: 8):
A.numRow() = 3; A.numCol() = 2
Sizes: [ unknown, %20 ]
Subscripts: [ i_1  + i_2  + 1, i_0 ]
Schedule Omega: [ 0, 1, 3, 0 ]
AffineLoopNest: alnb.getNumLoops() = 3
aln.getNumLoops() = 3
alnb.A =
[ -2  1  0 -1 -1  0
   0  0  0  0  0  1
  -1  0  1  0  0 -1
   0  0  0  1  0  0
   0  0  0  0  1  0 ]
Loop 2 lower bounds:
i_2 >= 0
Loop 2 upper bounds:
i_2 <= -2 + %23 - i_1
i = 0; getNumSymbols() = 3
Loop 1 lower bounds:
i_1 >= 0
Loop 1 upper bounds:
i_1 <= -2 + %23
i = 1; getNumSymbols() = 3
Loop 0 lower bounds:
i_0 >= 0
Loop 0 upper bounds:
i_0 <= -1 + %20

        Output:
Store:   store double %89, double addrspace(13)* %82, align 8, !dbg !97, !tbaa !61
ar.indexMatrix() =
[  1  0
   1  0
   0  1 ]
ArrayReference %55 (dim = 2, num loops: 3, element size: 8):
A.numRow() = 3; A.numCol() = 2
Sizes: [ unknown, %20 ]
Subscripts: [ i_1  + i_2  + 1, i_0 ]
Schedule Omega: [ 0, 1, 3, 2 ]
AffineLoopNest: alnb.getNumLoops() = 3
aln.getNumLoops() = 3
alnb.A =
[ -2  1  0 -1 -1  0
   0  0  0  0  0  1
  -1  0  1  0  0 -1
   0  0  0  1  0  0
   0  0  0  0  1  0 ]
Loop 2 lower bounds:
i_2 >= 0
Loop 2 upper bounds:
i_2 <= -2 + %23 - i_1
i = 0; getNumSymbols() = 3
Loop 1 lower bounds:
i_1 >= 0
Loop 1 upper bounds:
i_1 <= -2 + %23
i = 1; getNumSymbols() = 3
Loop 0 lower bounds:
i_0 >= 0
Loop 0 upper bounds:
i_0 <= -1 + %20

Schedule In:
nodeIndex = BitSet[2]; ref = ar.indexMatrix() =
[  1  0
   1  0
   0  1 ]
ArrayReference %55 (dim = 2, num loops: 3, element size: 8):
A.numRow() = 3; A.numCol() = 2
Sizes: [ unknown, %20 ]
Subscripts: [ i_1  + i_2  + 1, i_0 ]
s.getPhi()
[  1  1  0
   0  0  1
   0  1  0 ]
s.getFusionOmega() = [ 0, 0, 0, 0 ]
s.getOffsetOmega() = [ 0, 0, 0 ]

Schedule Out:
nodeIndex = BitSet[2]; ref = ar.indexMatrix() =
[  1  0
   1  0
   0  1 ]
ArrayReference %55 (dim = 2, num loops: 3, element size: 8):
A.numRow() = 3; A.numCol() = 2
Sizes: [ unknown, %20 ]
Subscripts: [ i_1  + i_2  + 1, i_0 ]
s.getPhi()
[  1  1  0
   0  0  1
   0  1  0 ]
s.getFusionOmega() = [ 0, 0, 0, 0 ]
s.getOffsetOmega() = [ 0, 0, 0 ]

        Edge = Dependence Poly x -> y:
v_5 <= -1 + %20
v_3 >= 0
v_4 >= 0
-v_1 - v_3 + v_4 <= 2
v_1 >= 0
v_3 + v_4 <= -2 + %23
v_5 >= 0
-v_0 - v_1 + v_3 + v_4 == 0
-2v_1 + 2v_4 == 2
-v_2 + v_5 == 0

A =
[ -1  0  1  0  0  0  0  0 -1
   0  0  0  0  0  0  1  0  0
   0  0  0  0  0  0  0  1  0
   2  0  0  0  1  0  1 -1  0
   0  0  0  0  1  0  0  0  0
  -2  1  0  0  0  0 -1 -1  0
   0  0  0  0  0  0  0  0  1 ]
E =
[  0  0  0  1  1  0 -1 -1  0
   2  0  0  0  2  0  0 -2  0
   0  0  0  0  0  1  0  0 -1 ]
Schedule Constraints:
Simplex; tableau =
[  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
   0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
   0  0  1 -1  0  0  2  0 -2  0  0  2  0  0 -2  0  0  0  0  0  0  0  1 -1 -1
   0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
   0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
   0  0  0  0  0  0  0  0  0  0  1  0  0 -1  0  0  1  0  0  0  0  0  0  0  0
   0  0  0  0  0  0  1  1  0  0  1  2  0 -1 -2  0  0  1  0  0  0  0  0  0  0
   0  0  0  0  0  0  0  0  0  0  0  0  1  0  0 -1  0  0  1  0  0  0  0  0  0
   0  0  0  0  1  0  1  0 -1  0 -1  0  0  1  0  0  0  0  0 -1  0  0  0  0  0
   0  0  0  0  0  1 -1  0 -1  0 -1 -2  0  1  2  0  0  0  0  0 -1  0  0  0  0
   0  0  0 -1  0  0  0  0  0  1  0  0 -1  0  0  1  0  0  0  0  0 -1  0  0  0 ]
Bounding Constraints:
Simplex; tableau =
[  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
   0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
   0  0  1 -1  0  0  2  0 -2  0  0  2  0  0 -2  0  0  0  0  0  0  0 -1  1 -1  0  0
   0  0  0  0  0  0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 -1  0
   0  0  0  1  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0 -1
   0  0  0  0  0  0  0  0  0  0  1  0  0 -1  0  0 -1  0  0  0  0  0  0  0  0  0  0
   0  0  0  0  0  0  1  1  0  0  1  2  0 -1 -2  0  0 -1  0  0  0  0  0  0  0  0  0
   0  0  0  0  0  0  0  0  0  0  0  0  1  0  0 -1  0  0 -1  0  0  0  0  0  0  0  0
   0  0  0  0  1  0  1  0 -1  0 -1  0  0  1  0  0  0  0  0  1  0  0  0  0  0  0  0
   0  0  0  0  0  1 -1  0 -1  0 -1 -2  0  1  2  0  0  0  0  0  1  0  0  0  0  0  0
   0  0  0 -1  0  0  0  0  0  1  0  0 -1  0  0  1  0  0  0  0  0  1  0  0  0  0  0 ]
        Input:
Store:   store double %89, double addrspace(13)* %82, align 8, !dbg !97, !tbaa !61
ar.indexMatrix() =
[  1  0
   1  0
   0  1 ]
ArrayReference %55 (dim = 2, num loops: 3, element size: 8):
A.numRow() = 3; A.numCol() = 2
Sizes: [ unknown, %20 ]
Subscripts: [ i_1  + i_2  + 1, i_0 ]
Schedule Omega: [ 0, 1, 3, 2 ]
AffineLoopNest: alnb.getNumLoops() = 3
aln.getNumLoops() = 3
alnb.A =
[ -2  1  0 -1 -1  0
   0  0  0  0  0  1
  -1  0  1  0  0 -1
   0  0  0  1  0  0
   0  0  0  0  1  0 ]
Loop 2 lower bounds:
i_2 >= 0
Loop 2 upper bounds:
i_2 <= -2 + %23 - i_1
i = 0; getNumSymbols() = 3
Loop 1 lower bounds:
i_1 >= 0
Loop 1 upper bounds:
i_1 <= -2 + %23
i = 1; getNumSymbols() = 3
Loop 0 lower bounds:
i_0 >= 0
Loop 0 upper bounds:
i_0 <= -1 + %20

        Output:
Load:   %83 = load double, double addrspace(13)* %82, align 8, !dbg !91, !tbaa !61
ar.indexMatrix() =
[  1  0
   1  0
   0  1 ]
ArrayReference %55 (dim = 2, num loops: 3, element size: 8):
A.numRow() = 3; A.numCol() = 2
Sizes: [ unknown, %20 ]
Subscripts: [ i_1  + i_2  + 1, i_0 ]
Schedule Omega: [ 0, 1, 3, 0 ]
AffineLoopNest: alnb.getNumLoops() = 3
aln.getNumLoops() = 3
alnb.A =
[ -2  1  0 -1 -1  0
   0  0  0  0  0  1
  -1  0  1  0  0 -1
   0  0  0  1  0  0
   0  0  0  0  1  0 ]
Loop 2 lower bounds:
i_2 >= 0
Loop 2 upper bounds:
i_2 <= -2 + %23 - i_1
i = 0; getNumSymbols() = 3
Loop 1 lower bounds:
i_1 >= 0
Loop 1 upper bounds:
i_1 <= -2 + %23
i = 1; getNumSymbols() = 3
Loop 0 lower bounds:
i_0 >= 0
Loop 0 upper bounds:
i_0 <= -1 + %20

Schedule In:
nodeIndex = BitSet[2]; ref = ar.indexMatrix() =
[  1  0
   1  0
   0  1 ]
ArrayReference %55 (dim = 2, num loops: 3, element size: 8):
A.numRow() = 3; A.numCol() = 2
Sizes: [ unknown, %20 ]
Subscripts: [ i_1  + i_2  + 1, i_0 ]
s.getPhi()
[  1  1  0
   0  0  1
   0  1  0 ]
s.getFusionOmega() = [ 0, 0, 0, 0 ]
s.getOffsetOmega() = [ 0, 0, 0 ]

Schedule Out:
nodeIndex = BitSet[2]; ref = ar.indexMatrix() =
[  1  0
   1  0
   0  1 ]
ArrayReference %55 (dim = 2, num loops: 3, element size: 8):
A.numRow() = 3; A.numCol() = 2
Sizes: [ unknown, %20 ]
Subscripts: [ i_1  + i_2  + 1, i_0 ]
s.getPhi()
[  1  1  0
   0  0  1
   0  1  0 ]
s.getFusionOmega() = [ 0, 0, 0, 0 ]
s.getOffsetOmega() = [ 0, 0, 0 ]


LoopBlock schedule (#mem accesses = 8):

Ref = ar.indexMatrix() =
[  1  0
   0  1 ]
ArrayReference %52 (dim = 2, num loops: 2, element size: 8):
A.numRow() = 2; A.numCol() = 2
Sizes: [ unknown, %20 ]
Subscripts: [ i_1 , i_0 ]
nodeIndex = 0
s.getPhi()
[  1  0
   0  1 ]
s.getFusionOmega() = [ 0, 0, 0 ]
s.getOffsetOmega() = [ 1, 0 ]
Ref = ar.indexMatrix() =
[  1  0
   0  1 ]
ArrayReference %55 (dim = 2, num loops: 2, element size: 8):
A.numRow() = 2; A.numCol() = 2
Sizes: [ unknown, %20 ]
Subscripts: [ i_1 , i_0 ]
nodeIndex = 0
s.getPhi()
[  1  0
   0  1 ]
s.getFusionOmega() = [ 0, 0, 0 ]
s.getOffsetOmega() = [ 1, 0 ]
Ref = ar.indexMatrix() =
[  1  0
   0  1 ]
ArrayReference %55 (dim = 2, num loops: 2, element size: 8):
A.numRow() = 2; A.numCol() = 2
Sizes: [ unknown, %20 ]
Subscripts: [ i_1 , i_0 ]
nodeIndex = 1
s.getPhi()
[  1  0
   0  1 ]
s.getFusionOmega() = [ 0, 0, 0 ]
s.getOffsetOmega() = [ 1, 0 ]
Ref = ar.indexMatrix() =
[  1  1
   0  0 ]
ArrayReference %58 (dim = 2, num loops: 2, element size: 8):
A.numRow() = 2; A.numCol() = 2
Sizes: [ unknown, %23 ]
Subscripts: [ i_1 , i_1 ]
nodeIndex = 1
s.getPhi()
[  1  0
   0  1 ]
s.getFusionOmega() = [ 0, 0, 0 ]
s.getOffsetOmega() = [ 1, 0 ]
Ref = ar.indexMatrix() =
[  1  0
   0  1 ]
ArrayReference %55 (dim = 2, num loops: 2, element size: 8):
A.numRow() = 2; A.numCol() = 2
Sizes: [ unknown, %20 ]
Subscripts: [ i_1 , i_0 ]
nodeIndex = 1
s.getPhi()
[  1  0
   0  1 ]
s.getFusionOmega() = [ 0, 0, 0 ]
s.getOffsetOmega() = [ 1, 0 ]

nodeIndex = 2
s.getPhi()
[  1  1  0
   0  0  1
   0  1  0 ]
s.getFusionOmega() = [ 0, 0, 0, 0 ]
s.getOffsetOmega() = [ 0, 0, 0 ]
Ref = ar.indexMatrix() =
[  1  0
   1  0
   0  1 ]
ArrayReference %55 (dim = 2, num loops: 3, element size: 8):
A.numRow() = 3; A.numCol() = 2
Sizes: [ unknown, %20 ]
Subscripts: [ i_1  + i_2  + 1, i_0 ]
nodeIndex = 2
s.getPhi()
[  1  1  0
   0  0  1
   0  1  0 ]
s.getFusionOmega() = [ 0, 0, 0, 0 ]
s.getOffsetOmega() = [ 0, 0, 0 ]
Ref = ar.indexMatrix() =
[  1  0
   1  1
   0  0 ]
ArrayReference %58 (dim = 2, num loops: 3, element size: 8):
A.numRow() = 3; A.numCol() = 2
Sizes: [ unknown, %23 ]
Subscripts: [ i_1  + i_2  + 1, i_1 ]
nodeIndex = 2
s.getPhi()
[  1  1  0
   0  0  1
   0  1  0 ]
s.getFusionOmega() = [ 0, 0, 0, 0 ]
s.getOffsetOmega() = [ 0, 0, 0 ]
Ref = ar.indexMatrix() =
[  1  0
   1  0
   0  1 ]
ArrayReference %55 (dim = 2, num loops: 3, element size: 8):
A.numRow() = 3; A.numCol() = 2
Sizes: [ unknown, %20 ]
Subscripts: [ i_1  + i_2  + 1, i_0 ]
nodeIndex = 2
s.getPhi()
[  1  1  0
   0  0  1
   0  1  0 ]
s.getFusionOmega() = [ 0, 0, 0, 0 ]
s.getOffsetOmega() = [ 0, 0, 0 ]
```
The value `-9223372036854775808` is used as a sentinel value indicating that there are no dependencies between separate graph vertices at that level (but there may be within the graph); the idea is that we'll solve for values here while determining vectorization and unrolling; the solution must satisfy dependencies within the edge, as well as be linearly independent of the set of solutions of loops exterior to it. We don't have any examples of that here.

A few more notes on interpretation:
Each vertex of the graph corresponds to a single store, which is the last memory access listed; all other memory accesses are loads. So, for example, in `v_2`:
```
v_2:
mem =
  store double %76, double addrspace(13)* %70, align 8, !dbg !73, !tbaa !61
  %83 = load double, double addrspace(13)* %82, align 8, !dbg !91, !tbaa !61
  %87 = load double, double addrspace(13)* %86, align 8, !dbg !91, !tbaa !61
  store double %89, double addrspace(13)* %82, align 8, !dbg !97, !tbaa !61
inNeighbors = v_0, v_1, v_2,
outNeighbors = v_1, v_2,
```
This corresponds to `store double %89, double addrspace(13)* %82`.
The other store listed above is actually a reload of that particular store.

Additionally, all our analysis of loop bounds and array indexing are built on SCEV (Scalar Evolution). Things like calculating trip counts, or using [AddRecExpr](https://llvm.org/doxygen/classllvm_1_1SCEVAddRecExpr.html) for indexing make it natural to treat all loops as starting at `0`.
Array indexing of course is also naturally treated as pointer-offsets, and thus starts at `0` as well.

For this reason, the Julia loop we showed earlier using 1-based indexing is understood internally as something closer to:
```julia
using OffsetArrays
function triangular_solve0!(_A,_B,_U)
    A = OffsetArray(_A, OffsetArrays.Origin(0))
	B = OffsetArray(_B, OffsetArrays.Origin(0))
    U = OffsetArray(_U, OffsetArrays.Origin(0))
    M,N = size(A)
    @assert M == size(B,1)
    @assert N == size(B,2)
    @assert N == LinearAlgebra.checksquare(U)
    @inbounds for m = 0:M-1
        for n = 0:N-1
            A[m,n] = B[m,n]
        end
        for n = 0:N-1
            A[m,n] /= U[n,n]
            for k = 0:N-1-n
                A[m,k+n+1] -= A[m,n]*U[n,k+n+1]
            end
        end
    end
end
```
The original loop nest structure is generally represented as being inner<->outer.
This is to ease parsing: we want an affine loop nest, so we start parsing loops from the inner-most loop outward, until we either reach toplevel, or hit something non-affine.

The rows of a schedule correspond to the schedules of loops from outer-most to inner-most, however. So, the interpretation of a schedule matrix is as follows:
```
s.getPhi()
[  1  1  0   # Outer most loop: k + n
   0  0  1   # Middle loop: m
   0  1  0 ] # Inner most loop: n
```
Note that this recovers our non-zero start, as we recreate a loop induction variable that equals the sum `k + n` as we originally had. In this case, that is desirable as we can now hoist the loads and stores `A[m,k+n+1]` out of the inner-most loop. Through unrolling both the `k+n` and `m` loop (and vectorizing the `m` loop), we can perform register tiling, where we use a large number of vector registers as the accumulators for the inner most loop, updating them on each iteration of the inner loop, maximizing the number of arithmetic operations relative to loads and stores. 

