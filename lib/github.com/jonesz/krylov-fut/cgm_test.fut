import "cgm"

module S = mk_cgm_impl f32

-- ==
-- entry: cgm_sym_test
-- input { [[ 4f32, 1f32 ], [ 1f32, 3f32 ]] [[ 1f32 ], [ 2f32 ]] [[ 2f32 ], [ 1f32 ]] 1i64 }
-- output { [[ 0.23564959f32 ], [ 0.3383686f32 ]] }
entry cgm_sym_test [n] (A: [n][n]f32) (b: [n][1]f32) (x: [n][1]f32) (r: i64): [n][1]f32 =
	let mat = S.sym.sym A
	in S.sym.cgm mat b x r

-- ==
-- entry: cgm_dense_test
-- input { [[ 4f32, 1f32 ], [ 1f32, 3f32 ]] [[ 1f32 ], [ 2f32 ]] [[ 2f32 ], [ 1f32 ]] 1i64 }
-- output { [[ 0.23564959f32 ], [ 0.3383686f32 ]] }
entry cgm_dense_test [n] (A: [n][n]f32) (b: [n][1]f32) (x: [n][1]f32) (r: i64): [n][1]f32 =
	S.dense.cgm A b x r
