import "cgm"

module S = mk_cgm_real f32

-- ==
-- entry: cgm_sparse_test
-- input { [ [0f32, 0f32, 4f32], [1f32, 0f32, 1f32], [0f32, 1f32, 1f32], [1f32, 1f32, 3f32] ] [ 1f32,  2f32 ] [ 2f32,  1f32 ] 1i64 }
-- output { [ 0.23564959f32, 0.3383686f32 ] }
entry cgm_sparse_test [n][nnz] (A: [nnz][3]f32) (b: [n]f32) (x: [n]f32) (r: i64): [n]f32 =
	let c = #iterations r
	let mat = map (\i -> (i64.f32 A[i, 0], i64.f32 A[i, 1], A[i, 2])) (iota nnz)
	let mat = S.sparse n mat
	in S.cgm mat x b c
