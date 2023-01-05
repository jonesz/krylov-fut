import "symmetric"

module S = mk_symmetric_matrix_def f32

-- ==
-- entry: mul_vec_test
-- input { [ 1f32, 2f32 ] [[ 4f32, 1f32 ], [ 1f32, 3f32 ]] }
-- output { [6f32, 7f32 ] }
entry mul_vec_test [n] (x: [n]f32) (A: [n][n]f32): [n]f32 =
	let mat = S.sym A
	in S.mul_vec mat x

-- ==
-- entry: smm_test
-- input { [[ 4f32, 1f32 ], [ 1f32, 3f32 ]] }
-- output { [[ 17f32, 7f32 ], [7f32, 10f32 ]] }
entry smm_test [n] (A: [n][n]f32): [n][n]f32 =
	let mat = S.sym A
	in S.smm mat mat

-- entry: smm_dense_test
-- input { [[ 4f32, 1f32 ], [ 1f32, 3f32 ]] }
-- output { [[ 17f32, 7f32 ], [7f32, 10f32 ]] }
entry smm_dense_test [n] (A: [n][n]f32): [n][n]f32 =
	let mat = S.sym A
	in S.smm_dense mat A

-- entry: dense_smm_test
-- input { [[ 4f32, 1f32 ], [ 1f32, 3f32 ]] }
-- output { [[ 17f32, 7f32 ], [ 7f32, 10f32 ]] }
entry dense_smm_test [n] (A: [n][n]f32): [n][n]f32 =
	let mat = S.sym A
	in S.dense_smm A mat