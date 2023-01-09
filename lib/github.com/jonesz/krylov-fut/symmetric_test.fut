import "symmetric"
module S = mk_symmetric_matrix_def f32


-- ==
-- entry: ssmm_test
-- input { [[ 4f32, 1f32 ], [ 1f32, 3f32 ]] }
-- output { [[ 17f32, 7f32 ], [7f32, 10f32 ]] }
entry ssmm_test [n] (A: [n][n]f32): [n][n]f32 =
	let mat = S.sym A
	in S.ssmm mat mat

-- entry: sdmm_dense_test
-- input { [[ 4f32, 1f32 ], [ 1f32, 3f32 ]] }
-- output { [[ 17f32, 7f32 ], [7f32, 10f32 ]] }
entry sdmm_test [n] (A: [n][n]f32): [n][n]f32 =
	let mat = S.sym A
	in S.sdmm mat A

-- entry: dsmm_test
-- input { [[ 4f32, 1f32 ], [ 1f32, 3f32 ]] }
-- output { [[ 17f32, 7f32 ], [ 7f32, 10f32 ]] }
entry dsmm_test [n] (A: [n][n]f32): [n][n]f32 =
	let mat = S.sym A
	in S.dsmm A mat