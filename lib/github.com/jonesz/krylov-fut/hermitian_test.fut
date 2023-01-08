import "hermitian"
import "../../diku-dk/complex/complex"

module c32 = mk_complex f32
module H = mk_hermitian_matrix_def c32

-- ==
-- entry: idx_test
-- input { [[[ 1f32, 2f32 ], [ 2f32, 1f32 ]], [[ 3f32, 3f32 ], [ 2f32, -1f32 ]]] }
-- output { [ 2f32, -1f32 ] }
entry idx_test [n] (A: [n][n][2]f32): [2]f32 =
	let f = map (\a -> c32.mk a[0] a[1])
	let idx =  map f A |> H.herm |> H.idx (1, 0)
	in [idx.0, idx.1]
