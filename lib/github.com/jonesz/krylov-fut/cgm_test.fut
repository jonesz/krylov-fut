import "cgm"

module S = mk_cgm f32

-- ==
-- entry: cgm_dense_real
-- input { [[4f32, 1f32], [1f32, 3f32]] [1f32, 2f32] [21f32, 1f32] }
-- output { [0.0909f32, 0.6364f32] }
entry cgm_dense_real [n] (A: [n][n]f32) (b: [n]f32) (x: [n]f32): [n]f32 =
	S.cgm A b x 0.001f32 100
