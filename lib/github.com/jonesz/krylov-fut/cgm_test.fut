import "cgm"

module S = mk_cgm f32

-- ==
-- entry: cgm_dense_real
-- input { [[2f32, 3f32], [4f32, 5f32]] [1f32, 1f32] [01f32, 0f32] }
-- output { [-1f32, 1f32] }
entry cgm_dense_real [n] (A: [n][n]f32) (b: [n]f32) (x: [n]f32): [n]f32 =
	S.cgm A b x 0.001f32 100
