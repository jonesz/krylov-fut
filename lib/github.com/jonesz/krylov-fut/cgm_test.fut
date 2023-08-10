import "cgm"

module S = mk_cgm { open f32 def conj a : f32 = a }

-- ==
-- entry: cgm_dense_real
-- input { [ [4f32, 1f32], [1f32, 3f32] ] [ [1f32], [2f32] ] [ [2f32], [1f32] ] 1i64 }
-- output { [ [0.23564959f32], [0.3383686f32] ] }
entry cgm_dense_real [n] (A: [n][n]f32) (b: [n][1]f32) (x: [n][1]f32) (i: i64): [n][1]f32 =
	S.cgm (#iter 1i64) A x b
