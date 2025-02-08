import "cgm"
import "../sparse/symmetric"

module S_d = mk_cgm f32
module Sym = mk_symmetric {open f32 def conj x = x}
module S_s = mk_cgm_sparse f32 Sym

-- ==
-- entry: cgm_dense_real
-- input { [[4f32, 1f32], [1f32, 3f32]] [1f32, 2f32] [21f32, 1f32] }
-- output { [0.0909f32, 0.6364f32] }
entry cgm_dense_real [n] (A: [n][n]f32) (b: [n]f32) (x: [n]f32): [n]f32 =
	S_d.cgm A b x 0.001f32 100

-- ==
-- entry: cgm_sparse_real
-- input { [4f32, 1f32, 3f32] [1f32, 2f32] [21f32, 1f32] }
-- output { [0.0909f32, 0.6364f32] }
entry cgm_sparse_real [n] A (b: [n]f32) (x: [n]f32): [n]f32 =
	let A = Sym.symmetric A
	in S_s.cgm A b x 0.001f32 100
