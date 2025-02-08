import "cgm"
import "preconditioner"
import "../sparse/symmetric"

module S = mk_symmetric {open f32 def conj x = x}
module D = mk_cgm_dense f32
module E = mk_cgm_symmetric f32 S

module J = mk_jacobi f32 {
  open S
  type~ m [n] = mat [n]

  -- TODO: These `m`/`mat` types should be better named.
  def diag [n] A = map (\i -> idx (i, i) A) (iota n)
}

module EJ = mk_cgm_symmetric_preconditioner f32 S J

-- ==
-- entry: cgm_dense
-- input { [[4f32, 1f32], [1f32, 3f32]] [1f32, 2f32] [21f32, 1f32] }
-- output { [0.0909f32, 0.6364f32] }
entry cgm_dense [n] (A: [n][n]f32) (b: [n]f32) (x: [n]f32) : [n]f32 =
  D.cgm A b x 0.001f32 100

-- ==
-- entry: cgm_symmetric_real
-- input { [4f32, 1f32, 3f32] [1f32, 2f32] [21f32, 1f32] }
-- output { [0.0909f32, 0.6364f32] }
entry cgm_symmetric_real [n] A (b: [n]f32) (x: [n]f32) : [n]f32 =
  let A = S.symmetric A
  in E.cgm A b x 0.001f32 100

-- ==
-- entry: cgm_symmetric_jacobi_real
-- input { [4f32, 1f32, 3f32] [1f32, 2f32] [21f32, 1f32] }
-- output { [0.0909f32, 0.6364f32] }
entry cgm_symmetric_jacobi_real [n] A (b: [n]f32) (x: [n]f32) : [n]f32 =
  let A = S.symmetric A
  in EJ.cgm A b x 0.001f32 1000
