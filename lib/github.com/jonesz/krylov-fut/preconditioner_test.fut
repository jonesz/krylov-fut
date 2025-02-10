import "preconditioner"
import "../sparse/symmetric"

module S = mk_symmetric {open f32 def conj x = x}

module J = mk_jacobi f32 {
  open S
  type~ A [n] = mat [n]
  def diag [n] A = map (\i -> idx (i, i) A) (iota n)
}

-- ==
-- entry: symmetric_jacobi_preconditioner
-- input { [4f32, 1f32, 10f32] }
-- output { [4f32, 10f32] }
entry symmetric_jacobi_preconditioner A : [2]f32 =
  S.symmetric A |> J.preconditioner

-- ==
-- entry: symmetric_jacobi_solve
-- input { [4f32, 1f32, 10f32] [20f32, 30f32] }
-- output { [5f32, 3f32] }
entry symmetric_jacobi_solve A b: [2]f32 =
  S.symmetric A |> J.preconditioner |> flip J.solve b
