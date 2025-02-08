import "preconditioner"
import "../sparse/symmetric"

module S = mk_symmetric {open f32 def conj x = x}
module J = mk_jacobi f32 { open S
			type~ m [n] = mat [n] -- TODO: These `m`/`mat` types should be better named.
			def diag [n] A = map (\i -> idx (i, i) A) (iota n)
		   }

-- ==
-- entry: jacobi_diag_symmetric
-- input { [4f32, 1f32, 10f32] }
-- output { [0.25f32, 0.1f32] }
entry jacobi_diag_symmetric A : [2]f32 =
	let A = S.symmetric A
	in J.M A
