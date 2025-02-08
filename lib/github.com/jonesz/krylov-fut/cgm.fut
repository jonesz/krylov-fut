-- | CGM implementation.
import "../../diku-dk/linalg/linalg"

module type cgm = {
  -- | The scalar type.
  type t
  type m [n]
  val cgm [n] : (A: [n]m[n]) -> (b: [n]t) -> (x_0: [n]t) -> (residual: t) -> (max_iter: i64) -> [n]t
}

module mk_cgm (R: real) : cgm with t = R.t with m [n] = [n]R.t = {
  type t = R.t
  type m [n] = [n]t
  module L = mk_linalg R

  def cgm A b x_0 residual max_iter =
    let r_0 = L.matvecmul_row A x_0 |> map2 (R.-) b
    let (x_star, _, _, _) =
      loop (x_k, r_k, p_k, i) = (x_0, r_0, r_0, 0i64)
      while ((R.>) (L.vecnorm r_k) residual) && (i < max_iter) do
        let a_k = L.matvecmul_row A p_k |> L.dotprod p_k |> (R./) (L.dotprod r_k r_k)
        let x_l = L.vecscale a_k p_k |> map2 (R.+) x_k
        let r_l = L.matvecmul_row A p_k |> L.vecscale a_k |> map2 (R.-) r_k
        let B_k = L.dotprod r_k r_k |> (R./) (L.dotprod r_l r_l)
        let p_l = L.vecscale B_k p_k |> map2 (R.+) r_l
        in (x_l, r_l, p_l, i + 1i64)
    in x_star
}
