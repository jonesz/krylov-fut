import "../../diku-dk/linalg/linalg"
import "../../jonesz/sparse/symmetric"
import "preconditioner"

-- | The basic conjugate gradient method.
module type cgm = {
  -- | The scalar type.
  type t
  type~ m [n]
  val cgm [n] : (A: m [n]) -> (b: [n]t) -> (x_0: [n]t) -> (residual: t) -> (max_iter: i64) -> [n]t
}

local module mk_cgm_inner
  (R: real)
  (I: {
    type t
    type~ m [n]
    val matvecmul_row [n] : m [n] -> [n]t -> [n]t
  }
  with t = R.t)
  : cgm with t = I.t with m [n] = I.m [n] = {
  type t = R.t
  type~ m [n] = I.m [n]
  module L = mk_linalg R

  def cgm [n] (A: m [n]) b x_0 residual max_iter =
    let r_0 = -- b - A * x_0
      I.matvecmul_row A x_0 |> map2 (R.-) b

    -- For the below, `*_l` are the values for the next iteration after `x_k`
    let (x_star, _, _, _) =
      loop (x_k, r_k, p_k, k) = (x_0, r_0, r_0, 0i64)
      while (L.vecnorm r_k |> (R.<) residual) && k < max_iter do

        -- And so begins the CGM iteration...
        -- TODO: A lot of this code is shared with a non-preconditioner CGM; is there
        -- way to refactor this?
        let a_k = -- a_k = (r_k^T * r_k) / (p_k^T * A * p_k)
          I.matvecmul_row A p_k |> L.dotprod p_k |> (R./) (L.dotprod r_k r_k)
        let x_l = -- x_l = x_k + a_k * p_k
          L.vecscale a_k p_k |> map2 (R.+) x_k
        let r_l = -- r_l = = r_k - a_k * A * p_k
          I.matvecmul_row A p_k |> L.vecscale a_k |> map2 (R.-) r_k
        let B_k = -- B_k = (r_l^T * r_l) / (r_k^T * r_k)
          L.dotprod r_k r_k |> (R./) (L.dotprod r_l r_l)
        let p_l = -- p_l = r_l + B_k * p_k
          L.vecscale B_k p_k |> map2 (R.+) r_l

        in (x_l, r_l, p_l, k + 1i64)
    in x_star
}

local module mk_cgm_inner_preconditioner
  (R: real)
  (I: {
    type t
    type~ m [n]
    val matvecmul_row [n] : m [n] -> [n]t -> [n]t
  }
  with t = R.t)
  (P: preconditioner with t = R.t with A [n] = I.m [n])
  : cgm with t = I.t with m [n] = I.m [n] = {
  type t = R.t
  type~ m [n] = I.m [n]
  module L = mk_linalg R

  def cgm [n] (A: m [n]) b x_0 residual max_iter =
    let r_0 = I.matvecmul_row A x_0 |> map2 (R.-) b -- r_0 = b - A * x_0
    let z_0 = P.preconditioner A |> flip P.solve b  -- z_0 = M_inv * r_0

    -- For the below, `*_l` are the values for the next iteration after `x_k`
    let (x_star, _, _, _, _) =
      loop (x_k, r_k, p_k, z_k, k) = (x_0, r_0, z_0, z_0, 0i64)
      while (L.vecnorm r_k |> (R.<) residual) && k < max_iter do

        -- And so begins the CGM iteration...
        -- TODO: A lot of this code is shared with a non-preconditioner CGM; is there
        -- way to refactor this?
        let a_k = -- a_k = (r_k^T * z_k)/(p_k^T * A * p_k)
          I.matvecmul_row A p_k |> L.dotprod p_k |> (R./) (L.dotprod r_k z_k)
        let x_l = -- x_l = x_k + a_k * p_k
          L.vecscale a_k p_k |> map2 (R.+) x_k
        let r_l = -- r_l = r_k - a_k * A * p_k
          I.matvecmul_row A p_k |> L.vecscale a_k |> map2 (R.-) r_k
        let z_l = -- z_l = M_inv r_l
          P.preconditioner A |> flip P.solve r_l
        let B_k = -- (r_l^T * z_l)/(r_k^T * z_k)
          L.dotprod r_k z_k |> (R./) (L.dotprod r_l z_l)
        let p_l = -- p_l = r_l + B_k * p_k
          L.vecscale B_k p_k |> map2 (R.+) z_l

        in (x_l, r_l, p_l, z_l, k + 1i64)
    in x_star
}

module mk_cgm_dense (R: real)
  : cgm with t = R.t
        with m [n] = [n][n]R.t = mk_cgm_inner R {
  module L = mk_linalg R
  type t = R.t
  type m [n] = [n][n]t
  def matvecmul_row = L.matvecmul_row
}

module mk_cgm_symmetric (R: real)
                        (S: symmetric_matrix with t = R.t)
  : cgm with t = R.t
        with m [n] = S.mat [n] = mk_cgm_inner R {
  type t = R.t
  type~ m [n] = S.mat [n]
  def matvecmul_row = S.smvm
}

module mk_cgm_symmetric_preconditioner (R: real)
                                       (S: symmetric_matrix with t = R.t)
                                       (P: preconditioner with t = R.t with A [n] = S.mat [n])
  : cgm with t = R.t
        with m [n] = S.mat [n] = mk_cgm_inner_preconditioner R {
  type t = R.t
  type~ m [n] = S.mat [n]
  def matvecmul_row = S.smvm
} P
