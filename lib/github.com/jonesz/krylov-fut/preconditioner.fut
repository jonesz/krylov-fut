module type preconditioner = {
  -- | The underlying scalar type.
  type t

  -- | The representation of the input matrix `A`.
  type~ A [n]

  -- | The output representation of the preconditioner.
  type~ M [n]

  -- | Given some matrix `A`, return the preconditioner `M` inverse.
  val preconditioner [n] : A [n] -> M [n]

  -- | Solve Mx = b for x.
  val solve [n] : M [n] -> [n]t -> [n]t
}

-- | A Jacobi preconditioner, which is simply the diagonal of A; suitable
-- for use in diagonally dominant matrices.
module mk_jacobi
  (R: real)
  (I: {
    type t
    type~ mat [n]
    val diag [n] : mat [n] -> [n]t
  }
  with t = R.t)
  : preconditioner
    with A [n] = I.mat [n]
    with t = R.t
    with M [n] = [n]R.t = {
  type t = R.t
  type~ A [n] = I.mat [n]
  type~ M [n] = [n]t

  -- The Jacobi preconditioner is the diagonal.
  def preconditioner A = I.diag A

  -- Solve a linear system Jx = b.
  def solve J b = map (R.recip) J |> map2 (R.*) b
}
