module type preconditioner = {
	type t 
	type~ m [n]
	type~ P [n] 

	val M [n] : m[n] -> P[n]
	val solve [n] : P[n] -> (b: [n]t) -> [n]t
}

module mk_jacobi(R: real)
(I: {
	type t
	type~ m [n]
	val diag [n] : m[n] -> [n]t
} with t = R.t): preconditioner with t = R.t with m [n] = I.m[n] with P [n] = [n]R.t = {
	type t = R.t
	type~ m [n] = I.m [n]
	type~ P [n] = [n]t

	-- The Jacobi preconditioner is simply the diagonal of `A`; the inverse of a diagonal
	-- is the reciprocal across the diagonal.
	def M A = I.diag A |> map (R.recip)
	def solve J b = map2 (R.*) J b
}
