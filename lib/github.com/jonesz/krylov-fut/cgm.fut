import "../../diku-dk/linalg/linalg"
import "../../diku-dk/sparse/compressed"

module type cgm = {
	type t
	type~ mat[n]

	type end_criteria = #iterations i64 | #residual t 

	val sparse [nnz] : (n: i64) -> [nnz](i64, i64, t) -> mat[n]
	val cgm [n] : mat[n] -> [n]t -> [n]t -> end_criteria -> [n]t
}

module mk_cgm_real(T: field): cgm with t = T.t = {
	type t = T.t
	type end_criteria = #iterations i64 | #residual t 

	open (mk_compressed T)
	type~ mat[n] = sr.mat[n][n]
	
	def sparse [nnz] (n: i64) (coord: [nnz](i64, i64, t)): mat[n] =
		sr.sparse n n coord

	def cgm [n] (A: mat[n]) (x: [n]t) (b: [n]t) (c: end_criteria): [n]t =

		-- b - Ax.
		let explicit_residual (A: mat[n]) (x: [n]t) (b: [n]t): [n]t =		
			sr.smvm A x |> map2 (T.-) b

		-- (rk1^T * rk1) / (rk^T * rk)
		let beta rk rk1 = 
			let dp x y = map2 (T.*) x y |> reduce (T.+) (T.i64 0)
			in (dp rk1 rk1) T./ (dp rk rk)

		let proj_krylov (rk) (pk) (xk) =
			-- A * pk
			let Apk = sr.smvm A pk
		
			-- alpha_k := (r^T * r) / p^T * A * p
			let alpha_k =
				let n = map2 (T.*) rk rk |> reduce (T.+) (T.i64 0)
				let d = Apk -- A is real symmetric: xA = (Ax)^T.
					|> map2 (T.*) pk 
					|> reduce (T.+) (T.i64 0)
				in n T./ d

			-- xk1 = xk + alpha * pk.
			let xk1 = map (T.* alpha_k) pk |> map2 (T.+) xk
			-- rk1 = rk - alpha_k * A * pk.
			let rk1 = map (T.* alpha_k) Apk |> map2 (T.-) rk

			in (xk1, rk1)

		-- Each loop is required to have the same (_, _, _, _) parameter
		-- structure; hopefully (T.i64 0) will compile down to a no-op.

		-- Compute r iterations.
		let loop_iter r0 p0 x0 r = 
			loop (rk, pk, xk, _rk_mag) = (r0, p0, x0, (T.i64 0)) for _i < r do
				let (xk1, rk1) = proj_krylov rk pk xk
				let bk = beta rk rk1
				-- pk1 = rk1 + bk * pk
				let pk1 = map (T.* bk) pk |> map2 (T.+) rk1
				in (rk1, pk1, xk1, (T.i64 0))

		
		-- Loop until a residual < r is found.
		let loop_residual r0 p0 x0 r =
			loop (rk, pk, xk, rk_mag) = (r0, p0, x0, reduce (T.+) (T.i64 0) r0) while (T.<) r rk_mag do
				let (xk1, rk1) = proj_krylov rk pk xk
				let bk = beta rk rk1
				-- pk1 = rk1 + bk * pk
				let pk1 = map (T.* bk) pk |> map2 (T.+) rk1
				let rk_mag1 = reduce (T.+) (T.i64 0) rk1
				in (rk1, pk1, xk1, rk_mag1)
					
		let r0 = explicit_residual A x b
		let p0 = r0
		let x0 = x

		let (_, _, xn, _) = match c
			case #iterations r -> loop_iter r0 p0 x0 r
			case #residual r -> loop_residual r0 p0 x0 r
		in xn
}