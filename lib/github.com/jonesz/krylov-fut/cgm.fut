-- 2022, Ethan Jones <etn.jones@gmail.com>.
import "symmetric"
import "hermitian"
import "../../diku-dk/linalg/linalg"
import "../../diku-dk/complex/complex"

-- https://futhark-lang.org/examples/matrix-multiplication.html
local def matmul [n][m][p] 'a
           (add: a -> a -> a) (mul: a -> a -> a) (zero: a)
           (A: [n][m]a) (B: [m][p]a) : [n][p]a =
  map (\A_row ->
         map (\B_col ->
                reduce add zero (map2 mul A_row B_col))
             (transpose B))
      A

module type cgm = {
	type t

	module dense : {
		type~ mat[n] = [n][n]t
		val cgm [n] : mat[n] -> [n][1]t -> [n][1]t -> i64 -> [n][1]t
	}

	module sparse : {
		include symmetric_matrix with t = t
		val cgm [n] : mat[n] -> [n][1]t -> [n][1]t -> i64 -> [n][1]t
	}
}

module mk_cgm_real(T: field): cgm with t = T.t = {
	type t = T.t

	-- dense * dense matrix multiply.
	let ddmm = matmul (T.+) (T.*) (T.i64 0)

	let beta [n] (rk: [n][1]t) (rk1: [n][1]t): t =
    	let n = ddmm (transpose rk1) rk1
    	let d = ddmm (transpose rk) rk
    	in n[0, 0] T./ d[0, 0]

	module dense = {
		type~ mat[n] = [n][n]t

		def residual [n] (A: mat[n]) (b: [n][1]t) (x: [n][1]t): [n][1]t =
			ddmm A x |> map2 (map2 (T.-)) b

		def cgm [n] (A: mat[n]) (b: [n][1]t) (x: [n][1]t) (r: i64): [n][1]t = 
			let inner_loop (rk) (pk) (xk) =
				-- ak := (r^T * r) / p^T * A * p
  				let alpha_k = 
    				let n = ddmm (transpose rk) rk
    				let d = ddmm (ddmm (transpose pk) A) pk
    				in n[0, 0] T./ d[0, 0]
		
				-- TODO: According to https://www.cs.cmu.edu/~quake-papers/painless-conjugate-gradient.pdf
				-- we want to outright calculate the residual every so-odd iterations; see page 14, Eq. 10,
				-- Eq. 13. 
  				let xk = map2 (map2 (T.+)) xk <| map (map (T.* alpha_k)) pk
				let rk = map2 (map2 (T.-)) rk <| ddmm (map (map (T.* alpha_k)) A) pk
				in (xk, rk)

  			let outer_loop r0 p0 x0 = 
				loop (rk, pk, xk) = (r0, p0, x0) for _i < r do
  					let (xk1, rk1) = inner_loop rk pk xk
					let bk = beta rk rk1
					let pk1 = map (map (T.* bk)) pk |> map2 (map2 (T.+)) rk1 
					in (rk1, pk1, xk1)

  			-- r0 := b - Ax
  			let r0 = residual A b x
  			let p0 = r0
	
  			let (_, _, xn) = outer_loop r0 p0 x
  			in xn
	}

	module sparse = {
		open (mk_symmetric_matrix_def T)

		def residual [n] (A: mat[n]) (b: [n][1]t) (x: [n][1]t): [n][1]t =
			sdmm A x |> map2 (map2 (T.-)) b

		def cgm [n] (A: mat[n]) (b: [n][1]t) (x: [n][1]t) (r: i64): [n][1]t = 
			let inner_loop (rk) (pk) (xk) =
			-- ak := (r^T * r) / p^T * A * p
				let alpha_k =
					let n = ddmm (transpose rk) rk
					let d = ddmm (dsmm (transpose pk) A) pk
    				in n[0, 0] T./ d[0, 0]

				-- TODO: According to https://www.cs.cmu.edu/~quake-papers/painless-conjugate-gradient.pdf
				-- we want to outright calculate the residual every so-odd iterations; see page 14, Eq. 10,
				-- Eq. 13. 
  				let xk = map2 (map2 (T.+)) xk <| map (map (T.* alpha_k)) pk
				let rk = map2 (map2 (T.-)) rk <| sdmm (scale alpha_k A) pk

				in (xk, rk)

			let outer_loop_a r0 p0 x0 =
				loop (rk, pk, xk) = (r0, p0, x0) for _i < r do
					let (xk1, rk1) = inner_loop rk pk xk
					let bk = beta rk rk1
					let pk1 = map (map (T.* bk)) pk |> map2 (map2 (T.+)) rk1 
					in (rk1, pk1, xk1)

			let r0 = residual A b x
			let p0 = r0

			let (_, _, xn) = outer_loop_a r0 p0 x
			in xn
	}
}

module mk_cgm_complex(T: complex): cgm with t = T.t = {
	type t = T.t

	-- dense * dense matrix multiply.
	let ddmm = matmul (T.+) (T.*) (T.i64 0)

	module dense = {
		type~ mat[n] = [n][n]t

		let conj_transpose A = map (map T.conj) A |> transpose

		def beta [n] (rk: [n][1]t) (rk1: [n][1]t): t =
    		let n = ddmm (conj_transpose rk1) rk1
    		let d = ddmm (conj_transpose rk) rk
    		in n[0, 0] T./ d[0, 0]

		def residual [n] (A: mat[n]) (b: [n][1]t) (x: [n][1]t): [n][1]t =
			ddmm A x |> map2 (map2 (T.-)) b

		def cgm [n] (A: mat[n]) (b: [n][1]t) (x: [n][1]t) (r: i64): [n][1]t = 
			let inner_loop (rk) (pk) (xk) =
				-- ak := (r^T * r) / p^T * A * p
  				let alpha_k = 
    				let n = ddmm (conj_transpose rk) rk
    				let d = ddmm (ddmm (conj_transpose pk) A) pk
    				in n[0, 0] T./ d[0, 0]
		
				-- TODO: According to https://www.cs.cmu.edu/~quake-papers/painless-conjugate-gradient.pdf
				-- we want to outright calculate the residual every so-odd iterations; see page 14, Eq. 10,
				-- Eq. 13. 
  				let xk = map2 (map2 (T.+)) xk <| map (map (T.* alpha_k)) pk
				let rk = map2 (map2 (T.-)) rk <| ddmm (map (map (T.* alpha_k)) A) pk
				in (xk, rk)

  			let outer_loop r0 p0 x0 = 
				loop (rk, pk, xk) = (r0, p0, x0) for _i < r do
  					let (xk1, rk1) = inner_loop rk pk xk
					let bk = beta rk rk1
					let pk1 = map (map (T.* bk)) pk |> map2 (map2 (T.+)) rk1 
					in (rk1, pk1, xk1)

  			-- r0 := b - Ax
  			let r0 = residual A b x
  			let p0 = r0
	
  			let (_, _, xn) = outer_loop r0 p0 x
  			in xn
	}
}
