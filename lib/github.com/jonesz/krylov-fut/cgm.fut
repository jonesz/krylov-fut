-- 2022, Ethan Jones <etn.jones@gmail.com>.

import "symmetric"

-- https://futhark-lang.org/examples/matrix-multiplication.html
def matmul [n][m][p] 'a
           (add: a -> a -> a) (mul: a -> a -> a) (zero: a)
           (A: [n][m]a) (B: [m][p]a) : [n][p]a =
  map (\A_row ->
         map (\B_col ->
                reduce add zero (map2 mul A_row B_col))
             (transpose B))
      A

entry matmul_f32 = matmul (+) (*) 0f32

-- Conjugate Gradient Method
-- https://en.wikipedia.org/wiki/Conjugate_gradient_method
-- Returns the value 'x' after 'r' iterations.
entry cgm [n] (A: [n][n]f32) (b: [n][1]f32) (x: [n][1]f32) (r: i64): [n][1]f32 =

	-- Compute the residual of b - Ax.
	let residual [n] (A: [n][n]f32) (b: [n][1]f32) (x: [n][1]f32): [n][1]f32 =
		matmul_f32 A x |> map2 (map2 (-)) b

	-- Calculation of beta: (rk1^T * rk1) / (rk^T * rk)
	let beta [n] (rk: [n][1]f32) (rk1: [n][1]f32): f32 =
    	let n = matmul_f32 (transpose rk1) rk1
    	let d = matmul_f32 (transpose rk) rk
    	in n[0, 0] / d[0, 0] -- Both of these matrices should be [1][1].

	let inner_loop (rk) (pk) (xk) =
		-- ak := (r^T * r) / p^T * A * p
  		let alpha_k = 
    		let n = matmul_f32 (transpose rk) rk
    		let d = matmul_f32 (matmul_f32 (transpose pk) A) pk
    		in n[0, 0] / d[0, 0] -- Both of these matrices should be [1][1].
		
		-- See note in "Explicit residual calculation". 
		-- TODO: Implement the explicit/implicit residual calculations with
		-- a match on a sum type.
  		let xk = map2 (map2 (+)) xk <| map (map (* alpha_k)) pk
		let rk = map2 (map2 (-)) rk <| matmul_f32 (map (map (* alpha_k)) A) pk
		in (xk, rk)

  let outer_loop r0 p0 x0 = 
	loop (rk, pk, xk) = (r0, p0, x0) for _i < r do
  		let (xk1, rk1) = inner_loop rk pk xk
		let bk = beta rk rk1
		let pk1 = map (map (* bk)) pk |> map2 (map2 (+)) rk1 
		in (rk1, pk1, xk1)

  -- r0 := b - Ax
  let r0 = residual A b x
  let p0 = r0
	
  let (_, _, xn) = outer_loop r0 p0 x
  in xn

module type cgm_impl = {
	--| Scalar type.
	type t 

	type~ mat[n]

	module sym : {
		include symmetric_matrix with t = t with mat [n] = mat[n]
	
		val cgm [n] : mat[n] -> [n]t -> [n]t -> i64 -> [n]t
	}
}

module mk_cgm_impl (T: numeric) = {
	type t = T.t

	module sym = {
		open (mk_symmetric_matrix_def T)

	-- Compute the residual of b - Ax.
		local def residual [n] (A: mat[n]) (b: [n]t) (x: [n]t): [n]t =
			mul_vec A x |> map2 (T.-) b

		local def beta [n] (rk: [n]t) (rk1: [n]t): t =
			let n = reduce_comm (T.+) (T.i64 0) (map2 (T.*) rk1 rk1)
			let d = reduce_comm (T.+) (T.i64 0) (map2 (T.*) rk rk)
			in n T./ d

		def cgm [n] (A: mat[n]) (b: [n]t) (x: [n]t) (r: i64): [n]t = 
			let inner_loop (rk) (pk) (xk) =
			-- ak := (r^T * r) / p^T * A * p
				let alpha_k =
					let n = map (\x -> (T.*) x x) rk |> reduce (T.+) (T.i64 0)
					let d = vec_mul pk A |> map2 (T.*) pk |> reduce (T.+) (T.i64 0)
					in n T./ d

				let xk = map2 (T.+) xk <| map (T.* alpha_k) pk
				let rk = map2 (T.-) rk <| mul_vec (scale A alpha_k) pk

				in (xk, rk)

			let outer_loop_a r0 p0 x0 =
				loop (rk, pk, xk) = (r0, p0, x0) for _i < r do
					let (xk1, rk1) = inner_loop rk pk xk
					let bk = beta rk rk1
					let pk1 = map (T.* bk) pk |> map2 (T.+) rk1
					in (rk1, pk1, xk1)

			let r0 = residual A b x
			let p0 = r0

			let (_, _, xn) = outer_loop_a r0 p0 x
			in xn
	}

	type~ mat[n] = sym.mat[n]
}
