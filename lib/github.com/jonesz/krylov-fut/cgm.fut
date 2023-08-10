-- | CGM implementation.
import "../../diku-dk/linalg/linalg"

-- https://futhark-lang.org/examples/matrix-multiplication.html
local def matmul [n][m][p] 'a
           (add: a -> a -> a) (mul: a -> a -> a) (zero: a)
           (A: [n][m]a) (B: [m][p]a) : [n][p]a =
	map (\A_row ->
         map (\B_col ->
                reduce add zero (map2 mul A_row B_col))
             (transpose B))
		A

-- Entrywise matrix operations (add, sub).
local def mat_entrywise [n][m] 'a
			(op: a -> a -> a) (x: [n][m]a) (y: [n][m]a) : [n][m]a =
	map2 (map2 (op)) x y

local def mat_map [n][m] 'a 'b (f: a -> b) (x: [n][m]a) : [n][m]b =
	map (map (f)) x

module type cgm = {
	-- | The scalar type.
	type t
	type criteria = #iter i64 | #residual t
	val cgm [n] : criteria -> [n][n]t -> [n][1]t -> [n][1]t -> [n][1]t
}

local module type field_plus_conj_transpose = {
	include field
	val conj : t -> t
	val < : t -> t -> bool -- TODO: Fix some of the linalg field stuff.
}

module mk_cgm(T: field_plus_conj_transpose): cgm with t = T.t = {
	type t = T.t
	type criteria = #iter i64 | #residual t

	let matmul_t = matmul (T.+) (T.*) (T.i64 0)
	let matsub_t = mat_entrywise (T.-)
	let matadd_t = mat_entrywise (T.+)
	let matscale_t a = mat_map (T.* a)

	let conj_transpose [n] (x: [n][1]t): [1][n]t =
		transpose x |> mat_map (T.conj)

	-- b - Ax
	let explicit_residual [n] (A: [n][n]t) (x: [n][1]t) (b: [n][1]t) = 
		matmul_t A x |> matsub_t b
				
	def cgm [n] (c: criteria) (A: [n][n]t) (x: [n][1]t) (b: [n][1]t): [n][1]t =

		-- Single iteration of the CGM method.
		let inner r_curr p_curr x_curr = 

			-- (r^T * r) / (p^T A p)
			let alpha [n] (A: [n][n]t) (r: [n][1]t) (p: [n][1]t): t =
				let top = matmul_t (conj_transpose r) r |> flatten |> head
				let bot = matmul_t (matmul_t (conj_transpose p) A) p |> flatten |> head
				in (T./) top bot

			-- (r1^T * r1) / (r0^T / r0)
			let beta [n] (r0: [n][1]t) (r1: [n][1]t): t =
				let top = matmul_t (conj_transpose r1) r1 |> flatten |> head
				let bot = matmul_t (conj_transpose r0) r0 |> flatten |> head
				in (T./) top bot

			let a = alpha A r_curr p_curr
			let x_next = matscale_t a p_curr |> matadd_t x_curr              -- xk + ak * pk.
			let r_next = matsub_t r_curr <| matmul_t (matscale_t a A) p_curr -- rk - ak * A * pk.

			let b = beta r_curr r_next
			let p_next = matadd_t r_next <| matscale_t b p_curr              -- rk1 + bk * pk.
			in (r_next, p_next, x_next)

		let loop_residual r p x (residual: t) =
			let sum_error e = flatten e |> reduce (T.+) (T.i64 0)
			let loop_condition r k = ((T.<) residual <| sum_error r) && k < 1000 -- TODO: What should the `MAX_ITER` value be?
			in loop (r_curr, p_curr, x_curr, k) = (r, p, x, 0i64) while loop_condition r_curr k do
				let (r_next, p_next, x_next) = inner r_curr p_curr x_curr
				in (r_next, p_next, x_next, k + 1)
			
		let loop_iter r p x iterations =  
			loop (r_curr, p_curr, x_curr, _) = (r, p, x, 0) for _i < iterations do
				let (r_next, p_next, x_next) = inner r_curr p_curr x_curr
				in (r_next, p_next, x_next, 0)
	
		let r_curr = explicit_residual A x b
		let p_curr = r_curr
		let x_curr = x

		let (_, _, x, _) = match c 
			case #iter k -> loop_iter r_curr p_curr x_curr k
			case #residual r -> loop_residual r_curr p_curr x_curr r
		in x
}