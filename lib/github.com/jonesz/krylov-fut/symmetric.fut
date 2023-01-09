-- | Sparse symmetric matrices, tailored for use in krylov methods.
--
-- A symmetric matrix is a square matrix where A = A^T.

import "../../diku-dk/linalg/linalg"

-- | The module type of a symmetric matrix.
module type symmetric_matrix = {
	-- | The scalar type.
	type t
	-- | The type of 'n' times 'n' symmetric matrices.
	type~ mat[n]

	-- | `idx (i, j) m` produces the element at logical position
	-- (i, j) in the symmetric matrix `m`.
	val idx [n] : (i64, i64) -> mat[n] -> t

	-- | Construct a symmetric matrix from a dense array.
	val sym [n] : [n][n]t -> mat[n]

	-- | Scale elements. Given a matrix and a scale value 'v', the
	-- function returns a new matrix with the elements scaled by 'v'.
	val scale [n] : t -> mat[n] -> mat[n]

	-- | Map a function across the elements of the matrix.
	val matmap [n] : (t -> t) -> mat[n] -> mat[n]

	-- | symmetric, symmetric matrix multiplication.
	val ssmm [n] : mat[n] -> mat[n] -> [n][n]t

	-- | symmetric, dense matrix multiplication.
	val sdmm [n][p] : mat[n] -> [n][p]t -> [n][p]t

	-- | dense, symmetric matrix multiplication.
	val dsmm [n][p] : [p][n]t -> mat[n] -> [p][n]t
}

-- The number of unique elements for a symmetric `n by n` array.
local def elements (n: i64) =
	(n * (1 + n)) / 2

-- See note within `diku-dk/sparse/triangular.fut`
local def row (i: i64) =
  i64.f64 (f64.ceil ((f64.sqrt(f64.i64(9+8*i))-1)/2))-1

-- Conversions between 1D <-> 2D.
local module type ranking = {
	val rank : (i64, i64) -> i64
	val unrank : i64 -> (i64, i64)
}

module mk_symmetric_matrix (T: field) (R: ranking): symmetric_matrix with t = T.t = {
	type t = T.t

	type~ mat [n] =
		?[nnz]. {
			size: [0][n](),
			data: [nnz]t
		}

	def idx [n] (i, j) (sym: mat[n]) =
		let (i, j) = if i > j then (j, i) else (i, j)
		in #[unsafe] sym.data[R.rank (i, j)]

	def sym [n] (arr: [n][n]t): mat[n] =
		{ size = [],
		  data = tabulate (elements n) (\p -> let (i, j) = R.unrank p in #[unsafe] arr[i, j])
		}

	def scale [n] (x: t) (A: mat[n]): mat[n] =
		A with data = map (T.* x) A.data

	def ssmm [n] (A: mat[n]) (B: mat[n]): [n][n]t =
		let row M i = map (\j -> idx (i, j) M) (iota n)
		let col M i = map (\j -> idx (j, i) M) (iota n)

		in map (\i -> 
			let r = row A i
			in map (\j -> 
				let c = col B j
				in reduce (T.+) (T.i64 0) <| map2 (T.*) r c) 
			(iota n))
		(iota n)	

	def sdmm [n][p] (A: mat[n]) (B: [n][p]t): [n][p]t =
		let row M i = map (\j -> idx (i, j) M) (iota n)
		let col M i = map (\j -> M[j, i]) (iota n)

		in map (\i -> 
			let r = row A i
			in map (\j -> 
				let c = col B j
				in reduce (T.+) (T.i64 0) <| map2 (T.*) r c) 
			(iota p))
		(iota n)

	def dsmm [n][p] (A: [p][n]t) (B: mat[n]): [p][n]t =
		let row M i = map (\j -> M[i, j]) (iota n)
		let col M i = map (\j -> idx (j, i) M) (iota n)

		in map (\i -> 
			let r = row A i
			in map (\j ->
				let c = col B j
				in reduce (T.+) (T.i64 0) <| map2 (T.*) r c)
					
				(iota n))
			(iota p)

	def matmap f (sym: mat[]) =
		sym with data = map f sym.data
}

module mk_symmetric_matrix_def (T: field) = 
	mk_symmetric_matrix T {
		def rank (i, j) =
			elements i + j

		def unrank (p : i64) =
			let i = row p
			let j = p - elements i
			in (j, i)
	}
