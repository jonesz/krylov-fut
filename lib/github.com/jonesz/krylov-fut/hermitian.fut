import "../../diku-dk/complex/complex"

module type hermitian_matrix = {
	-- | The scalar type.
	type t
	-- | The type of 'n' times 'n' hermitian matrices.
	type~ mat[n]

	-- | `idx (i, j) m` produces the element at logical position
	-- (i, j) in the symmetric matrix `m`.
	val idx [n] : (i64, i64) -> mat[n] -> t
		
	-- | Construct a hermitian matrix from a dense array.
	val herm [n] : [n][n]t -> mat[n]

	-- | Scale elements. Given a matrix and a scale value 'v', the
	-- function returns a new matrix with the elements scaled by 'v'.
	val scale [n] : t -> mat[n] -> [n][n]t

	-- | Map a function across the elements of the matrix.
	val matmap [n] : (t -> t) -> mat[n] -> [n][n]t

	-- | hermitian, hermitian matrix multiplication.
	val ssmm [n] : mat[n] -> mat[n] -> [n][n]t
	
	-- | hermitian, dense matrix multiplication.
	val sdmm [n][p] : mat[n] -> [n][p]t -> [n][p]t

	-- | dense, hermitian matrix multiplication.
	val dsmm [n][p] : [p][n]t -> mat[n] -> [p][n]t

	-- | Compute the conjugate transpose of a hermitian matrix.
	val conj_transpose [n] : mat[n] -> mat[n]
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

module mk_hermitian_matrix (T: complex) (R: ranking): hermitian_matrix with t = T.t = {
	type t = T.t

	type~ mat [n] = 
		?[nnz]. {
			size: [0][n](),
			data: [nnz]t
		}

	def idx [n] (i, j) (herm: mat[n]) =
		let (a, b) = if i > j then (j, i) else (i, j)
		let dt = #[unsafe] herm.data[R.rank (a, b)]
		in if i > j then T.conj dt else dt

	def herm [n] (arr: [n][n]t): mat[n] =
		{ size = [],
		  data = tabulate (elements n) (\p -> let (i, j) = R.unrank p in #[unsafe] arr[i, j])
		}

	def scale [n] (x: t) (A: mat[n]): [n][n]t = 
		map (\i -> map (\j -> (idx (i, j) A) |> (T.*) x) (iota n)) (iota n)

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

	def matmap [n] f (A: mat[]): [n][n]t =
		map (\i -> map (\j -> (idx (i, j) A) |> f) (iota n)) (iota n)

	def conj_transpose [n] (A: mat[n]): mat[n] =
		A with data = map T.conj A.data
}

module mk_hermitian_matrix_def (T: complex) =
	mk_hermitian_matrix T {
		def rank (i, j) = 
			elements i + j

		def unrank (p : i64) =
			let i = row p
			let j = p - elements i
			in (j, i)
	}
