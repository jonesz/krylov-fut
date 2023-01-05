-- | Sparse symmetric matrices, tailored for use in krylov methods.
--
-- A symmetric matrix is a square matrix where A = A^T.

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
	val map [n] : (t -> t) -> mat[n] -> mat[n]

		-- A * x
	val mul_vec [n] : mat[n] -> [n]t -> [n]t

	-- | Matrix multiplication.
	val smm [n] : mat[n] -> mat[n] -> [n][n]t

	-- | Matrix multiplication with a dense matrix.
	val smm_dense [n][p] : mat[n] -> [n][p]t -> [n][p]t
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

module mk_symmetric_matrix (T: numeric) (R: ranking): symmetric_matrix with t = T.t = {
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

	-- A * x, x * A
	-- NOTE: xA and Ax are just transposes of each other.
	def mul_vec [n] (A: mat[n]) (x: [n]t): [n]t =
		map (\i -> 
			map (\j -> 
				T.((idx (i, j) A) * x[j])) 
			(iota n) |> reduce_comm (T.+) (T.i64 0)) 
			(iota n)

	def smm [n] (A: mat[n]) (x: mat[n]): [n][n]t =
		???
	
	def smm_dense [n][p] (A: mat[n]) (x: [n][p]t): [n][p]t =
		???

	def map f (sym: mat[]) =
		sym with data = map f sym.data
}

module mk_symmetric_matrix_def (T: numeric) = 
	mk_symmetric_matrix T {
		def rank (i, j) =
			elements i + j

		def unrank (p : i64) =
			let i = row p
			let j = p - elements i
			in (j, i)
	}
