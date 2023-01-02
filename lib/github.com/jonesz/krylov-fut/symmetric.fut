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

module mk_symmetric_matrix (T: real) (R: ranking) = {
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
}

module mk_symmetric_matrix_def (T: real) = 
	mk_symmetric_matrix T {
		def rank (i, j) =
			elements i + j

		def unrank (p : i64) =
			let i = row p
			let j = p - elements i
			in (j, i)
	}
