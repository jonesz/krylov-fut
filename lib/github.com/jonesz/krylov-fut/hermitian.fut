import "symmetric"
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
