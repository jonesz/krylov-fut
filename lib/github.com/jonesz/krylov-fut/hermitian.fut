import "symmetric"
import "complex"

module type hermitian_matrix = {
	include symmetric_matrix
	
	val transpose [n] : mat[n] -> mat[n]
}

module mk_hermitian_matrix (T: complex) (R: ranking): hermitian_matrix with t = T.t = {
	type t = T.t

	def transpose [n] (A: mat[n]): mat[n] =
		A with data = map T.conj A.data
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

module mk_hermitian_matrix_def (T: complex) =
	mk_hermitian_matrix T {
		def rank (i, j) = 
			elements i + j

		def unrank (p : i64) =
			let i = row p
			let j = p - elements i
			in (j, i)
	}
