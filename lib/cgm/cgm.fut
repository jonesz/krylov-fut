-- 2022, Ethan Jones <etn.jones@gmail.com>.

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

-- ==
-- entry: cgm
-- input { [[ 4f32, 1f32 ], [ 1f32, 3f32 ]] [[ 1f32 ], [ 2f32 ]] [[ 2f32 ], [ 1f32 ]] }
-- output { [[ 0.23564959f32, 0.3383686f32 ]] }
entry cgm [n] (A: [n][n]f32) (b: [n][1]f32) (x: [n][1]f32): [n][1]f32 =

  -- residual := b - Ax
  let residual = map2 (map2 (-)) b <| matmul_f32 A x
  let p0 = residual

  -- ak := (r^T * r) / p^T * A * p
  let alpha_k = 
    let numerator = matmul_f32 (transpose residual) residual
    let denominator = matmul_f32 (matmul_f32 (transpose p0) A) p0
    in numerator[0, 0] / denominator[0, 0] -- Both of these matrices should be [1][1].

  let xk = map2 (map2 (+)) x <| map (map (* alpha_k)) p0
  -- let rk = map2 (map2 (-)) residual <| matmul_f32 (map (map (* alpha_k)) A) p0

  in xk

  -- TODO: This is a single iteration: we can calculate the beta vector to the next direction.
  -- Ultimately, we need to determine what "residual" is and introduce a loop here.
