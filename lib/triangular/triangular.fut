-- Taken from: https://futhark-lang.org/examples/triangular.html

module type triangular = {
  type~ triangular [n] 'a

  val get [n] 'a :
    a -> (i64, i64) -> triangular [n] a -> a

  val from_array [n] 'a :
    [n][n]a -> triangular [n]a

  val to_array [n] 'a :
    a -> triangular [n] a -> [n][n]a

  val map [n] 'a 'b :
    (a -> b) -> triangular [n] a -> triangular [n] b
}

module triangular : triangular = {
  type~ triangular [n] 'a =
    { 
      size: [n][0](),
      data: []a
    }

    def elements (n: i64) = (n * (1 + n)) / 2

    def row (i: i64) =
      i64.f64 (f64.ceil ((f64.sqrt(f64.i64(9+8*i))-1)/2))-1

    def rank i j = elements i + j

    def unrank (p: i64) =
      let i = row p
      let j = p - elements i
      in (i, j)

    def get [n] 'a zero (i, j) (tri: triangular [n] a) =
      if j > i then zero else tri.data[rank i j]

    def from_array arr = 
      { 
        size = map (const []) arr,
        data = map (\p -> let (i, j) = unrank p
                          in arr[i, j])
                   (iota (elements (length arr)))
      }

    def to_array [n] 'a zero (tri: triangular [n] a) =
      tabulate_2d n n (\i j -> get zero (i, j) tri)

    def map f {size, data} =
      {size, data = map f data}
}
