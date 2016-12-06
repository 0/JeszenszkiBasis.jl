# JeszenszkiBasis.jl

Bosonic occupation basis using algorithms from [Szabados et al., 2012](http://dx.doi.org/10.1016/j.chemphys.2011.10.003) ([preprint](http://coulson.chem.elte.hu/surjan/PREPRINTS/181.pdf)).

Tested with Julia 0.5.


## Installation

1. `Pkg.clone("https://github.com/0/JeszenszkiBasis.jl.git")`


## Examples

```julia
using JeszenszkiBasis
```

```julia
### 2 sites, 3 particles
basis = Szbasis(2, 3)
println(join([join(v, " ") for v in basis], ", "))
#-> 3 0, 2 1, 1 2, 0 3
println(length(basis))
#-> 4
```

```julia
### 4 sites, 4 particles
basis = Szbasis(4, 4)
v = basis.vectors[:, 8]
println(join(v, " "))
#-> 1 2 1 0
println(serial_num(basis, v))
#-> 8
println(sub_serial_num(basis, v[1:2]))
#-> 9
```

```julia
### 3 sites, 3 particles, 2 maximum
basis = RestrictedSzbasis(3, 3, 2)
println(join([join(v, " ") for v in basis], ", "))
#-> 2 1 0, 1 2 0, 2 0 1, 1 1 1, 0 2 1, 1 0 2, 0 1 2
println([2, 1, 0] in basis)
#-> true
println([3, 0, 0] in basis)
#-> false
```


## Caveats

* Iteration reuses the same vector for each step:

  ```julia
  basis = Szbasis(2, 1)
  println(join([v for v in basis], ", "))
  #-> [0,1], [0,1]
  ```


## Testing

Run all the tests:
```
JULIA_LOAD_PATH=src julia --color=yes test/runtests.jl
```


## Acknowledgements

Thanks to [Roger Melko](http://www.science.uwaterloo.ca/~rgmelko/) for getting me up to speed and providing a reference implementation!


## License

Provided under the terms of the MIT license.
See `LICENSE` for more information.
