# JeszenszkiBasis.jl

Bosonic occupation basis using algorithms from [Szabados et al., 2012](http://dx.doi.org/10.1016/j.chemphys.2011.10.003) ([preprint](http://coulson.chem.elte.hu/surjan/PREPRINTS/181.pdf)).


## Installation

1. `Pkg.clone("https://github.com/0/JeszenszkiBasis.jl.git")`


## Examples

```julia
using JeszenszkiBasis
```

```julia
### 2 sites, 3 particles
basis = Szbasis(2, 3)
println(join([join(basis.vectors[:, i], " ") for i=1:basis.D], ", "))
#-> 3 0, 2 1, 1 2, 0 3
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


## Acknowledgements

Thanks to [Roger Melko](http://www.science.uwaterloo.ca/~rgmelko/) for getting me up to speed and providing a reference implementation!


## License

Provided under the terms of the MIT license.
See `LICENSE` for more information.
