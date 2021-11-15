# BZIntegral

This package is a Julia implementation of the recursive hybrid tetrahedron method (to be published) for Brillouin-zone (BZ) integration. The package is intended to provide systematic solutions to BZ integrals encountered in practical condensed matter calculation that converge slowly due to known singularites. 



A typical BZ integral may looks like

<!-- $ \text{Int} = \int d^d k W(\mathbf{k}) F(\mathbf{k})$ --> <img style="transform: translateY(0.1em); background: white;" src="svg\nVEwBt5r48.svg">

with a regular part <!-- $F(\mathbf{k})$ --> <img style="transform: translateY(0.1em); background: white;" src="svg\0XntFncjaW.svg"> and a singular weight function <!-- $W(\mathbf{k})$ --> <img style="transform: translateY(0.1em); background: white;" src="svg\JzfadJb3gv.svg">. <!-- $W(\mathbf{k})$ --> <img style="transform: translateY(0.1em); background: white;" src="svg\KLcG9sLLIN.svg"> may be discontinous or divergent in certain place in BZ, which causes poor convergence in integration. Typical example includes 
<!-- $\begin{aligned}
W(\mathbf{k}) &= \Theta(\epsilon_F-\varepsilon(\mathbf{k}))\\
W(\mathbf{k}) &= \delta(\epsilon_F-\varepsilon(\mathbf{k}))\\
W(\mathbf{k}) &= \Theta(\epsilon_F-\varepsilon(\mathbf{k}))\frac{1}{D(\mathbf{k})}\\
W(\mathbf{k}) &= \Theta(\epsilon_F-\varepsilon(\mathbf{k}))\delta(D(\mathbf{k})).
\end{aligned}$ --> <img style="transform: translateY(0.1em); background: white;" src="svg\ZQbznJ5iXz.svg">
BZIntegral.jl takes in values of some interpolable function, e.g. <!-- $\varepsilon(\mathbf{k})$ --> <img style="transform: translateY(0.1em); background: white;" src="svg\vvgZ1nqylK.svg"> and <!-- $D(\mathbf{k})$ --> <img style="transform: translateY(0.1em); background: white;" src="svg\WWNYzbQ5WE.svg">, on a **regular grid** filling the BZ and returns the integral weights on the same <!-- $\mathbf{k}$ --> <img style="transform: translateY(0.1em); background: white;" src="svg\edaQHJxNsx.svg">-grid which can be used to integrate any regular <!-- $F(\mathbf{k})$ --> <img style="transform: translateY(0.1em); background: white;" src="svg\JuYSjlFSX4.svg"> over the BZ by forming a weighted sum on the regular grid.

Currently three types of singular factors in <!-- $W(\mathbf{k})$ --> <img style="transform: translateY(0.1em); background: white;" src="svg\c8NcJYlWRs.svg"> are supported: <!-- $\Theta(X(\mathbf{k}))$ --> <img style="transform: translateY(0.1em); background: white;" src="svg\dUCcDpKnpX.svg">,<!-- $\delta(X(\mathbf{k}))$ --> <img style="transform: translateY(0.1em); background: white;" src="svg\aGX8b962A9.svg">, and <!-- $1/D(\mathbf{k})$ --> <img style="transform: translateY(0.1em); background: white;" src="svg\F11quYO4A4.svg">. Routines that handles <!-- $W(\mathbf{k})$ --> <img style="transform: translateY(0.1em); background: white;" src="svg\eZNVK1oF42.svg"> with multiple singular factors are available.
## Installation
At this point, the package has not been registered with General Registry, and can be installed by entering

```
] add https://github.com/SelimLin/BZIntegral.jl
```

in julia REPL.

## Naming convention
The routines in this package are named as 
`Quad*DRule***`
For example, the routine for 
<!-- $W(\mathbf{k}) = \Theta(\epsilon_F-\varepsilon(\mathbf{k}))$ --> <img style="transform: translateY(0.1em); background: white;" src="svg\NNva7pArbi.svg">
in 3D is named 
`Quad3DRuleΘ` 

Symbol "Θ","δ",and "𝔇" stand for <!-- $\Theta(X(\mathbf{k}))$ --> <img style="transform: translateY(0.1em); background: white;" src="svg\q482B73qPV.svg">,<!-- $\delta(X(\mathbf{k}))$ --> <img style="transform: translateY(0.1em); background: white;" src="svg\wQnxzkoOBp.svg">, and <!-- $1/D(\mathbf{k})$ --> <img style="transform: translateY(0.1em); background: white;" src="svg\iWuCF7eKVN.svg"> respectively. And the routine that handles a weight function with multiple singular factor,
<!-- $W(\mathbf{k}) = \Theta(\epsilon_F-\varepsilon(\mathbf{k}))\frac{1}{D(\mathbf{k})},$ --> <img style="transform: translateY(0.1em); background: white;" src="svg\0sj115mHVJ.svg">
 in 2D is named 
`Quad2DRuleΘ𝔇`

"𝔇" is typed as `\frakD` in Julia, which reminds me of `\frac{1}{D}` in Latex. 

## Usage
Currently, the package provides the following routines:
```
Quad*DRuleΘ,Quad*DRuleδ,Quad*DRuleΘ𝔇,Quad*DRuleΘδ,Quad*DRuleΘΘ,Quad*DRuleΘΘ𝔇,Quad*DRuleδδ,Quad*DRuleΘΘδ
```
along with some broadening methods contained in submodule `BZInt3D` and `BZInt2D`.

A typical routine with function header
```
function Quad3DRuleΘ(Emesh,eF,iter=2)
```
is described as 
```
Recursive tetrahedron method for weight function W(k) = Θ(eF-E(k))
``` 
in documented description, and we call it as
```
Wmesh = Quad3DRuleΘ(Emesh,eF,iter)
```
The varible `Emesh` and `Wmesh` should be 3D-Arrays with odd number of entres in each dimension. They correspond to values on a mesh filling the BZ with open boundary condition (OBC), points on one boundary are identified with points on the oppsite boundary in case of periodic boundary condition (PBC).

And with <!-- $F(\mathbf{k})$ --> <img style="transform: translateY(0.1em); background: white;" src="svg\yvruFDAsO5.svg"> sampled on the same grid, `Fmesh`, the integral can be calculated as 
```
res = sum(Fmesh.*Wmesh)*vol
```
`vol` is the volume of BZ.

`iter` is the variable that controls the number of grid refinements. If the algorithm is working properly, we should expect similar accuracy of results using <!-- $N\times N \times N$ --> <img style="transform: translateY(0.1em); background: white;" src="svg\GzGOSrMLAB.svg"> grid and `iter=n` with that using <!-- $2^nN\times 2^nN \times 2^nN$ --> <img style="transform: translateY(0.1em); background: white;" src="svg\JWW98BHZ4Y.svg"> and `iter=0`. Thus with larger `iter`, we can reproduce results on a dense grid without our recursive refinement using a considerably coarse grid.

## Parallelism
The package has been parallelled using `@threads`, the internal multi-threading technique provided by Julia. Package `StaticArrays` is used to avoid massive dynamical allocation on heaps during grid refinement, so the efficiency of the threading will not be impaired even in case of a large `iter`. 