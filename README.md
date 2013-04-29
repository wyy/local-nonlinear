# Usage

    matrix([OPTIONS])
    matrix('post', NODE, DOF)

## Description

- OPTIONS
  - `'irs'` or `'lr'`
    Model reduction using IRS method (default), or model reduction
    based on local rigid body mode.
  - `'linear'` or `'nonlinear'`
    Linear transient analysis (default), or local nonlinear transient
    analysis.
- NODE
  - Node number.
- DOF
  - (1, 2, 3, 4, 5, 6) for (UX, UY, UZ, ROTX, ROTY, ROTZ).

## Examples

`matrix()` 'irs' 'linear'

`matrix('nonlinear')` 'irs' 'nonlinear'

`matrix('lr')` 'lr' 'linear'

`matrix('lr', 'nonlinear')` 'lr' 'nonlinear'

`matrix('post', 191, 3)` specifies nodal data from the results file
