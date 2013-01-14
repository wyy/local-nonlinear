# Usage

    matrix([OPTIONS])
    matrix('post', NODE, DOF)

## Description

- OPTIONS
  - `'linear'` or `'nonlinear'`
    Linear transient analysis (default), or local nonlinear transient
    analysis.
  - `'irs'` or `'lr'`
    Model reduction using IRS method (default), or model reduction
    based on local rigid body mode.
- NODE
  - Node for which data are to be stored.
- DOF
  - (1, 2, 3, 4, 5, 6) for (UX, UY, UZ, ROTX, ROTY, ROTZ).

## Examples

`matrix()` 'linear' 'irs'

`matrix('nonlinear')` 'nonlinear' 'irs'

`matrix('lr')` 'linear' 'lr'

`matrix('nonlinear', 'lr')` 'nonlinear' 'lr'

`matrix('post', 1646, 3)` specifies nodal data from the results file
