# An asymptotic study of nonlinear instability to Langmuir circulation in stratified shallow layers

**Author:** D. T. Hayes  
**Date:** July 20, 2020

> Converted from PDF to Markdown using text-preserving extraction. Display equations are kept in fenced text blocks where helpful to preserve structure for LLM parsing.

<!-- Page 1 -->

An asymptotic study of nonlinear instability to Langmuir circulation in stratified shallow layers D. T. Hayes

July 20, 2020

## Abstract

The CL equations governing instability to Langmuir circulation (LC) are solved by three approximate methods, viz: a small-l asymptotic expansion where l is the spanwise wavenumber, a power series method, and a Galerkin method. Interest is focused on the CL2 instability mechanism to LC and how it is influenced by stratification throughout the layer in which LC live. Some results are provided to illustrate the CL2 instability and how it is affected by nonlinearities.

## 1 Introduction

Langmuir circulation (LC) is a system of counter-rotating vortices that forms below wind driven waves in the upper ocean when the wind speed exceeds 3 m/s (Leibovich, 1983) and occupies the region of fluid that is sheared by the wind. Moreover LC is made visible by its surface footprints as almost parallel streaks or windrows on the ocean surface, with spacings of up to hundreds of meters (Plueddemann et al., 1996; Thorpe, 2004) and can extend for several kilometers in the direction of the wind (Thorpe, 2004). LC helps mix and form a region called the mixed layer (Langmuir, 1938) and in doing so this alters the variation with depth of density and temperature (Smith, 1992), on occasion to such a degree that the bottom of the layer is defined by a sharp change in temperature (density) known as a thermocline (pycnocline). Of interest in the present study is the role stratification plays on the evolution of LC in layers bounded by a thermocline. The prevailing theory for LC is due to Craik & Leibovich (1976), who derived a set of evolution equations to describe them known as the CL equations. Two instability mechanisms to excite LC follow from the CL equations (Leibovich, 1980) and both rest upon the interaction between shear U ′ in the surface layer resulting from the wind and differential Lagrangian drift D′ that results from the wave field. They are denoted CL1 and CL2. However CL2, which assumes that the drift does not vary cross stream to the wind, is considered the more likely instability to occur in Nature and is the mechanism studied in this paper. Of course, to ensure the problem is well posed, boundary conditions must be specified at the free surface and some distance below it. Neumann conditions are an obvious choice but, when imposed on finite layers as opposed to infinite ones (in the sense of deep water waves), the linear least stable wavenumber lc in usual circumstances is zero. This oddity was explained by Cox & Leibovich (1993), who noted that Neumann conditions ignore coupling between the perturbation flow and the extra stress it produces, implying that mixed boundary conditions that reflect that extra stress of magnitude γ, should be imposed. In doing so they found that lc is nonzero when γ is nonzero and also that lc ≪ 1 when γ ≪ 1. In view of that they chose to use perturbation methods to study the instability to CL2 of the simplest case U ′ = D′ = 1 in the small l limit, followed by Hayes & Phillips (2016) who allowed D′ and U ′ to be arbitrary functions of depth, while Hayes & Phillips (2017) studied the role of nonlinearities. In fact Cox & Leibovich (1993) also allowed for thermal stratification of slope H ′ = 1 and magnitude S whereas Hayes & Phillips (2016, 2017) set

```text
S = 0.
```

Our object here is to consider the role of nonlinearities when S is nonzero and D′ , U ′ , H ′ are arbitrary functions of depth. The governing equations are stated in §2. Then in §3 and §4 we outline solution methods for the respective linearised and nonlinear problems. In §5 we discuss our results on how the parameters and nonlinearities affect the CL2 instability to LC in the small l limit when S is nonzero. This includes a revisit of Figure 3 of Cox & Leibovich (1993). In §6 we conclude this paper and discuss some possibilities for further work on LC.

## 2 Problem description

### 2.1 CL2 equations

The CL2 equations follow from perturbations to the CL equations, where the perturbation velocity u = (u, v, w) and perturbation temperature θ are each defined for position x = (x, y, z) and time t. We take the x axis to be in the direction of the imposed shear, the y axis is in the spanwise direction, and the z axis is in the vertical direction. The flow is assumed to be independent of x. In dimensionless form, the CL2 equations are then (Craik & Leibovich, 1976)

```text
                                         ∂                  ∂u   ∂θ
                                     (      - ∇2 )∇2 ψ = RD′ - S    + J(ψ, ∇2 ψ),                                  (2.1.1)
                                         ∂t                 ∂y   ∂y
                                                  ∂            ∂ψ ′
                                              (      - ∇2 )u =    U + J(ψ, u),                                     (2.1.2)
                                                  ∂t           ∂y
                                                  ∂             ∂ψ ′
                                              (      - τ∇2 )θ =    H + J(ψ, θ)                                     (2.1.3)
                                                  ∂t            ∂y
```

<!-- Page 2 -->

where J is the Jacobian J(a, b) = ay bz - az by . To satisfy the continuity equation the stream function ψ is defined by

```text
v = ψz and w = -ψy . We further have that U ′ and H ′ must satisfy
```

∂ ∂

```text
                                              (      - ∇2 )U = F,              (      - ∇2 )H = G                                       (2.1.4)
                                                  ∂T                               ∂T
```

where T and t are disparate time scales and F, G are due to body forces and heat sources respectively. The differential drift D′ results from the Stokes drift whose details depend on the wavefield. We can thus take D′ , U ′ , H ′ to be arbitrary functions of z. Herein we allow D′ , U ′ , H ′ to each be arbitrary polynomials of z N X N X N X

```text
                                      D′ =         A n zn , U ′ =          Bn zn , and H ′ =                   Cn zn                    (2.1.5)
                                             n=0                     n=0                                 n=0
```

where An , Bn and Cn are arbitrary constant coefficients. The Rayleigh number is denoted by R, the magnitude of the stratification is denoted by S , and τ , 0 is an inverse Prandtl number. Nonlinearities are accounted for through the Jacobian J. When nonlinearities are assumed to be small we discard J to yield the linearised CL2 equations. When

```text
S = 0, equations (2.1.1), (2.1.2) are those used in Hayes & Phillips (2017).
```

### 2.2 Boundary conditions

We will use mixed boundary conditions on the top and bottom of the layer of fluid that are similar to those introduced by Cox & Leibovich (1993)

```text
                                 ∂2 ψ      ∂ψ ∂u            ∂θ
                                      + γ1    =    + γ2 u =    + β1 θ = ψ = 0 on z = 0,                                                 (2.2.1)
                                 ∂z2       ∂z   ∂z          ∂z

                                    ∂2 ψ      ∂ψ ∂u              ∂θ
                                         + γ3     =     + γ4 u =     + β2 θ = ψ = 0 on z = -1                           (2.2.2)
                                     ∂z2       ∂z    ∂z          ∂z
```

where γi , β j for i = 1, 2, 3, 4, j = 1, 2 are constants. We set z = 0 at the top of the layer and z = -1 at the bottom of the layer.

## 3 Linear methods

### 3.1 Linear perturbation solution

We seek a perturbation solution to the linearised version of the CL2 equations (2.1.1), (2.1.2), (2.1.3) with boundary conditions (2.2.1) and (2.2.2) using l ≪ 1 as a small parameter. This calculation is an extension of the work of Cox & Leibovich (1993) and Hayes & Phillips (2016). We assume ∞ X ∞ X ∞ X

```text
                         ψ=i          l2k+1 ψ2k+1 eily eσt , u =               l2k u2k eily eσt , θ =              l2k θ2k eily eσt ,   (3.1.1)
                                k=0                                   k=0                                    k=0
```

∞ X ∞ X

```text
                                                  σ=         l2k+2 σ2k+2 , R =                   l2k R2k ,                              (3.1.2)
                                                       k=0                                 k=0
                                                                                ∞
                                                                                X
                                                       γi = l4 γi , βi =                 βi,2k l2k .                                    (3.1.3)
                                                                                   k=0
```

Here ψ2k+1 , u2k , and θ2k are each functions of z. To proceed we substitute the above expansions (3.1.1-3.1.3) into the linear CL2 equations and boundary conditions, equate like powers of l, and solve the resulting equations at successive orders in l. To equate like powers of l we can use the Cauchy product formula (Hardy, 1949) ∞ X ∞ X m ∞ X X

```text
                                                    am x m         bn xn =                   am-n bn xm .                               (3.1.4)
                                             m=0             n=0               m=0 n=0
```

In §3.1.1 we take the calculation to O(l4 ). In §3.1.2 an algorithm is derived so we can take the perturbation solution to O(lP ) for any integer P > 0 within computational limits. The algorithm can then be coded into Maple. The solutions we obtain may then be used to validate more general numerical calculations. The ci and ci, j appearing here are given in the Appendix. It turns out that the linear perturbation solution separates into two separate cases. We have case I:

```text
-β1,0 + β2,0 β1,0 + β2,0 , 0 and case II: -β1,0 + β2,0 β1,0 + β2,0 = 0. This becomes evident when applying the boundary
```

conditions for θ0 .

<!-- Page 3 -->

#### 3.1.1 The first few orders

At O(l0 ) we have u′′

```text
                                                              0 =0                                                (3.1.5)
```

with boundary conditions

```text
                                                     u′0 = 0 on z = 0, -1.                                        (3.1.6)
```

So

```text
                                                             u0 = c0                                              (3.1.7)
```

where c0 is an arbitrary constant. Without affecting σ we may let

```text
                                                             u0 = 1.                                              (3.1.8)
```

Also at O(l0 ) we have

```text
                                                             τθ0′′ = 0                                            (3.1.9)
```

with boundary conditions

```text
                                                   θ0′ + β1,0 θ0 = 0 on z = 0                                   (3.1.10)
```

and

```text
                                                  θ0′ + β2,0 θ0 = 0 on z = -1.                                  (3.1.11)
```

For case I we find

```text
                                                          θ0 = c3 = 0.                                          (3.1.12)
```

For case II we find

```text
                                                  θ0 = c3 (-β1,0 z + 1) = c3 θ̂0                                (3.1.13)
```

where c3 is an arbitrary constant. At O(l1 ) we have ψ′′′′ ′

```text
                                                     1 = -D R0 u0 + S θ0                                        (3.1.14)
```

with boundary conditions ψ′′

```text
                                                   1 = ψ1 = 0 on z = 0, -1.                                     (3.1.15)
```

Solving (3.1.14), (3.1.15) gives & z & z z3 z2

```text
              ψ1 = -R0 u0               D′ dz dz dz dz + S     θ0 dz dz dz dz + c4 + c5 + c6 z + c7
                                                                                  6    2
                        = -R0 u0 ψ̃1 + S ψ̂1 c3                                                                 (3.1.16)
```

where ψ̃1 , ψ̂1 are devoid of R0 , S , and c3 . Notice here that since multiple layers of LC occurred with D′ as a linear function of z in Hayes & Phillips (2016), this then means that here even with D′ = U ′ = H ′ = 1 we can have multiple layers of LC due to θ0 being linear in z. At O(l2 ) we have u′′ ′

```text
                                                    2 = U ψ1 + u0 σ2 + u0                                        (3.1.17)
```

with boundary conditions

```text
                                                     u′2 = 0 on z = 0, -1.                                      (3.1.18)
```

Solving (3.1.17), (3.1.18) gives " z z2

```text
                                     u2 =         U ′ ψ1 dz dz + u0 (σ2 + 1)        + c8 z + c9 .               (3.1.19)
                                                                                 2
```

The constant of integration c9 is chosen so that there is no net flux of fluid due to the perturbation flow Z 0

```text
                                                           u2 dz = 0.                                           (3.1.20)
                                                            -1
```

Also at O(l2 ) we have

```text
                                                   τθ2′′ = ψ1 H ′ + θ0 (τ + σ2 )                                (3.1.21)
```

with boundary conditions

```text
                                              θ2′ + β1,0 θ2 + β1,2 θ0 = 0 on z = 0                              (3.1.22)
```

and

```text
                                             θ2′ + β2,0 θ2 + β2,2 θ0 = 0 on z = -1.                             (3.1.23)
```

We find " z " z H′ σ2

```text
                              θ2 =          ψ1    dz dz +         θ0 (      + 1) dz dz + c10 z + c11 .          (3.1.24)
                                               τ                         τ
```

<!-- Page 4 -->

For case II we have θ2 = θ̃2 + c11 θ̂2 where θ˜2 and θˆ2 are independent of c11 . The boundary conditions lead to equations for σ2 . For case I we find Z 0

```text
                                                   σ2 = R0           U ′ ψ̃1 dz - 1.                                    (3.1.25)
                                                                -1
```

For case II there are two equations involving σ2 which can be written as the matrix equation for v = (u0 , c3 )T

```text
                                                             Mv = 0                                                     (3.1.26)
```

where the elements of the 2 × 2 matrix M are given in the Appendix. For a nontrivial solution the determinant of M must be zero, M1,1 M2,2 - M1,2 M2,1 = 0. This leads to a quadratic equation

```text
                                                      aσ22 + bσ2 + c = 0                                                (3.1.27)
```

for σ2 where a, b, c are in the Appendix. The quadratic formula gives √ -b ± b2 - 4ac

```text
                                               σ2 =                   .                                                 (3.1.28)
                                                            2a
```

The matrix equation then yields M1,1 u0

```text
                                                        c3 = -               .                                          (3.1.29)
                                                                      M1,2
```

At O(l3 ) we have ψ′′′′ ′′ ′

```text
                                       3 = (2 + σ2 )ψ1 - (R2 u0 + R0 u2 )D + S θ2                                       (3.1.30)
```

with boundary conditions ψ′′

```text
                                                  3 = ψ3 = 0 on z = 0, -1.                                              (3.1.31)
```

Solving (3.1.30), (3.1.31) gives " z & z

```text
                          ψ3 = (2 + σ2 )     ψ1 dz dz -          D′ (R2 u0 + u2 R0 ) dz dz dz dz
                                   & z
                                                              z3        z2
                               +S         θ2 dz dz dz dz + c12 + c13 + c14 z + c15 .                                    (3.1.32)
                                                              6         2
```

For case I we have ψ3 = ψ̂3 - R2 ψ̃3 + S ψ̆3 where ψ̂3 , ψ̃3 , ψ̆3 are each independent of R2 and S . For case II we have ψ3 = ψ̂3 - R2 ψ̃3 + c11 ψ̆3 where ψ̂3 , ψ̃3 , ψ̆3 are each independent of R2 and c11 . The dependence on S appears too complicated to isolate for case II. At O(l4 ) we have u′′

```text
                                                  4 = u0 σ4 + u2 (1 + σ2 ) + ψ3 U
                                                                                  ′
                                                                                                                  (3.1.33)
```

with boundary conditions

```text
                                                   u′4 + γ2 u0 = 0 on z = 0                                             (3.1.34)
```

and

```text
                                                  u′4 + γ4 u0 = 0 on z = -1.                                            (3.1.35)
```

Solving (3.1.33), (3.1.34), (3.1.35) gives " z " z z2

```text
                         u4 =        u2 (1 + σ2 ) dz dz +     ψ3 U ′ dz dz + u0 σ4 + c16 z + c17                        (3.1.36)
                                                                                  2
```

where the constant of integration c17 is chosen so to exclude net mass transfer as before Z 0

```text
                                                         u4 dz = 0.                                                     (3.1.37)
                                                           -1
```

Also at O(l4 ) we have

```text
                                               τθ4′′ = ψ3 H ′ + (τ + σ2 )θ2 + θ0 σ4                                     (3.1.38)
```

with boundary conditions

```text
                                          θ4′ + β1,0 θ4 + β1,2 θ2 + β1,4 θ0 = 0 on z = 0                                (3.1.39)
```

and

```text
                                      θ4′ + β2,0 θ4 + β2,2 θ2 + β2,4 θ0 = 0 on z = -1.                                  (3.1.40)
```

We find " z " z " z " z H′ σ2 σ4

```text
          θ4 =         ψ3    dz dz +            θ2 dz dz +                θ0 dz dz +         θ2 dz dz + c18 z + c19 .   (3.1.41)
                          τ          τ                     τ
```

<!-- Page 5 -->

The boundary conditions lead to equations for σ4 . For case I we find Z 0 U′

```text
                                        σ4 = -            (ψ̂3 - R2 ψ̃3 + S ψ̆3 ) dz + γ4 - γ2 .                (3.1.42)
                                                    -1 u0
```

For case II there are two equations for σ4 which can be written as the matrix equation for w = (σ4 , c11 )T

```text
                                                                       Nw = q                                   (3.1.43)
```

where the elements of the 2 × 2 matrix N and the elements of the 2 × 1 vector q are given in the Appendix. Solving this matrix equation then yields N1,1 q2 - N2,1 q1

```text
                                            c11 =                                                          (3.1.44)
                                                   N2,2 N1,1 - N2,1 N1,2
```

and q1 N1,2 N1,1 q2 - N2,1 q1

```text
                                            σ4 =         -     (                     ).                         (3.1.45)
                                                     N1,1 N1,1 N2,2 N1,1 - N2,1 N1,2
```

This calculation recovers Cox & Leibovich (1993) on setting U ′ = D′ = H ′ = 1 and recovers Hayes & Phillips (2016) on setting S = 0.

#### 3.1.2 Linear perturbation solution algorithm

At O(l2 j ) for integer j ≥ 2 we have j-1 X u′′

```text
                                          2j =            u2 j-(2m+2) σ2m+2 + u2 j-2 + U ′ ψ2 j-1               (3.1.46)
                                                   m=0
```

with boundary conditions

```text
                                                     u′2 j + γ2 u2 j-4 = 0 on z = 0                             (3.1.47)
```

and

```text
                                                    u′2 j + γ4 u2 j-4 = 0 on z = -1.                            (3.1.48)
```

Consistent with our progression above we choose the constant of integration so that Z 0

```text
                                                                     u2 j dz = δ0, j u0                         (3.1.49)
                                                                -1
```

where δi, j is the Kronecker delta 

```text
                                                             δi, j =        1,   i= j                           (3.1.50)
                                                                            0,   i , j.
```

On solving for u2 j we find " zX j-1 " z

```text
                               u2 j =                       u2 j-(2m+2) σ2m+2 dz dz +            u2 j-2 dz dz
                                                 m=0
                                                " z
                                            +             U ′ ψ2 j-1 dz dz + c0, j z + c1, j .                  (3.1.51)
```

Also at O(l2 j ) for integer j ≥ 2 we have j-1

```text
                                                   X                       σ2m+2            H′
                                         θ2′′j =         θ2 j-(2m+2)             + θ2 j-2 +    ψ2 j-1           (3.1.52)
                                                   m=0
                                                                             τ              τ
```

with boundary conditions j X

```text
                                                θ2′ j +         β1,2m θ2 j-2m = 0 on z = 0                      (3.1.53)
                                                          m=0
```

and j X

```text
                                              θ2′ j +          β2,2m θ2 j-2m = 0 on z = -1.                     (3.1.54)
                                                         m=0
```

<!-- Page 6 -->

On solving for θ2 j we find " zX j-1 " z σ2m+2

```text
                                   θ2 j =                   θ2 j-(2m+2)       dz dz +                  θ2 j-2 dz dz
                                                        m=0
                                                                          τ
                                                  " z
                                                             H′
                                             +                  ψ2 j-1 dz dz + c2, j z + c3, j .                                         (3.1.55)
                                                             τ
```

For case II we have θ2 j = θ̃2 j + c3, j θ̂2 j where θ̃2 j and θ̂2 j are independent of c3, j . At O(l2 j+1 ) we have j X j-1 X ψ′′′′

```text
                                    2 j+1 = -                D′ u2 j-2m R2m + S θ2 j +           ψ′′
                                                                                                  2 j-2m-1 σ2m+2
                                                       m=0                                 m=0
                                                       j-2
                                                       X
                                                  -          ψ2 j-2m-3 σ2m+2 + 2ψ′′
                                                                                 2 j-1 - ψ2 j-3                                          (3.1.56)
                                                       m=0
```

with boundary conditions ψ′′ ′

```text
                                               2 j+1 + γ1 ψ2 j-3 = ψ2 j+1 = 0 on z = 0                                                   (3.1.57)
```

and ψ′′ ′

```text
                                              2 j+1 + γ3 ψ2 j-3 = ψ2 j+1 = 0 on z = -1.                                                  (3.1.58)
```

On solving for ψ2 j+1 we find & zX j " zX j-1 ′

```text
                ψ2 j+1 = -                        D u2 j-2m R2m dz dz dz dz +                         ψ2 j-2m-1 σ2m+2 dz dz
                                     m=0                                                       m=0
                                  & zX
                                     j-2                                                  & z
                              -                   ψ2 j-2m-3 σ2m+2 dz dz dz dz -                        ψ2 j-3 dz dz dz dz
                                            m=0
                                   " z                           & z
                                                                                                        z3        z2
                              +2         ψ2 j-1 dz dz + S                   θ2 j dz dz dz dz + c4, j       + c5, j + c6, j z + c7, j .   (3.1.59)
                                                                                                        6         2

For case II we have ψ2 j+1 = ψ̂2 j+1 - R2 j ψ̃2 j+1 + c3, j ψ̆2 j+1 where ψ̂2 j+1 , ψ̃2 j+1 , and ψ̆2 j+1 are independent of R2 j and
```

c3, j . The boundary conditions lead to equations for σ2 j . For case I we find Z 0 γ γ U′

```text
                                  σ2 j = - 2 u2 j-4 |z=0 + 4 u2 j-4 |z=-1 -                  ψ2 j-1 dz - δ0, j-1 .                       (3.1.60)
                                          u0              u0                           -1 u0
```

For case II there are two equations for σ2 j which can be written as a matrix equation for w = (σ2 j , c3, j-1 )T as

```text
                                                                       Nw = q                                                            (3.1.61)
```

where the elements of the 2 × 2 matrix N and the elements of the 2 × 1 vector q which here depend on j are given in the Appendix. Solving this matrix equation then yields N1,1, j q2, j - N2,1, j q1, j

```text
                                                    c3, j-1 =                                                                            (3.1.62)
                                                                  N2,2, j N1,1, j - N2,1, j N1,2, j
```

and q1, j N1,2, j N1,1, j q2, j - N2,1, j q1, j

```text
                                         σ2 j =           -        (                                 ).                                  (3.1.63)
                                                   N1,1, j N1,1, j N2,2, j N1,1, j - N2,1, j N1,2, j
```

The above calculation can then be coded in Maple within a loop. This calculation recovers Hayes & Phillips (2016) on setting S = 0.

### 3.2 Linear power series method

Consistent with the perturbation method above, we assume

```text
                                          ψ = iAeily eσt , u = Beily eσt , and θ = Ceily eσt                                              (3.2.1)
```

<!-- Page 7 -->

where A = A(z), B = B(z), and C = C(z). We substitute (3.2.1) into the linearised CL2 equations which leads to

```text
                                A′′′′ - (2l2 + σ)A′′ + (l4 + l2 σ)A + RD′ Bl - S Cl = 0,                         (3.2.2)
                                                  B′′ - (l2 + σ)B - lAU ′ = 0,                                   (3.2.3)
```

and

```text
                                            τC ′′ - (τl2 + σ)C - lAH ′ = 0.                                      (3.2.4)
```

Substituting (3.2.1) into the boundary conditions leads to

```text
                                A′′ + γ1 A′ = B′ + γ2 B = C ′ + β1C = A = 0 on z = 0                             (3.2.5)
```

and

```text
                               A′′ + γ3 A′ = B′ + γ4 B = C ′ + β2C = A = 0 on z = -1.                            (3.2.6)
```

In the linear power series method we let 4+M X 2+M X 2+M X

```text
                                 A=          am zm , B =             bm zm , and C =             cm zm .         (3.2.7)
                                      m=0                      m=0                         m=0
```

We then substitute (3.2.7) into the differential equations (3.2.2), (3.2.3), (3.2.4) and boundary conditions (3.2.5), (3.2.6). Equating like powers of z in accordance with Theorem A in the Appendix then leads to a set of algebraic equations. These algebraic equations can then be solved numerically. To produce some of the results, numerical methods were combined. All of our linear power series method codes used adaptive Newton’s method for systems of algebraic equations. When the minimum turning point on the neutral curve was required the Golden section algorithm was used. When finding the point where Re σ = 0 the bisection method was used. Since the numerical methods used here are iterative, rapid convergence depended upon initial guesses and the perturbation solution results shined light on appropriate initial guesses. For case II where σ has two branches, different complex valued initial guesses must be used to find both of the branches of σ. We here set b0 = 1.

### 3.3 Linear Galerkin method

As in the linear power series method, we seek solutions of the form

```text
                                    ψ = iAeily eσt , u = Beily eσt , and θ = Ceily eσt                           (3.3.1)
```

which leads, as above, to (3.2.2) through (3.2.6). However, we here express A, B, and C in terms of orthogonal basis functions premultiplied by coefficients M X M-2 X M-2 X

```text
                                A=          am Pm , B =              bm Pm , and C =             cm Pm .         (3.3.2)
                                     m=0                       m=0                         m=0
```

Here Pm = Pm (z) are shifted Legendre basis functions on z ∈ [-1, 0] defined by 1 dm

```text
                                            Pm =         {    [(x2 - 1)m ]}| x=2z+1                              (3.3.3)
                                                    2m m! dxm
```

and satisfy Z 0

```text
                                                               Pi P j dz ∝ δi, j .                               (3.3.4)
                                                          -1
```

We substitute expansions (3.3.1) into the linearised CL2 equations to obtain equations in z whose residuals r1 , r2 , r3 can be expanded as X∞ ∞ X ∞ X

```text
                                    r1 =     a⋆i Pi , r2 =   b⋆i Pi , r3 =   c⋆i Pi                            (3.3.5)
                                            i=0                      i=0                 i=0
```

where

```text
                                      a⋆i ∝ hr1 , Pi i, b⋆i ∝ hr2 , Pi i, c⋆i ∝ hr3 , Pi i.                      (3.3.6)
```

In the Galerkin method we require

```text
                                              hr1 , Pi i = hr2 , Pi i = hr3 , Pi i = 0.                          (3.3.7)
```

That is, we require Z 0 Z 0 Z 0

```text
                                            r1 P j dz =         r2 P j dz =          r3 P j dz = 0               (3.3.8)
                                       -1                  -1                   -1
```

for j = 0, 1, 2, . . . , M - 4, which yield algebraic equations for the unknown coefficients. The basis functions do not inherently satisfy the boundary conditions and extra equations are found by substituting into the boundary conditions. This technique is known as the tau-method. The resulting algebraic equations are then treated much the same as in the linear power series method. For consistency with the linear power series solutions we herein choose b0 so that

```text
B|z=0 = 1.
```

<!-- Page 8 -->

## 4 Nonlinear methods

### 4.1 Nonlinear perturbation solution

```text
We seek a nonlinear perturbation solution to the CL2 equations (2.1.1), (2.1.2), (2.1.3) with boundary conditions (2.2.1)
```

and (2.2.2) in the small l limit. This calculation is an extension of the work of Hayes & Phillips (2017). Consistent with the linear perturbation solution we write

```text
                 Y = ly, T = l2 t, ψ(y, z, t) = lΨ̄(Y, z, T ), u(y, z, t) = ū(Y, z, T ), θ(y, z, t) = Θ̄(Y, z, T ),      (4.1.1)
```

and ∞ X

```text
                                                             γi = l4 γi , βi =           βi,2k l2k .                      (4.1.2)
                                                                                   k=0
```

Equation (2.1.1) then becomes

```text
                                                     ∂ 2 ∂2 Ψ̄ ∂2 Ψ̄
                                                     l3 [l         + 2]
                                                    ∂T ∂Y 2            ∂z
                                                        ∂ 4 Ψ̄           ∂4 Ψ̄      ∂4 Ψ̄
                                                  -l[l4 4 + 2l2 2 2 + 4 ]
                                                        ∂Y             ∂Y ∂z        ∂z
                                                       ∂ū         ∂ Θ̄        ∂ Ψ̄ ∂     ∂2 Ψ̄ ∂2 Ψ̄
                                                = RD′ l - S             l + l3        [l2 2 + 2 ]
                                                       ∂Y          ∂Y          ∂Y ∂z ∂Y         ∂z
                                                      ∂Ψ̄    ∂     ∂ 2 Ψ̄    ∂ 2 Ψ̄
                                                  -l3          [l2        + 2]                                            (4.1.3)
                                                      ∂z ∂Y ∂Y 2              ∂z
```

while equation (2.1.2) becomes

```text
                                              ∂ū     ∂2 ū ∂2 ū ∂Ψ̄ ′ 2 ∂Ψ̄ ∂ū 2 ∂Ψ̄ ∂ū
                                         l2       - l2 2 - 2 = l2    U +l        -l                                       (4.1.4)
                                              ∂T      ∂Y    ∂z    ∂Y      ∂Y ∂z     ∂z ∂Y
```

and equation (2.1.3) becomes

```text
                                         ∂Θ̄      ∂2 Θ̄ ∂2 Θ̄   ∂Ψ̄ ′ 2 ∂Ψ̄ ∂Θ̄ 2 ∂Ψ̄ ∂Θ̄
                                    l2       - τl2 2 - τ 2 = l2    H +l        -l         .                               (4.1.5)
                                         ∂T       ∂Y    ∂z      ∂Y      ∂Y ∂z     ∂z ∂Y
```

The boundary conditions become ∞

```text
                              ∂2 Ψ̄ 5 ∂Ψ̄ ∂ū 4                 ∂Θ̄ X
                          l        + l γ 1    =    + l γ 2 ū =    +   β1,2m l2m Θ̄ = lΨ̄ = 0 on z = 0,                   (4.1.6)
                              ∂z2          ∂z   ∂z              ∂z m=0
                                                                 ∞
                          ∂2 Ψ̄ 5 ∂Ψ̄ ∂ū 4                 ∂Θ̄ X
                      l        + l γ 3    =    + l γ 4 ū =    +   β2,2m l2m Θ̄ = lΨ̄ = 0 on z = -1.                      (4.1.7)
                          ∂z2          ∂z   ∂z              ∂z m=0
```

We let ∞ X ∞ X ∞ X

```text
                                              Ψ̄ =         Ψ2k l2k , ū =         u2k l2k , Θ̄ =             Θ2k l2k ,    (4.1.8)
                                                     k=0                    k=0                        k=0
```

and ∞ X

```text
                                                                    R=            R2k l2k                                 (4.1.9)
                                                                            k=0
```

where Ψ2k , u2k , and Θ2k , are functions of Y, z, and T . These expansions are consistent with those from the linear perturbation solution. We substitute (4.1.8) and (4.1.9) into (4.1.3 - 4.1.7) and equate like powers of l using the Cauchy product formula. At O(l2k ) we have

```text
                    ∂                           ∂2             ∂2
                      u2(k-1) -                     u2(k-1) -      u2k                                                   (4.1.10)
                   ∂T                          ∂Y 2            ∂z2
                                                                 k-1                     k-1
                                                ∂                X     ∂          ∂      X   ∂            ∂
                                     =             Ψ2(k-1) U ′ +         Ψ2(k-m-1) u2m -        Ψ2(k-m-1) u2m
                                               ∂Y                m=0
                                                                     ∂Y           ∂z     m=0
                                                                                             ∂z          ∂Y
```

with boundary conditions ∂

```text
                                                              u2k + γ2 u2(k-2) = 0 on z = 0,                             (4.1.11)
                                                           ∂z
```

<!-- Page 9 -->

∂

```text
                                                        u2k + γ4 u2(k-2) = 0 on z = -1,                          (4.1.12)
                                                     ∂z
```

and

```text
                     ∂             ∂2              ∂2
                       Θ2(k-1) - τ 2 Θ2(k-1) - τ 2 Θ2k                                                           (4.1.13)
                    ∂T            ∂Y               ∂z
                                                  k-1                    k-1
                                  ∂               X    ∂          ∂      X   ∂            ∂
                               =    Ψ2(k-1) H ′ +        Ψ2(k-m-1) Θ2m -        Ψ2(k-m-1) Θ2m
                                 ∂Y               m=0
                                                      ∂Y          ∂z     m=0
                                                                             ∂z          ∂Y
```

with boundary conditions k ∂ X

```text
                                                   Θ2k +     β1,2m Θ2k-2m = 0 on z = 0,                          (4.1.14)
                                                ∂z       m=0
                                                         k
                                               ∂        X
                                                  Θ2k +     β2,2m Θ2k-2m = 0 on z = -1.                          (4.1.15)
                                               ∂z       m=0
```

At O(l2k+1 ) we have

```text
                                    ∂ ∂2                ∂ ∂2
                                          2
                                            Ψ2(k-2) +         Ψ2(k-1)
                                   ∂T ∂Y               ∂T ∂z2
                                      ∂4                 ∂4            ∂4
                                   - 4 Ψ2(k-2) - 2 2 2 Ψ2(k-1) - 4 Ψ2k
                                     ∂Y               ∂Y ∂z            ∂z
                                    k                                  k-2
                                   X               ∂           ∂       X    ∂           ∂ ∂2
                                 =     R2(k-m) D′ u2m - S        Θ2k +        Ψ2(k-m-2)         Ψ2m
                                   m=0
                                                  ∂Y          ∂Y       m=0
                                                                           ∂Y           ∂z ∂Y 2
                                        k-1                               k-2
                                        X   ∂                  ∂3         X   ∂              ∂3
                                    +              Ψ2(k-m-1)       Ψ 2m -        Ψ 2(k-m-2)      Ψ2m
                                        m=0
                                              ∂Y               ∂z3        m=0
                                                                              ∂z            ∂Y 3
                                        k-1
                                        X   ∂             ∂ ∂2
                                    -          Ψ2(k-m-1)        Ψ2m                                              (4.1.16)
                                        m=0
                                            ∂z           ∂Y ∂z2
```

with boundary conditions ∂2 ∂

```text
                                                   Ψ2k + γ1 Ψ2(k-2) = Ψ2k = 0 on z = 0,                          (4.1.17)
                                              ∂z           ∂z
                                      ∂2          ∂
                                        2
                                          Ψ2k + γ3 Ψ2(k-2) = Ψ2k = 0 on z = -1.                                   (4.1.18)
                                     ∂z           ∂z
```

The equations above are to be solved for every integer k ⩾ 0. As in the linear perturbation solution, the nonlinear perturbation solution separates into two separate cases, that is case I: -β1,0 + β2,0 β1,0 + β2,0 , 0 and case II: -β1,0 +

```text
β2,0 β1,0 + β2,0 = 0.
```

#### 4.1.1 The first few orders

At O(l0 ) we have ∂2 u0

```text
                                                                          =0                                     (4.1.19)
                                                                     ∂z2
```

with boundary conditions ∂u0

```text
                                                              = 0 on z = 0, -1.                                  (4.1.20)
                                                           ∂z
```

The solution to this problem is

```text
                                                                 u0 = u0 (Y, T )                                 (4.1.21)
```

where u0 (Y, T ) is arbitrary. Also at O(l0 ) we have ∂2 Θ0

```text
                                                                  τ         =0                                   (4.1.22)
                                                                       ∂z2
```

with boundary conditions ∂Θ0

```text
                                                           + β1,0 Θ0 = 0 on z = 0,                               (4.1.23)
                                                        ∂z
```

<!-- Page 10 -->

∂Θ0

```text
                                                  + β2,0 Θ0 = 0 on z = -1.                                           (4.1.24)
                                               ∂z
```

For case I the solution to this problem is

```text
                                                     Θ0 = c3 (Y, T ) = 0.                                            (4.1.25)
```

For case II the solution to this problem is

```text
                                          Θ0 = c3 (Y, T )(-β1,0 z + 1) = c3 (Y, T )Θ̂0                               (4.1.26)
```

where c3 (Y, T ) is arbitrary. At O(l1 ) we have

```text
                                                ∂4 Ψ0          ∂u0    ∂Θ0
                                                      = -R0 D′     +S                                                (4.1.27)
                                                 ∂z4           ∂Y     ∂Y
```

with boundary conditions ∂2 Ψ0

```text
                                                   = Ψ0 = 0 on z = 0, -1.                                            (4.1.28)
                                              ∂z2
```

The solution to this problem can be found to be & z & z

```text
                                      ∂u0             ′                          ∂Θ0
                         Ψ0 = -R0                   D dz dz dz dz + S                dz dz dz dz
                                      ∂Y                                          ∂Y
                                            z3            z2
                                 +c4 (Y, T ) + c5 (Y, T ) + c6 (Y, T )z + c7 (Y, T )
                                            6              2
                                         ∂u0         ∂c3 (Y, T )
                                 = -R0       Ψ̃0 + S             Ψ̂0                                                 (4.1.29)
                                         ∂Y             ∂Y
```

where Ψ̃0 , Ψ̂0 are independent of R0 , u0 , S , and c3 (Y, T ). Note that Ψ̃0 = ψ̃1 and Ψ̂0 = ψ̂1 where ψ̃1 and ψ̂1 are from the linear problem (3.1.16). At O(l2 ) we have

```text
                                 ∂2 u2 ∂u0 ∂2 u0 ∂Ψ0 ′ ∂Ψ0 ∂u0 ∂Ψ0 ∂u0
                                        =        -        -       U +             -                                   (4.1.30)
                                  ∂z2       ∂T       ∂Y 2     ∂Y         ∂z ∂Y         ∂Y ∂z
```

with boundary conditions ∂u2

```text
                                                          = 0 on z = 0, -1.                                           (4.1.31)
                                                      ∂z
```

We find " z " z

```text
                                                  ∂Ψ0 ′                  ∂Ψ0 ∂u0
                                 u2 = -                U dz dz +                    dz dz
                                                   ∂Y                     ∂z ∂Y
                                    ∂u0 ∂2 u0 z2
                                 +(      -          ) + c8 (Y, T )z + c9 (Y, T ) = ũ2 + û2 c9 (Y, T )               (4.1.32)
                                    ∂T        ∂Y 2 2
```

where c9 (Y, T ) is an arbitrary function of Y and T . Here ũ2 and û2 are independent of c9 (Y, T ). Also at O(l2 ) we have

```text
                                   ∂2 Θ2 ∂Θ0         ∂2 Θ0 ∂Ψ0 ′ ∂Ψ0 ∂Θ0 ∂Ψ0 ∂Θ0
                                 τ 2 =           -τ         -       H +              -                               (4.1.33)
                                    ∂z      ∂T       ∂Y 2       ∂Y           ∂z ∂Y      ∂Y ∂z
```

with boundary conditions ∂Θ2

```text
                                                 + β1,0 Θ2 + β1,2 Θ0 = 0 on z = 0,                                   (4.1.34)
                                             ∂z
                                           ∂Θ2
                                                + β2,0 Θ2 + β2,2 Θ0 = 0 on z = -1.                                   (4.1.35)
                                            ∂z
```

We find " z " z 2 " z " z

```text
                               ∂Θ0 1                ∂ Θ0                    ∂Ψ0 H ′             1 ∂Ψ0 ∂Θ0
              Θ2 =                    dz dz -              dz dz -                  dz dz -               dz dz
                               ∂T τ                 ∂Y 2                     ∂Y τ               τ ∂Y ∂z
                          " z
                                 1 ∂Ψ0 ∂Θ0
                        +                    dz dz + c10 (Y, T )z + c11 (Y, T ).                                     (4.1.36)
                                 τ ∂z ∂Y
```

For case II we have Θ2 = Θ̂2 c11 (Y, T ) + Θ̃2 where Θ̂2 and Θ̃2 are independent of c11 (Y, T ). The boundary conditions lead to further equations which differ for the separate cases. For case I we find Z 0

```text
                                               ∂u0 ∂2 u0          ∂Ψ0 ′
                                                  -        =          U dz.                                     (4.1.37)
                                               ∂T    ∂Y 2      -1 ∂Y
```

<!-- Page 11 -->

If we now use (4.1.29) for Ψ0 we obtain an equation for u0 Z 0

```text
                                              ∂u0 ∂2 u0
                                                                                             !
                                                 -      1 - R0                     ψ̃1 U ′ dz = 0                   (4.1.38)
                                              ∂T   ∂Y 2                       -1
```

which on making use of equation (3.1.25) becomes

```text
                                                           ∂u0     ∂2 u0
                                                               + σ2 2 = 0.                                          (4.1.39)
                                                           ∂T      ∂Y
```

A periodic in y Fourier cosine solution to equation (4.1.39) is (Hayes & Phillips, 2017) ∞ X 2

```text
                                                   u0 (Y, T ) =          h p eσ2 p T cos pY                         (4.1.40)
                                                                   p=0
```

where h p are constant coefficients. For case II, we have two coupled nonlinear partial differential equations for u0 and c3 (Y, T ) as (4.1.37) and a further lengthy equation in the Appendix. In special cases such as βi,2m = 0 for m , 1 these partial differential equations are then linear and exact solutions can be found. Moreover when βi,2m = 0 for m , 1 there are two equations in terms of u0 and c3 (Y, T ) as Z 0

```text
                                      ∂u0 ∂2 u0                     ∂2 u0             ∂ 2 c3
                                          -        =        (-R   0       Ψ̃ 0 +  S          Ψ̂0 )U ′ dz,            (4.1.41)
                                      ∂T     ∂Y 2        -1         ∂Y 2              ∂Y 2
                                            Z 0
                           1 ∂c3 ∂2 c3                   ∂2 u0           ∂2 c3         H′
                                  -       =      (-R   0       Ψ̃ 0 +  S        Ψ̂ 0 )      dz - β1,2 c3 + β2,2 c3 . (4.1.42)
                           τ ∂T      ∂Y 2     -1         ∂Y 2             ∂Y 2          τ
```

We assume X∞ X∞

```text
                                        u0 =      f p (T ) cos pY, c3 =            g p (T ) cos pY.                  (4.1.43)
                                                p=0                                p=0
```

Substituting into the two coupled partial differential equations for u0 and c3 and equating like harmonics yields

```text
                                                  f˙p (T ) + a p f p (T ) + b p g p (T ) = 0,                       (4.1.44)
                                                  ġ p (T ) + c p f p (T ) + d p g p (T ) = 0                       (4.1.45)
```

where the constants a p , b p , c p , d p are given in the Appendix. If b p , 0 we find

- f˙p (T ) - a p f p (T )

```text
                                                      g p (T ) =                             ,                      (4.1.46)
                                                                              bp

                                      f¨p (T ) + (a p + d p ) f˙p (T ) + (d p a p - c p b p ) f p (T ) = 0.         (4.1.47)
```

The latter is a simple second order differential equation. We will omit the expressions for f p (T ), g p (T ). Note for this case that (a p + d p )/τ = b and (d p a p - c p b p )/τ = c when p = 1 where b and c appear in (3.1.28). At O(l3 ) we have

```text
                              ∂4 Ψ2          ∂3 Ψ0       ∂4 Ψ0           ∂u0         ∂u2    ∂Θ2
                                       =            - 2          - R2 D′     - R0 D′     +S
                               ∂z4          ∂T ∂z 2        2
                                                        ∂Y ∂z  2         ∂Y          ∂Y     ∂Y
                                                    3
                                              ∂Ψ0 ∂ Ψ0 ∂Ψ0 ∂ Ψ0    3
                                            -            +                                                          (4.1.48)
                                               ∂Y ∂z3        ∂z ∂Y∂z2
```

with boundary conditions ∂2 Ψ2

```text
                                                    = Ψ2 = 0 on z = 0, -1.                                          (4.1.49)
                                               ∂z2
```

The solution to this problem can be found to be " z " z 2 & z

```text
                                    ∂Ψ0                   ∂ Ψ0                         ∂u0
                     Ψ2 =                 dz dz -       2        dz dz -         D′ R2      dz dz dz dz
                                     ∂T                    ∂Y 2                        ∂Y
                                & z                                  & z
                                                ∂u2                          ∂Θ2
                              -           D′ R0     dz dz dz dz + S                dz dz dz dz
                                                ∂Y                            ∂Y
                                & z                                  & z
                                          ∂Ψ0 ∂3 Ψ0                          ∂Ψ0 ∂3 Ψ0
                              +                       dz  dz dz dz -                     dz dz dz dz
                                           ∂z ∂Y∂z2                           ∂Y ∂z3
                                          z3            z2
                              +c12 (Y, T ) + c13 (Y, T ) + c14 (Y, T )z + c15 (Y, T ).                              (4.1.50)
                                          6              2
```

<!-- Page 12 -->

For case I we have Ψ2 = Ψ̂2 - R2 Ψ̃2 + S Ψ̆2 where Ψ̂2 , Ψ̃2 , and Ψ̆2 are each independent of R2 and S . For case II we

```text
have Ψ2 = Ψ̆2 + Ψ̃2 ∂c9∂Y
                       (Y,T )
                              + Ψ̂2 ∂c11∂Y(Y,T ) where Ψ̂2 , Ψ̃2 , and Ψ̆2 are each independent of c9 (Y, T ) and c11 (Y, T ).
```

At O(l ) we have

```text
                                      ∂ 2 u4        ∂u2 ∂2 u2 ∂Ψ2 ′
                                                =        -      -    U
                                       ∂z2          ∂T     ∂Y 2   ∂Y
                                                      ∂Ψ2 ∂u0 ∂Ψ0 ∂u2 ∂Ψ2 ∂u0 ∂Ψ0 ∂u2
                                                    -          -       +       +                                      (4.1.51)
                                                       ∂Y ∂z     ∂Y ∂z   ∂z ∂Y   ∂z ∂Y
```

with boundary conditions as ∂u4

```text
                                                            + γ2 u0 = 0 on z = 0,                                     (4.1.52)
                                                        ∂z
                                                      ∂u4
                                                           + γ4 u0 = 0 on z = -1.                                     (4.1.53)
                                                       ∂z
```

We find " z " z 2 " z

```text
                                        ∂u2                   ∂ u2             ∂Ψ2 ′
                        u4 =                  dz dz -            2
                                                                   dz dz -          U dz dz
                                        ∂T                    ∂Y                ∂Y
                                    " z                        " z                  " z
                                           ∂Ψ0 ∂u2                  ∂Ψ2 ∂u0              ∂Ψ0 ∂u2
                                  -                  dz dz +                dz dz +              dz dz
                                            ∂Y ∂z                     ∂z ∂Y               ∂z ∂Y
                                  +c16 (Y, T )z + c17 (Y, T )                                                         (4.1.54)
```

where c17 (Y, T ) is arbitrary. Also at O(l4 ) we have

```text
                                      ∂2 Θ4         ∂Θ2     ∂2 Θ2 ∂Ψ2 ′
                                  τ             =        -τ      -    H
                                       ∂z2           ∂T     ∂Y 2   ∂Y
                                                      ∂Ψ2 ∂Θ0 ∂Ψ0 ∂Θ2 ∂Ψ2 ∂Θ0 ∂Ψ0 ∂Θ2
                                                    -          -        +       +                                     (4.1.55)
                                                       ∂Y ∂z     ∂Y ∂z    ∂z ∂Y   ∂z ∂Y
```

with boundary conditions as ∂Θ4

```text
                                               + β1,0 Θ4 + β1,2 Θ2 + β1,4 Θ0 = 0 on z = 0,                            (4.1.56)
                                            ∂z
                                          ∂Θ4
                                               + β2,0 Θ4 + β2,2 Θ2 + β2,4 Θ0 = 0 on z = -1.                           (4.1.57)
                                           ∂z
```

We find " z " z 2 " z

```text
                               1 ∂Θ2               ∂ Θ2                    ∂Ψ2 H ′
                   Θ4    =           dz dz -              dz dz  -                 dz dz
                               τ ∂T                 ∂Y 2                    ∂Y τ
                             " z                    " z                          " z
                                 1 ∂Ψ2 ∂Θ0                 1 ∂Ψ0 ∂Θ2                   1 ∂Ψ2 ∂Θ0
                           -                dz dz -                      dz dz +                 dz dz
                                 τ ∂Y ∂z                   τ ∂Y ∂z                     τ ∂z ∂Y
                             " z
                                 1 ∂Ψ0 ∂Θ2
                           +                dz dz + c18 (Y, T )z + c19 (Y, T ).                                       (4.1.58)
                                 τ ∂z ∂Y
```

The boundary conditions lead to equations for c9 (Y, T ) and c11 (Y, T ). For case I we have a single partial differential equation for c9 (Y, T ) appearing as

```text
                                              Z 0       Z 0 2          Z 0
                                               ∂u2           ∂ u2          ∂Ψ2 ′
                                         -         dz +         2
                                                                  dz +         U dz
                                            -1 ∂T         -1 ∂Y         -1 ∂Y
                                           Z 0               Z 0
                                               ∂Ψ0 ∂u2            ∂Ψ0 ∂u2
                                         +              dz -              dz - γ2 u0 + γ4 u0 = 0.                     (4.1.59)
                                            -1 ∂Y ∂z           -1 ∂z ∂Y
```

For case II we have two coupled partial differential equations in terms of c9 (Y, T ) and c11 (Y, T ) as (4.1.59) and a further very lengthy equation in the Appendix. At higher orders the complexity of the calculation becomes unwieldy. This calculation recovers Hayes & Phillips (2017) on setting S = 0. For time varying solutions it may be more convenient to use numerical methods such as those in §4.2, §4.3.

<!-- Page 13 -->

#### 4.1.2 Nonlinear perturbation solution algorithm

In light of the nonlinear perturbation solution above we let L X X ∞ L X X ∞ L X X ∞ 2k 2k

```text
          Ψ̄ =             Ψ2k,m sin(mY)l , ū =              u2k,m cos(mY)l , Θ̄ =                    Θ2k,m cos(mY)l2k       (4.1.60)
                 m=0 k=0                            m=0 k=0                                  m=0 k=0
```

with ∞ X ∞ X

```text
                                           R=         R2k l2k , γi = l4 γi , βi =         βi,2k l2k                           (4.1.61)
                                                k=0                                 k=0
where Ψ2k,m , u2k,m , and Θ2k,m are functions of z and T . We substitute (4.1.60), (4.1.61) into equations (4.1.3) to (4.1.7)
```

and discard harmonics in Y in the residuals that are of higher order than in the expansion of the solution in accordance with Theorem B in the Appendix. We then equate like harmonics in Y and like powers of l and then need to solve the resulting equations for Ψ2k,m , u2k,m , and Θ2k,m at each order in l. With the nonlinear perturbation solution we are particularly interested in the nonlinear steady states, for which we set ∂/∂T = 0. In this case, arbitrary constants of integration will appear in the nonlinear perturbation solution. We choose them so that ui, j |z=0 = δi,0 δ j,1 . Note that while this choice is dissimilar to that in the linear perturbation solutions it is similar to that in the linear numerical solutions. We found that the nonlinear steady states appear to require restrictions on the boundary conditions at O(l6 ) such as γ3 = γ4 = 0. This may be related to observations where LC tend to curl up near the bottom of the mixed layer. This calculation recovers Hayes & Phillips (2017) on setting S = 0.

### 4.2 Nonlinear power series method

Here we look for solutions of the form L 4+M X X L 2+M X X L 2+M X X m m

```text
                   ψ=              am,k sin(kly)z , u =             bm,k cos(kly)z , θ =                   cm,k cos(kly)zm     (4.2.1)
                        k=0 m=0                           k=0 m=0                                k=0 m=0
```

where the coefficients am,k , bm,k , and cm,k are unknown functions of t. Here ψ is a Fourier sine series in y while u and θ are both Fourier cosine series in y; each are Maclaurin series in z. Substituting into the governing equations and equating the appropriate like coefficients in accordance with Theorem A and Theorem B leads to a system of nonlinear ordinary differential equations for am,k , bm,k , and cm,k which can be numerically solved for by using methods such as the Runge-Kutta method. In the case for which dam,k /dt = dbm,k /dt = dcm,k /dt = 0 this leads to a system of algebraic equations. These algebraic equations are treated much the same as in the linear power series method. We here set b0,1 = 1 for consistency with the nonlinear perturbation solutions.

### 4.3 Nonlinear Galerkin method

In this method we look for solutions of the form X M L X L M-2 X X L M-2 X X

```text
                  ψ=              am,k Pm sin(kly), u =             bm,k Pm cos(kly), θ =                  cm,k Pm cos(kly)    (4.3.1)
                        k=0 m=0                           k=0 m=0                                k=0 m=0
```

where am,k , bm,k , and cm,k are unknown functions of t. Here ψ is a Fourier sine series in y while u and θ are Fourier cosine series’ in y. Different are the basis functions. Pm (z) are shifted Legendre basis functions on z ∈ [-1, 0]. We substitute these expansions into the CL2 equations, discard the higher order harmonics, and collect like trigonometrical terms in accordance with Theorem B to obtain a set of equations in z and t whose residuals we call r1,i (z, t), r2,i (z, t), and r3,i (z, t). In the Galerkin method we require Z 0 Z 0 Z 0

```text
                                          r1,i P j dz =     r2,i P j dz =     r3,i P j dz = 0                    (4.3.2)
                                           -1                 -1               -1
```

for i = 0, 1, . . . , L, and j = 0, 1, . . . , M - 4. We obtain the further equations required to close the system by substitution of (4.3.1) into the boundary conditions. This results in a system of nonlinear ordinary differential equations which can be solved numerically by using methods such as the Runge-Kutta method. For the case of nonlinear steady states they reduce to a set of algebraic equations. Once again, these algebraic equations are treated much the same as in the linear power series method. Herein we choose b0,1 so that the coefficient of cos(ly) in u|z=0 is unity for consistency with the nonlinear perturbation solutions and nonlinear power series solutions.

## 5 Results

In this section we are interested in how the parameters and nonlinearities affect the CL2 instability to LC over a restricted parameter range. For case I we let βi,2m = 0 for m , 0 and for case II we let βi,2m = 0 for m , 1. Herein ǫ = 0 is for the linear case and ǫ = 1 is for its nonlinear counterpart. We here choose L = 1 and 15 ⩽ M ⩽ 20.

<!-- Page 14 -->

### 5.1 Growth rate

We consider first the growth rate σ. For case I we find that σ is real valued and for case II σ is complex valued where we see that in accord with (3.1.28) there are two solutions. When Re σ < 0 the motion is stable and when Re σ > 0 the motion is unstable. When σ = 0 there is neutral instability. The instability is oscillatory when Im σ , 0. The growth rate σ from the linear perturbation solution to O(l4 ) for case I where D′ , U ′ , H ′ are constants is

1 691((β2 - 2077 2077 5544 ′ ′ ′ 691 )β1 + 691 β2 - 691 )RD U S H 4

```text
            σ=            (-530R2 D′2 U ′2 - 67320RD′ U ′ -                                              )l
                79833600                                                  ((β2 - 1)β1 + β2 )τ
                   RD′ U ′ 2
            +(-1 +        )l + γ4 - γ2                                                                      (5.1.1)
                    120
```

and the growth rate σ from the linear perturbation solution to O(l2 ) for case II where D′ , U ′ , H ′ are constants is

## 1 RD′ U ′ S H′

```text
               σ= (          - (β1,2 - β2,2 + 1)τ -       -1                                                          (5.1.2)
                r 2 120                             120
                   RD′ U ′ 2    RD′ U ′ S H ′                                S H′
               ± (        ) -2          (      - τ(β1,2 - β2,2 + 1) + 1) + (      + τ(β1,2 - β2,2 + 1) - 1)2 )l2 .
                    120           120 120                                    120
```

There can be uncertainty in deciding when case I or case II is appropriate. What happens is either the case I result or the case II result will converge or both case I and case II result will converge each for separate parts of the domain of discourse, and the appropriate case is that which converges. This is the competition between case I and case II as mentioned in Cox & Leibovich (1993). A good indication of whether the instability is case I or case II is whenever σ4 > σ2 for case I then σ for case I is likely to diverge and so the appropriate instability is then case II. Plots of σ which do illustrate this competition are shown in Figures 1, 2. Plots of Re σ vs R and Im σ vs R for case II are shown

Figure 1: Plots of linear growth rate (left) Re σ vs R and (right) Im σ vs R for β1 = 1/100. Here D′ = U ′ =

```text
H ′ = 1, S = 100, l = 1/10, γ1 = 1/20000, γ2 = 1/10000, γ3 = γ4 = 0, β2 = 0, and τ = 1/10.
```

Figure 2: Plots of linear growth rate σ vs R for (left) β1 = 1/10 and (right) β1 = 1. Here D′ = U ′ = H ′ = 1,

```text
S = 100, l = 1/10, γ1 = 1/20000, γ2 = 1/10000, γ3 = γ4 = 0, β2 = 0, and τ = 1/10.
```

in Figure 1 and plots of σ vs R for case I and case II are shown in Figure 2. In these plots we see that the instability changes from case II to case I as β1 increases. We also see that the fluid motion switches from stable to unstable as R increases and so here increasing R is destabilising. It is then quite obvious from (5.1.1) and (5.1.2) how the parameters would affect σ in the small l limit where the expressions are valid. For example, increasing D′ or U ′ is destabilising whenever increasing R is destabilising, and increasing H ′ is stabilising whenever increasing S is stabilising. For the boundary conditions of Cox & Leibovich (1993) we see in case I that increasing R or τ is destabilising and increasing S is stabilising. We also see in case I that increasing γ2 - γ4 is stabilising. For case II with the boundary conditions of Cox & Leibovich (1993) and on assuming σ remains complex we see that increasing R, D′ , or U ′ is destabilising and increasing S , H ′ , or τ is stabilising.

<!-- Page 15 -->

### 5.2 Neutral instability

For case II, we see from (3.1.28) that neutral instability for which σ = 0 is seldom possible. Linear neutral curves and nonlinear steady states do exist for case I. From the case I linear and nonlinear perturbation solution for neutral instability at O(l2 ) we have

```text
                                                   R0 = R 0            .                                        (5.2.1)
                                                             ψ̃ U ′ dz
                                                          -1 1
```

From the case I linear perturbation solution for neutral instability at O(l4 ) we have R0 (ψ̂3 + S ψ̆3 )U ′ dz + u0 (γ2 - γ4 )

```text
                          R2 = -1           R0                          = R̃2 + R̂2 (γ2 - γ4 ) + R̆2 S .             (5.2.2)
                                                 ψ̃  U ′ dz
                                              -1   3
```

From the case I nonlinear perturbation solution for neutral instability at O(l4 ) we have

```text
                          Z 0                               Z 0 2          Z 0
                                ′ ∂Ψ̆2                          ∂ u2             ∂Ψ̂2 ′
              R2 = (S         U        dz - γ2 u0 + γ4 u0 +         2
                                                                      dz +           U dz
                           -1      ∂Y                        -1 ∂Y           -1 ∂Y
                         Z 0                Z 0                 Z 0
                             ∂Ψ0 ∂u2              ∂Ψ0 ∂u2              ∂Ψ̃2
                       +               dz -                dz)/     U′      dz = R̃2 + R̂2 (γ2 - γ4 ) + R̆2 S .      (5.2.3)
                          -1 ∂Y ∂z            -1 ∂z ∂Y           -1     ∂Y
```

In these equations R̃2 , R̂2 , and R̆2 are each independent of γi and S . Also note that nonlinear R2 is here projected onto a mode in Y. We here choose L = 1 in the nonlinear expansions. For both the case I linear and nonlinear problems the expression for R appears as R̂2

```text
                                          R = R0 + (R̃2 + R̆2 S )l2 + 2 (γ2 - γ4 ) + . . . .                          (5.2.4)
                                                                      l
```

The neutral curve from the perturbation solution to O(l4 ) for case I where D′ , U ′ , H ′ are constants is 2077 2077 5544 ′ 2 5455 l2 691 ((β2 - 691 )β1 + 691 β2 - 691 )S H l 1550 ǫl2 120 120(γ2 - γ4 )

```text
      R=        ′  ′
                     +                ′  ′
                                                                  +             + ′ ′+               .               (5.2.5)
           231 D U     5544        τD U ((β2 - 1)β1 + β2 )               ′
                                                                     21 D U  ′3  DU      D′ U ′ l 2
```

Here it is evident that 5455 1 1550 ǫ

```text
                                                R̃2 =            +             ,                                     (5.2.6)
                                                      231 D′ U ′    21 D′ U ′3
                                                           2077      2077     5544   ′
                                                691 ((β2 - 691 )β1 + 691 β2 - 691 )H
                                        R̆2 =                                       ,                               (5.2.7)
                                                5544 τD′ U ′ ((β2 - 1)β1 + β2 )
                                                         120
                                                  R̂2 = ′ ′ = R0 .                                                  (5.2.8)
                                                        DU
```

We plot linear neutral curves and nonlinear steady states as R vs l in the small l limit in Figure 3 (left). For this case,

Figure 3: (left) Plots of neutral curve R vs l, linear (top) and nonlinear (bottom). (right) Plot of d = RL - RNL

```text
vs l. Here D′ = U ′ = H ′ = 1, S = 100, γ1 = 1/20000, γ2 = 1/10000, γ3 = γ4 = 0, β1 = 1, β2 = 0, and
τ = 1/10.
```

since the fluid motion switches from stable to unstable as σ passes through zero with increasing R, any point above the neutral curve is unstable, while any point below the neutral curve is stable. In Figure 3 (right) we see that nonlinearities are small when l ≪ 1 similarly to as shown in Hayes & Phillips (2017) for the case S = 0. In the small l limit we see that nonlinearities have a stabilising effect. It is also quite obvious from (5.2.5) how the parameters and nonlinearities would affect neutral instability in the small l limit where this expression is valid. From (5.2.5) we see for the boundary conditions of Cox & Leibovich (1993) that increasing S is stabilising. This is consistent with Langmuir (1938) in that temperature is thought to be secondary to the formation of LC. Increasing H ′ or decreasing τ has a similar effect as increasing S . Also we see from (5.2.5) that nonlinearities are stabilising. When using the boundary conditions of Cox & Leibovich (1993) we find for the nonlinear problem that increasing D′ or U ′ is destabilising and increasing U ′ is more effective in destabilising the flow than increasing D′ . The effect of increasing γ2 - γ4 is stabilising.

<!-- Page 16 -->

### 5.3 Figure 3 of Cox & Leibovich (1993) revisited

We are interested in the nonlinear counterpart to Figure 3 of Cox & Leibovich (1993). Figure 4 (left) is Figure 3 of Cox & Leibovich (1993). The flat parts of Figure 4 (left) represent the stability margin for oscillatory convection and the parabolic parts of Figure 4 (left) represent the stability margin to steady convection (Cox & Leibovich, 1993). The nonlinear version of Figure 3 of Cox & Leibovich (1993) is plotted within Figure 4 (right) for the corresponding nonlinear steady states only. In Figure 4, lc increases with increasing β1 . To obtain the flat parts of Figure 4 (left)

Figure 4: (left) Figure 3 of Cox & Leibovich (1993). (right) nonlinear steady states R vs l. Here

```text
D′ = U ′ = H ′ = 1, S = 100, γ1 = 1/40000, γ2 = -γ3 = -γ4 = 1/20000, β1 = -β2 ∈
{1/20000, 1/2000, 1/1000, 1/500, 1/200, 1/20, 1/2}, and τ = 10/67.
```

where Re σ = 0 with the case II linear perturbation solution, we solve b = 0 for R0 and then Re σ2 = 0 providing c ⩾ 0 where b and c appear in (3.1.28). Then Re σ2 j = 0 can be solved for R2 j-2 for j > 1 on assuming c ⩾ 0. Figure 4 (right) was obtained by using the Galerkin method as outlined in §4. The Galerkin method used for finding the steady states assumes ∂/∂t = 0 and so Figure 4 (right) represents the nonlinear counterpart of the stability margin to steady convection in Figure 4 (left) only up to a certain l value for each curve. Moreover, if we plot Re σ vs R and Im σ vs R for values consistent with Figure 4, then for l before the bifurcation point in Figure 4 (left) the growth rate is case I like that of Figure 2 (right) and for l after the bifurcation point in Figure 4 (left) the growth rate is case II like that of Figure 1. If we solve for σ = 0 in case II our Galerkin method will find where only one of the branches of σ is zero, but both branches must be included for the solution to be real valued. A similar idea must also apply to the nonlinear case because in the nonlinear perturbation solution it is found that the solution is linear up until O(l2 ). Finding a nonlinear counterpart to σ at higher orders appears to be quite difficult and is omitted. In Figure 4 (right) nonlinearities are small and stabilising in the small l limit. When Figure 4 (left) and (right) are overlayed the neutral curves for the linear case appear indistinguishable to the corresponding nonlinear steady states about their minimum turning points.

### 5.4 Onset

Onset occurs at the minimum point on the neutral curve R = R(l), which we denote by (lc , Rc ) where lc and Rc are called the critical wavenumber and critical Rayleigh number respectively. From the perturbation solutions onset is found by solving dR 4 dl = 0. We find from the O(l ) perturbation solutions for case I that

!1 (γ2 - γ4 )R̂2 4

```text
                                                         lc =                                                           (5.4.1)
                                                               R̃2 + R̆2 S
```

and thus that

```text
                            Rc = R0 + 2((γ2 - γ4 )R̂2 )1/2 (R̃2 + R̆2 S )1/2 = R0 + 2lc2 (R̃2 + R̆2 S ).                (5.4.2)
```

The critical wavenumber lc from the perturbation solution to O(l4 ) for case I where D′ , U ′ , H ′ are constants is   41 (γ2 - γ4 ) D120    ′U′ 

```text
                              lc =                                       2077      2077     5544    ′
                                                                                                                     (5.4.3)
                                              1
                                          5455      1550 ǫ             ((β -      )β +      β -      )H
                                         231 D′ U ′ + 21 D′ U ′3   + 691 2 691 1 691 2 691 S 
                                                                      5544       τD′ U ′ ((β2 -1)β1 +β2 )
```

and the corresponding critical Rayleigh number Rc from the perturbation solution to O(l4 ) for case I where D′ , U ′ , H ′ are constants is 2077 2077 5544 ′ 120 120 1/2 5455 1 1550 ǫ 691 ((β2 - 691 )β1 + 691 β2 - 691 )H 1/2

```text
 Rc =          + 2((γ 2 - γ4 )        ) (            +            +                                        S ) . (5.4.4)
        D′ U ′                 D′ U ′     231 D′ U ′    21 D′ U ′3 5544      τD′ U ′ ((β2 - 1)β1 + β2 )
```

It is here quite obvious from (5.4.3) and (5.4.4) how the parameters and nonlinearities would affect lc and Rc in the small l limit where these expression are valid. From (5.4.3) and (5.4.4) we see for the boundary conditions of Cox & Leibovich (1993) that increasing S reduces lc and increases Rc . Here we also see that increasing H ′ or decreasing τ has a similar effect as increasing S . Also we see that nonlinearities reduce the value of lc and increase Rc . Increasing

<!-- Page 17 -->

D′ has no effect on lc . Increasing U ′ only has an effect on lc in the presence of nonlinearities. When using the boundary conditions of Cox & Leibovich (1993) we find for the nonlinear problem that increasing U ′ increases lc , and increasing D′ or U ′ decreases Rc where increasing U ′ is more effective in decreasing Rc than increasing D′ . The effect of increasing γ2 - γ4 is to increase the value of lc and increase the value of Rc , and we see that γ2 = γ4 = 0 leads to unphysical results. Another peculiarity is that depending on the choice of β1 , β2 , and different from Cox & Leibovich (1993), we see that there can be a singularity of lc when S increases. In Figure 5 are plots of lc vs S and Rc vs S for the linear and nonlinear problems. In Figure 5 the effect of increasing S is to lower the value of lc and increase the

Figure 5: (left) Plots of lc vs S , linear (top) and nonlinear (bottom). (right) Plots of Rc vs S , linear (bottom)

```text
and nonlinear (top). Here D′ = U ′ = H ′ = 1, γ1 = 1/20000, γ2 = 1/10000, γ3 = γ4 = 0, β1 = 1, β2 = 0, and
τ = 1/10.
```

value of Rc . In Figure 6 are plots of the ratio of linear to nonlinear critical wavenumber κ = lc,linear /lc,nonlinear vs S and plots of the ratio of linear to nonlinear critical Rayleigh number ρ = Rc,linear /Rc,nonlinear vs S for parameters consistent with Figure 5. In Figure 6 we see that the nonlinearities appear to diminish as S increases. Figure 7 shows how lc and

```text
Figure 6: (left) Plots of κ vs S . (right) Plots of ρ vs S . Here D′ = U ′ = H ′ = 1, γ1 = 1/20000, γ2 = 1/10000,
γ3 = γ4 = 0, β1 = 1, β2 = 0, and τ = 1/10.
```

Rc varies with β1 for both the linear and nonlinear cases with S = 100 and other parameters consistent with Figure 5. In Figure 7 we see that lc increases with increasing β1 and Rc decreases with increasing β1 . In Figure 8 are plots

Figure 7: (left) Plots of lc vs β1 , linear (top) and nonlinear (bottom). (right) Plots of Rc vs β1 , linear (bottom)

```text
and nonlinear (top). Here D′ = U ′ = H ′ = 1, γ1 = 1/20000, γ2 = 1/10000, γ3 = γ4 = 0, S = 100, β2 = 0,
and τ = 1/10.
```

of the ratio of linear to nonlinear critical wavenumber κ = lc,linear /lc,nonlinear vs β1 and plots of the ratio of linear to nonlinear critical Rayleigh number ρ = Rc,linear /Rc,nonlinear vs β1 for S ∈ {0, 100, 200} and other parameters consistent with Figure 5. In Figure 8 the κ curves decrease for increasing S and the ρ curves increase for increasing S . We see for

```text
S = 0 that κ and ρ are independent of β1 as expected. For S = 100 and S = 200 we see that κ increases with increasing
```

β1 and ρ decreases with increasing β1 . Also, the nonlinearities appear to be small for small β1 . For S = 0 we find that

<!-- Page 18 -->

```text
Figure 8: (left) Plots of κ vs β1 . (right) Plots of ρ vs β1 . Here D′ = U ′ = H ′ = 1, γ1 = 1/20000,
γ2 = 1/10000, γ3 = γ4 = 0, S ∈ {0, 100, 200}, β2 = 0, and τ = 1/10.
```

κ ≈ 1.425 at O(l4 ) which is consistent with the value reported in Hayes & Phillips (2017). As shown in Figure 8 this is κ ≈ 1.433 at O(l6 ). In Figures 5 to 8 we see that nonlinearities reduce the value of lc and increase Rc .

## 6 Discussion

The methods used in this paper are very useful for the LC problem. The perturbation method is particularly useful in that the effect of altering parameters and of nonlinearities is evident in the small l limit by inspecting the simple expressions found from the perturbation solutions. Note that when γi = O(1) our perturbation solutions would require a more direct perturbation expansion. For example, in the linear perturbation solution σ would then need to have an O(1) term σ0 . This then leads to a messy calculation, especially for its nonlinear counterpart, with many separate subcases. The preferable strategy may then be to use the numerical methods such as the nonlinear power series and nonlinear Galerkin methods presented in this paper. I have also constructed animations of LC varying with time. In the nonlinear realm there is flexibility for animations of LC to show LC spacing changing with time due to the fact that the number of modes in y can increase as time increases. This is to be explored in further work on LC.

## 7 Appendix

### 7.1 Linear perturbation solution details

From ψ1 " z " z

```text
                               c4 = -R0 u0           D′ dz dz|z=-1 + S            θ0 dz dz|z=-1 + c5 ,
                                                  " z                        " z
                                                          ′
                                    c5 = R0 u0          D dz dz|z=0 - S               θ0 dz dz|z=0 ,
                             & z                                    & z
                                                                                      1       1
               c6 = -R0 u0           D′ dz dz dz dz|z=-1 + S    θ0 dz dz dz dz|z=-1 - c4 + c5 + c7 ,
                                                                                      6       2
                                     & z                        & z
                          c7 = R0 u0     D′ dz dz dz dz|z=0 - S          θ0 dz dz dz dz|z=0 .
```

From u2 Z z

```text
                                                     c8 = -         U ′ ψ1 dz|z=0 ,
                                            Z 0" z
                                                                             u0           1
                                   c9 = -               U ′ ψ1 dz dz dz -       (σ2 + 1) + c8 .
                                             -1                              6            2
```

From θ2 for case I " z Z z H′ H′

```text
                               c10 = -β1,0 (            ψ1 dz dz|z=0 + c11 ) -                 ψ1 dz|z=0 ,
                                                     τ                                      τ
```

Z 0 " z ′ H′ H

```text
                c11 = (            ψ1 dz + (β1,0 - β1,0 β2,0 )        ψ1 dz dz|z=0
                            -1 τ                                    τ
                                " z ′                           Z z ′
                                      H                            H
                          -β2,0          ψ1 dz dz|z=-1 - β2,0         ψ1 dz|z=0 )/(-β1,0 + β1,0 β2,0 + β2,0 ).
                                       τ                            τ
```

<!-- Page 19 -->

From θ2 for case II " z " z " z ′ θ 0 σ2 H

```text
                   c10   = -β1,0 (              dz dz|z=0 +       θ0 dz dz|z=0 +          ψ1 dz dz|z=0 + c11 )
                                            τ                                          τ
                             Z z                   Z z             Z z ′
                                   θ0 σ2                                H
                           -             dz|z=0 -      θ0 dz|z=0 -         ψ1 dz|z=0 - β1,2 θ0 |z=0 .
                                     τ                                   τ
```

From the matrix M for case II Z 0

```text
                                                    M1,1 = σ2 + 1 - R0              ψ̃1 U ′ dz,
                                                                               -1
                                                                     Z 0
                                                         M1,2 = S           ψ̂1 U ′ dz,
                                                                       -1
```

Z z " z Z z H′ H′ H′

```text
          M2,1    = R0     ψ̃1    dz|z=-1 - β1,0 R0      ψ̃1     dz dz|z=0 - R0      ψ̃1    dz|z=0
                               τ                              τ                          τ
                             " z                                   "   z                            Z z
                                      H′                                     H′                             H′
                    +β2,0 R0      ψ̃1    dz dz|z=-1 + β2,0 β1,0 R0       ψ̃1    dz dz|z=0 + β2,0 R0     ψ̃1    dz|z=0 ,
                                       τ                                     τ                              τ
```

Z z " z Z z θ̂0 θ̂0

```text
                     M2,2      = -σ2 (           dz|z=-1 + β2,0                dz dz|z=-1 ) -      θ̂0 dz|z=-1
                                             τ                             τ
                                    Z z                                "     z                         Z z
                                              H′                                   H′                          H′
                                 -S      ψ̂1       dz|z=-1 + β1,0 S            ψ̂1     dz dz|z=0 + S       ψ̂1    dz|z=0
                                               τ                                    τ                          τ
                                        " z                               " z
                                                                                       H′
                                 -β2,0         θ̂0 dz dz|z=-1 - β2,0 S             ψ̂1    dz dz|z=-1
                                                                                       τ
                                               " z                                    Z z
                                                          H′                                  H′
                                 -β2,0 β1,0 S         ψ̂1    dz dz|z=0 - β2,0 S           ψ̂1    dz|z=0
                                                          τ                                   τ
                                 +β1,2 (1 - β2,0 )θˆ0 |z=0 - β2,2 θˆ0 |z=-1 .
```

From the quadratic equation for σ2 for case II Z z " z θ̂0 θ̂0

```text
                                   a=-             dz|z=-1 - β2,0         dz dz|z=-1 ,
                                                τ                      τ
                         Z z                    Z z                            " z                            Z z
                                                        H′                                  H′                         H′
             b = -           θ̂0 dz|z=-1 - S        ψ̂1     dz|z=-1 + β1,0 S           ψ̂1      dz dz|z=0 + S      ψ̂1    dz|z=0
                                                        τ                                    τ                         τ
                            " z                              " z                                        " z
                                                                        H′                                        H′
                      -β2,0         θ̂0 dz dz|z=-1 - β2,0 S        ψ̂1       dz dz|z=-1 - β1,0 β2,0 S         ψ̂1     dz dz|z=0
                                                                         τ                                        τ
                               Z z
                                         H′
                      -β2,0 S       ψ̂1      dz|z=0 + β1,2 (1 - β2,0 )θˆ0 |z=0 - β2,2 θˆ0 |z=-1
                                          τ
                                  Z 0              Z z                        " z
                                             ′          θ̂0                         θ̂0
                      -(1 - R0          ψ̃1 U dz)(          dz|z=-1 + β2,0               dz dz|z=-1 ),
                                   -1                    τ                           τ
```

Z 0 Z z Z z " z ′ H′ H′

```text
          c = (1 - R0       ψ̃1 U dz)(-         θ̂0 dz|z=-1 - S        ψ̂1     dz|z=-1 + β1,0 S         ψ̂1     dz dz|z=0
                                                                            τ                                τ
                  Z z -1 ′                     " z                             " z
                           H                                                              H′
              +S       ψ̂1     dz|z=0 - β2,0          θ̂0 dz dz|z=-1 - β2,0 S         ψ̂1    dz dz|z=-1
                           τ                                                              τ
                           " z                                 Z z
                                       H′                               H′
              -β2,0 β1,0 S         ψ̂1    dz dz|z=0 - β2,0 S        ψ̂1     dz|z=0 + β1,2 (1 - β2,0 )θˆ0 |z=0 - β2,2 θˆ0 |z=-1 ),
                                       τ                                 τ
                  Z 0                  Z z                              " z                           Z z
                                               H′                                 H′                            H′
              -S       ψ̂1 U ′ dz(R0       ψ̃1      dz|z=-1 - β1,0 R0         ψ̃1     dz dz|z=0 - R0        ψ̃1    dz|z=0
                                               τ                                   τ                            τ
                    -1
                        " z                                        " z                                 Z z
                                   H′                                         H′                                H′
              +β2,0 R0         ψ̃1     dz dz|z=-1 + β2,0 β1,0 R0          ψ̃1     dz dz|z=0 + β2,0 R0       ψ̃1    dz|z=0 ).
                                    τ                                         τ                                 τ
```

From ψ3 " z " z ′

```text
                            c12 = -           D (R2 u0 + u2 R0 ) dz dz|z=-1 + S               θ2 dz dz|z=-1 + c13 ,
```

<!-- Page 20 -->

" z " z

```text
                                 c13 =         D′ (R2 u0 + u2 R0 ) dz dz|z=0 - S            θ2 dz dz|z=0 ,
```

" z & z

```text
                      c14 =        ψ1 (2 + σ2 ) dz dz|z=-1 -        D′ (R2 u0 + u2 R0 ) dz dz dz dz|z=-1
                                   & z
                                                               1       1
                                +S        θ2 dz dz dz dz|z=-1 - c12 + c13 + c15 ,
                                                               6       2
```

& z

```text
                             c15 =            D′ (R2 u0 + u2 R0 ) dz dz dz dz|z=0
                                         " z                             & z
                                       -     ψ1 (2 + σ2 ) dz dz|z=0 - S           θ2 dz dz dz dz|z=0 .
```

From u4 Z z Z z

```text
                                    c16 = -γ2 u0 -         u2 (1 + σ2 ) dz|z=0 -         U ′ ψ3 dz|z=0 ,
```

Z 0" z

```text
                                      c17 = -            u2 (1 + σ2 ) dz dz dz
                                                    -1
                                                                 Z 0" z
                                                 1     1
                                                + c16 - u0 σ4 -             U ′ ψ3 dz dz dz.
                                                 2     6          -1
```

From θ4 for case I " z ′ " z H σ2

```text
                           c18 = -β1,0 (         ψ3 dz dz|z=0 + (     + 1)       θ2 dz dz|z=0 + c19 )
                                              τ                    τ
                                   Z z ′                          Z z
                                         H               σ2
                                 -         ψ3 dz|z=0 - (    + 1)      θ2 dz|z=0 - β1,2 θ2 |z=0 ,
                                         τ               τ
```

Z z Z z H′ σ2

```text
     c19 = (-                ψ3 dz|z=-1 - (     + 1)      θ2 dz|z=-1
                          τ                  τ
                           " z ′                                " z
                                 H                    σ2
                 +(-β1,0 (           ψ3 dz dz|z=0 + (     + 1)       θ2 dz dz|z=0 )
                                  τ                    τ
                   Z z ′                           Z z
                         H                σ2
                 -          ψ3 dz|z=0 - (     + 1)      θ2 dz|z=0 - β1,2 θ2 |z=0 )(β2,0 - 1)
                         τ                τ
                         " z ′                                " z
                              H                     σ2
                 -β2,0 (          ψ3 dz dz|z=-1 + (      + 1)       θ2 dz dz|z=-1 ) - β2,2 θ2 |z=-1 )/(-β1,0 + β1,0 β2,0 + β2,0 ).
                               τ                     τ
```

From θ4 for case II " z " z " z H′ σ2 σ4

```text
           c18     = -β1,0 (         ψ3 dz dz|z=0 + (     + 1)       θ2 dz dz|z=0 +       θ0     dz dz|z=0 + c19 )
                                  τ                    τ                                     τ
                       Z z ′                          Z z             Z z
                             H               σ2                               σ4
                     -         ψ3 dz|z=0 - (    + 1)      θ2 dz|z=0 -      θ0    dz|z=0 - β1,2 θ2 |z=0 - β1,4 θ0 |z=0 .
                             τ               τ                                τ
```

From the matrix N for case II

```text
                                                              N1,1 = u0 ,
                                                                Z 0
                                                         N1,2 =     ψ̆3 U ′ dz,
                                                                  -1
```

Z z " z Z z θ0 θ0 θ0

```text
                    N2,1   =        dz|z=-1 - β1,0          dz dz|z=0 -         dz|z=0
                                  τ                      τ                   τ
                                   " z                             " z                      Z z
                                         θ0                             θ0                      θ0
                             +β2,0          dz dz|z=-1 + β2,0 β1,0         dz dz|z=0 + β2,0        dz|z=0 ,
                                         τ                              τ                       τ
```

<!-- Page 21 -->

Z z Z z Z z " z σ2 H′ σ2

```text
   N2,2 =         θ̂2     dz|z=-1 +         θ̂2 dz|z=-1 +         ψ̆3     dz|z=-1 - β1,0          θ̂2     dz dz|z=0
                       τ                                              τ                                τ
                      " z                            " z                         Z z                     Z z               Z z
                                                                H′                         σ2                                      H′
               -β1,0        θ̂2 dz dz|z=0 - β1,0            ψ̆3     dz dz|z=0 -        θ̂2    dz|z=0 -        θ̂2 dz|z=0 -     ψ̆3    dz|z=0
                                                                τ                          τ                                       τ
                      " z                                 " z                             " z
                                σ2                                                                  H′
               +β2,0        θ̂2      dz dz|z=-1 + β2,0           θ̂2 dz dz|z=-1 + β2,0          ψ̆3      dz dz|z=-1
                                 τ                                                                    τ
                          " z                                     " z                                " z
                                     σ2                                                                         H′
               +β1,0 β2,0        θ̂2     dz dz|z=0 + β1,0 β2,0           θ̂2 dz dz|z=0 + β1,0 β2,0         ψ̆3     dz dz|z=0
                                      τ                                                                          τ
                      Z z                         Z z                       Z z
                              σ2                                                    H′
               +β2,0      θ̂2     dz|z=0 + β2,0         θ̂2 dz|z=0 + β2,0       ψ̆3      dz|z=0
                               τ                                                     τ
               -β1,2 (1 - β2,0 )θˆ2 |z=0 + β2,2 θˆ2 |z=-1 .
```

From the vector q for case II Z 0 Z 0

```text
                                      q1 = -           ψ̂3 U ′ dz + R2         ψ̃3 U ′ dz + γ4 u0 - γ2 u0 ,
                                                  -1                      -1
```

Z z Z z Z z Z z σ2 H′ H′

```text
      q2 = -          θ̃2    dz|z=-1 -         θ̃2 dz|z=-1 -        ψ̂3    dz|z=-1 + R2         ψ̃3     dz|z=-1
                          τ                                             τ                            τ
                       " z                         " z                     " z                              " z
                                 σ2                                                   H′                            H′
               +β1,0 (       θ̃2      dz dz|z=0 +        θ̃2 dz dz|z=0 +          ψ̂3     dz dz|z=0 - R2        ψ̃3     dz dz|z=0 )
                                  τ                                                    τ                             τ
                 Z z                     Z z               Z z                            z
                                                                     H′                         H′
                                                                                       Z
                          σ2
               +      θ̃2    dz|z=0 +        θ̃2 dz|z=0 +        ψ̂3     dz|z=0 - R2        ψ̃3     dz|z=0
                          τ                                           τ                         τ
                      " z                            " z                      " z                              " z
                                 σ2                                                      H′                             H′
               -β2,0 (       θ̃2      dz dz|z=-1 +         θ̃2 dz dz|z=-1 +          ψ̂3     dz dz|z=-1 - R2        ψ̃3     dz dz|z=-1
                                  τ                                                      τ                               τ
                          " z                         " z                      " z                            " z
                                     σ2                                                  H′                            H′
               -(-β1,0 (         θ̃2    dz dz|z=0 +          θ̃2 dz dz|z=0 +         ψ̂3      dz dz|z=0 - R2       ψ̃3     dz dz|z=0 )
                                     τ                                                   τ                              τ
                 Z z                     Z z               Z z                            z
                                                                     H′                         H′
                                                                                       Z
                          σ2
               -      θ̃2    dz|z=0 -        θ̃2 dz|z=0 -        ψ̂3     dz|z=0 + R2        ψ̃3     dz|z=0 ))
                          τ                                           τ                         τ
               +β1,2 (1 - β2,0 )θ˜2 |z=0 + β1,4 (1 - β2,0 )θ0 |z=0 - β2,2 θ˜2 |z=-1 - β2,4 θ0 |z=-1 .
```

### 7.2 Linear perturbation solution algorithm details

From u2 j in the linear perturbation solution algorithm Z 0 Z z j-1 Z zX ′

```text
              c0, j = -γ2 u2 j-4 |z=0 -         U ψ2 j-1 dz|z=0 -         u2 j-2 dz|z=0 -               u2 j-(2m+2) σ2m+2 dz|z=0 ,
                                           -1                                                     m=0
```

Z 0 " zX j-1 " z " z c0, j

```text
        c1, j = -     [      u2 j-(2m+2) σ2m+2 dz dz +     u2 j-2 dz dz +     U ′ ψ2 j-1 dz dz] dz + δ0, j u0 +       .
                   -1    m=0
                                                                                                                 2
```

From θ2 j in the linear perturbation solution algorithm for case I " z " z ′ H

```text
              c2, j = -β1,0 (        θ2 j-2 dz dz|z=0 +           ψ2 j-1 dz dz|z=0
                                                               τ
                          " zX   j-1                                          Z z                 Z z ′
                                                 σ2m+2                                               H
                        +            θ2 j-(2m+2)        dz dz|z=0 + c3, j ) -     θ2 j-2 dz|z=0 -       ψ2 j-1 dz|z=0
                                m=0
                                                    τ                                                 τ
                          Z zX j-1                               j
                                               σ2m+2            X
                        -          θ2 j-(2m+2)        dz|z=0 -      β1,2m θ2 j-2m |z=0 ,
                              m=0
                                                  τ             m=1
```

<!-- Page 22 -->

Z z Z z Z zX j-1 H′ σ2m+2

```text
        c3, j   = (    θ2 j-2 dz|z=-1 +           ψ2 j-1 dz|z=-1 +            θ2 j-(2m+2)          dz|z=-1
                                               τ                         m=0
                                                                                              τ
                         " z                       " z ′                           " zX   j-1
                                                          H                                               σ2m+2
                  -β1,0 (       θ2 j-2 dz dz|z=0 +           ψ2 j-1 dz dz|z=0 +               θ2 j-(2m+2)        dz dz|z=0 )
                                                          τ                             m=0
                                                                                                             τ
                    Z z                   Z z ′                    Z zX j-1                                  j
                                               H                                        σ2m+2              X
                  -      θ2 j-2 dz|z=0 -          ψ2 j-1 dz|z=0 -           θ2 j-(2m+2)           dz|z=0 -     β1,2m θ2 j-2m |z=0
                                               τ                       m=0
                                                                                            τ              m=1
                          " z                       " z ′                            " zX    j-1
                                                           H                                                 σ2m+2
                  +β2,0 (       θ2 j-2 dz dz|z=-1 +            ψ2 j-1 dz dz|z=-1 +               θ2 j-(2m+2)       dz dz|z=-1
                                                            τ                               m=0
                                                                                                               τ
                         " z                       " z ′                           " zX   j-1
                                                          H                                               σ2m+2
                  +β1,0 (       θ2 j-2 dz dz|z=0 +           ψ2 j-1 dz dz|z=0 +               θ2 j-(2m+2)        dz dz|z=0 )
                                                          τ                             m=0
                                                                                                             τ
                    Z z                   Z z ′                    Z zX j-1                                  j
                                               H                                        σ2m+2              X
                  +      θ2 j-2 dz|z=0 +          ψ2 j-1 dz|z=0 +           θ2 j-(2m+2)           dz|z=0 +     β1,2m θ2 j-2m |z=0 )
                                               τ                       m=0
                                                                                            τ              m=1
                        j
                        X
                    +         β2,2m θ2 j-2m |z=-1 )/(β1,0 - β1,0 β2,0 - β2,0 ).
                        m=1
```

From θ2 j in the linear perturbation solution algorithm for case II " z " z ′ H

```text
  c2, j = -β1,0 (         θ2 j-2 dz dz|z=0 +           ψ2 j-1 dz dz|z=0
                                                    τ
               " zX   j-1                                           j                       Z z                 Z z ′
                                      σ2m+2                        X                                               H
            +             θ2 j-(2m+2)        dz dz|z=0 + c3, j ) -     β1,2m θ2 j-2m |z=0 -     θ2 j-2 dz|z=0 -       ψ2 j-1 dz|z=0
                     m=0
                                         τ                         m=1
                                                                                                                    τ
               Z zX j-1
                                    σ2m+2
            -           θ2 j-(2m+2)        dz|z=0 .
                   m=0
                                       τ
```

From ψ2 j+1 in the linear perturbation solution algorithm " zX j " z ′

```text
                c4, j = -                     D u2 j-2m R2m dz dz|z=-1 + S           θ2 j dz dz|z=-1
                                     m=0
                                  " zX
                                     j-2                                      " z
                              -               ψ2 j-2m-3 σ2m+2 dz dz|z=-1 -           ψ2 j-3 dz dz|z=-1 + c5, j + γ3 ψ′2 j-3 |z=-1 ,
                                        m=0
```

" zX j " z

```text
                     c5, j =                       D′ u2 j-2m R2m dz dz|z=0 - S         θ2 j dz dz|z=0
                                         m=0
                                        " zX
                                           j-2                                     " z
                                    +                ψ2 j-2m-3 σ2m+2 dz dz|z=0 +          ψ2 j-3 dz dz|z=0 - γ1 ψ′2 j-3 |z=0 ,
                                               m=0
```

& zX j & z ′

```text
                c6, j = -                          D u2 j-2m R2m dz dz dz dz|z=-1 + S              θ2 j dz dz dz dz|z=-1
                                        m=0
                                  " zX
                                     j-1                                      & zX
                                                                                 j-2
                              +               ψ2 j-2m-1 σ2m+2 dz dz|z=-1 -                     ψ2 j-2m-3 σ2m+2 dz dz dz dz|z=-1
                                   m=0                                                   m=0
                                  & z                                     " z
                                                                                                     1       1
                              -           ψ2 j-3 dz dz dz dz|z=-1 + 2             ψ2 j-1 dz dz|z=-1 - c4, j + c5, j + c7, j ,
                                                                                                     6       2
```

<!-- Page 23 -->

& zX j & z ′ c7, j = D u2 j-2m R2m dz dz dz dz|z=0 - S θ2 j dz dz dz dz|z=0

```text
                                    m=0
                                " zX
                                   j-1                                          & zX
                                                                                   j-2
                            -            ψ2 j-2m-1 σ2m+2 dz dz|z=0 +                            ψ2 j-2m-3 σ2m+2 dz dz dz dz|z=0
                                 m=0                                                      m=0
                                & z                                      " z
                            +           ψ2 j-3 dz dz dz dz|z=0 - 2                ψ2 j-1 dz dz|z=0 .
```

From the matrix N in the linear perturbation solution algorithm for case II

```text
                                                                 N1,1, j = u0 ,
                                                                     Z 0
                                                         N1,2, j =         ψ̆2 j-1 U ′ dz,
                                                                      -1
```

Z z " z Z z θ0 θ0 θ0

```text
                  N2,1, j   =        dz|z=-1 - β1,0          dz dz|z=0 -         dz|z=0
                                   τ                      τ                   τ
                                    " z                             " z                      Z z
                                          θ0                             θ0                      θ0
                              +β2,0          dz dz|z=-1 + β2,0 β1,0         dz dz|z=0 + β2,0        dz|z=0 ,
                                          τ                              τ                       τ
```

Z z Z z Z z H′ σ2

```text
             N2,2, j   =    θ̂2 j-2 dz|z=-1 +          ψ̆2 j-1     dz|z=-1 +         θ̂2 j-2    dz|z=-1
                                                                τ                             τ
                                 " z                          " z                               "    z
                                                                             H′                                 σ2
                         -β1,0 (        θ̂2 j-2 dz dz|z=0 +          ψ̆2 j-1     dz dz|z=0 +           θ̂2 j-2      dz dz|z=0 )
                                                                              τ                                  τ
                           Z z                     Z z                            z
                                                                H′
                                                                                Z
                                                                                             σ2
                         -      θ̂2 j-2 dz|z=0 -        ψ̆2 j-1     dz|z=0 -        θ̂2 j-2     dz|z=0
                                                                 τ                           τ
                                 " z                            " z                               " z
                                                                               H′                                  σ2
                         +β2,0 (        θ̂2 j-2 dz dz|z=-1 +           ψ̆2 j-1    dz dz|z=-1 +             θ̂2 j-2    dz dz|z=-1
                                                                               τ                                   τ
                                 " z                          " z                               "    z
                                                                             H′                                 σ2
                         +β1,0 (        θ̂2 j-2 dz dz|z=0 +          ψ̆2 j-1     dz dz|z=0 +           θ̂2 j-2      dz dz|z=0 )
                                                                              τ                                  τ
                           Z z                     Z z                            z
                                                                H′
                                                                                Z
                                                                                             σ2
                         +      θ̂2 j-2 dz|z=0 +        ψ̆2 j-1     dz|z=0 +        θ̂2 j-2     dz|z=0 )
                                                                 τ                           τ
                         -β1,2 (1 - β2,0 )θ̂2 j-2 |z=0 + β2,2 θ̂2 j-2 |z=-1 .
```

From the vector q in the linear perturbation solution algorithm for case II Z 0 Z 0

```text
                q1, j = -δ0, j-1 u0 -          ψ̂2 j-1 U ′ dz + R2 j-2          ψ̃2 j-1 U ′ dz + γ4 u2 j-4 |z=-1 - γ2 u2 j-4 |z=0 ,
                                          -1                               -1
```

<!-- Page 24 -->

Z z Z z Z zX j-2 H′ σ2m+2

```text
q2, j   = -      θ̃2 j-2 dz|z=-1 -       (ψ̂2 j-1 - R2 j-2 ψ̃2 j-1 ) dz|z=-1 -                  θ2 j-(2m+2)             dz|z=-1
                                                                       τ                   m=1
                                                                                                                  τ
            Z z                                " z                           " z
                         σ2                                                                                    H′
          -      θ̃2 j-2      dz|z=-1 + β1,0 (       θ̃2 j-2 dz dz|z=0 +          (ψ̂2 j-1 - R2 j-2 ψ̃2 j-1 ) dz dz|z=0
                          τ                                                                                     τ
            " zX    j-2                                     " z                                Z z
                                     σ2m+2                                  σ2
          +              θ2 j-(2m+2)          dz dz|z=0 +           θ̃2 j-2    dz dz|z=0 ) +        θ̃2 j-2 dz|z=0
                   m=1
                                        τ                                   τ
            Z z                                            Z zX  j-2                                   Z z
                                            H′                                    σ2m+2                             σ2
          +     (ψ̂2 j-1 - R2 j-2 ψ̃2 j-1 ) dz|z=0 +                  θ2 j-(2m+2)          dz|z=0 +         θ̃2 j-2      dz|z=0
                                            τ                   m=1
                                                                                     τ                               τ
                  " z                           " z                                                    " zX    j-2
                                                                                 H′                                              σ2m+2
          -β2,0 (        θ̃2 j-2 dz dz|z=-1 +        (ψ̂2 j-1 - R2 j-2 ψ̃2 j-1 ) dz dz|z=-1 +                       θ2 j-(2m+2)        dz dz|z=-1
                                                                                  τ                           m=1
                                                                                                                                   τ
            " z                                        " z                          " z
                           σ2                                                                                           H′
          +        θ̃2 j-2      dz dz|z=-1 - (-β1,0 (         θ̃2 j-2 dz dz|z=0 +          (ψ̂2 j-1 - R2 j-2 ψ̃2 j-1 ) dz dz|z=0
                            τ                                                                                            τ
            " zX    j-2                                     " z                                Z z
                                     σ2m+2                                  σ2
          +              θ2 j-(2m+2)          dz dz|z=0 +           θ̃2 j-2    dz dz|z=0 ) -        θ̃2 j-2 dz|z=0
                   m=1
                                        τ                                   τ
            Z z                                            Z zX  j-2                                   Z z
                                            H′                                    σ2m+2                             σ2
          -     (ψ̂2 j-1 - R2 j-2 ψ̃2 j-1 ) dz|z=0 -                  θ2 j-(2m+2)          dz|z=0 -         θ̃2 j-2      dz|z=0 ))
                                            τ                   m=1
                                                                                     τ                               τ
                 j
                 X                                                           j
                                                                             X
            +(       β1,2m θ2 j-2m |z=0 + β1,2 θ̃2 j-2 |z=0 )(1 - β2,0 ) -         β2,2m θ2 j-2m |z=-1 - β2,2 θ̃2 j-2 |z=-1 .
               m=2                                                           m=2
```

### 7.3 Nonlinear perturbation solution details

From Ψ0 " z " z

```text
                                                           ∂u0               ∂Θ0
                           c4 (Y, T ) = -          D′ R0       dz dz|z=-1 + S    dz dz|z=-1 + c5 (Y, T ),
                                                           ∂Y                 ∂Y
                                                  " z                        " z
                                                       ′   ∂u0                   ∂Θ0
                                     c5 (Y, T ) =     D R0     dz dz|z=0 - S         dz dz|z=0 ,
                                                           ∂Y                    ∂Y
                                           & z                                   & z
                                                       ∂u0
                                                       ′                             ∂Θ0
                      c6 (Y, T ) = -              D R0      dz dz dz dz|z=-1 + S          dz dz dz dz|z=-1
                                                       ∂Y                             ∂Y
                                        1            1
                                       - c4 (Y, T ) + c5 (Y, T ) + c7 (Y, T ),
                                        6            2
                                       & z                                     & z
                                                  ′  ∂u0                           ∂Θ0
                          c7 (Y, T ) =          D R0      dz dz dz dz|z=0 - S          dz dz dz dz|z=0 .
                                                      ∂Y                           ∂Y
```

From u2 Z z ∂Ψ0 ′

```text
                                                        c8 (Y, T ) =             U dz|z=0 .
                                                                              ∂Y
```

From Θ2 for case I Z z " z

```text
                                                 ∂Ψ0 H ′                               ∂Ψ0 H ′
                           c10 (Y, T ) =                 dz|z=0 - β1,0 (-                      dz dz|z=0 + c11 (Y, T )),
                                                  ∂Y τ                                  ∂Y τ
```

Z z Z z " z

```text
                                        ∂Ψ0 H ′                 ∂Ψ0 H ′                       ∂Ψ0 H ′
                     c11 (Y, T ) = (             dz|z=-1 -               dz|z=0 - β1,0                dz dz|z=0
                                         ∂Y τ                    ∂Y τ                          ∂Y τ
                                         " z                                 z
                                               ∂Ψ0 H ′                         ∂Ψ0 H ′
                                                                          Z
                                   +β2,0                dz dz|z=-1 + β2,0               dz|z=0
                                                ∂Y τ                            ∂Y τ
                                              " z
                                                   ∂Ψ0 H ′
                                   +β2,0 β1,0               dz dz|z=0 )/(-β1,0 + β2,0 β1,0 + β2,0 ).
                                                    ∂Y τ
```

<!-- Page 25 -->

From Θ2 for case II Z z Z z 2 Z z

```text
                                      1 ∂Θ0               ∂ Θ0                ∂Ψ0 H ′
             c10 (Y, T ) = -                dz|z=0 +             dz|z=0 +               dz|z=0
                                      τ ∂T                 ∂Y 2                ∂Y τ
                                Z z                        Z z
                                      1 ∂Ψ0 ∂Θ0                 1 ∂Ψ0 ∂Θ0
                              +                  dz|z=0 -                  dz|z=0
                                      τ ∂Y ∂z                   τ ∂z ∂Y
                                     " z                       " z 2                    " z
                                           1 ∂Θ0                     ∂ Θ0                      ∂Ψ0 H ′
                              -β1,0 (             dz dz|z=0 -            2
                                                                           dz dz| z=0 -                  dz dz|z=0
                                           τ ∂T                       ∂Y                        ∂Y τ
                                " z                             " z
                                       1 ∂Ψ0 ∂Θ0                      1 ∂Ψ0 ∂Θ0
                              -                    dz dz|z=0 +                     dz dz|z=0 + c11 (Y, T )) - β1,2 Θ0 |z=0 .
                                        τ ∂Y ∂z                       τ ∂z ∂Y
```

For case II, the second coupled nonlinear partial differential equation for u0 and c3 (Y, T ) at O(l2 ) is Z z Z z 2 Z z Z z

```text
                 1 ∂Θ0                   ∂ Θ0                    ∂Ψ0 H ′                     1 ∂Ψ0 ∂Θ0
                        dz|z=-1 -            2
                                                dz|z=-1 -                  dz|z=-1 -                       dz|z=-1
                 τ ∂T                    ∂Y                       ∂Y τ                       τ ∂Y ∂z
              Z z                                 " z                          " z 2                       " z
                   1 ∂Ψ0 ∂Θ0                             1 ∂Θ0                         ∂ Θ0                       ∂Ψ0 H ′
            +                   dz|z=-1 - β1,0 (                dz dz|z=0 -                   dz dz|z=0  -                 dz dz|z=0
                   τ ∂z ∂Y                               τ ∂T                           ∂Y 2                       ∂Y τ
              " z                               " z                                  Z z                    Z z 2
                     1 ∂Ψ0 ∂Θ0                          1 ∂Ψ0 ∂Θ0                         1 ∂Θ0                  ∂ Θ0
            -                    dz dz|z=0 +                         dz dz|z=0 ) -                dz|z=0 +             dz|z=0
                     τ ∂Y ∂z                            τ ∂z ∂Y                           τ ∂T                   ∂Y 2
              Z z                         z                             z                                "   z
                   ∂Ψ0 H ′
                                      Z                              Z
                                            1 ∂Ψ0 ∂Θ0                     1 ∂Ψ0 ∂Θ0                            1 ∂Θ0
            +               dz|z=0 +                      dz|z=0 -                       dz|z=0 + β2,0 (              dz dz|z=-1
                    ∂Y τ                    τ ∂Y ∂z                       τ ∂z ∂Y                              τ ∂T
              " z 2                        " z                             "     z
                     ∂ Θ0                         ∂Ψ0 H ′                          1 ∂Ψ0 ∂Θ0
            -              dz dz| z=-1  -                   dz dz| z=-1  -                       dz dz|z=-1
                     ∂Y 2                          ∂Y τ                            τ ∂Y ∂z
              " z                                      " z                           " z 2
                     1 ∂Ψ0 ∂Θ0                                1 ∂Θ0                         ∂ Θ0
            +                    dz dz|z=-1 + β1,0 (                   dz dz|z=0 -                 dz dz|z=0
                     τ ∂z ∂Y                                  τ ∂T                           ∂Y 2
              " z                           "    z                             "     z
                     ∂Ψ0 H ′                       1 ∂Ψ0 ∂Θ0                           1 ∂Ψ0 ∂Θ0
            -                dz dz|z=0 -                         dz dz|z=0 +                         dz dz|z=0 )
                      ∂Y τ                         τ ∂Y ∂z                             τ ∂z ∂Y
              Z z                    Z z 2                  Z z                        Z z
                   1 ∂Θ0                  ∂ Θ0                  ∂Ψ0 H ′                    1 ∂Ψ0 ∂Θ0
            +             dz|z=0 -               dz| z=0  -                dz| z=0 -                     dz|z=0
                   τ ∂T                    ∂Y 2                  ∂Y τ                      τ ∂Y ∂z
              Z z
                   1 ∂Ψ0 ∂Θ0
            +                   dz|z=0 ) - β1,2 (1 - β2,0 )Θ0 |z=0 + β2,2 Θ0 |z=-1 = 0.
                   τ ∂z ∂Y
```

The constants appearing in the differential equations (4.1.44), (4.1.45) are Z 0 2 2

```text
                                                    a p = p - R0 p           Ψ̃0 U ′ dz,
                                                                            -1
                                                                      Z 0
                                                         b p = S p2         Ψ̂0 U ′ dz,
                                                                       -1
                                                                       Z 0
                                                       c p = -R0 p2           Ψ̃0 H ′ dz,
                                                                         -1
                                                         Z 0
                                            d p = S p2         Ψ̂0 H ′ dz - τ(β2,2 - β1,2 - p2 ).
                                                          -1
```

From Ψ2 " z

```text
                                       ∂Ψ0           ∂ 2 Ψ0                         ∂u0
                   c12 (Y, T ) =           |z=-1 - 2        |z=-1  -          D′ R2      dz dz|z=-1
                                       ∂T             ∂Y 2                          ∂Y
                                        " z                                  " z
                                                    ∂u2                             ∂Θ2
                                      -       D′ R0       dz dz|z=-1 + S                  dz dz|z=-1
                                                     ∂Y                              ∂Y
                                        " z                                  " z
                                               ∂Ψ0 ∂3 Ψ0                            ∂Ψ0 ∂3 Ψ0
                                      +                      dz dz| z=-1   -                    dz dz|z=-1 + c13 (Y, T ),
                                                ∂z ∂Y∂z2                             ∂Y ∂z3
                                                                              " z
                                                ∂Ψ0            ∂2 Ψ0                       ∂u0
                            c13 (Y, T ) = -         |z=0 + 2       2
                                                                      | z=0 +        D′ R2      dz dz|z=0
                                                ∂T              ∂Y                          ∂Y
                                                " z                                 " z
                                                              ∂u2                          ∂Θ2
                                              +       D′ R0        dz dz|z=0 - S                dz dz|z=0
                                                              ∂Y                           ∂Y
                                                " z                                 "    z
                                                      ∂Ψ0 ∂3 Ψ0                            ∂Ψ0 ∂3 Ψ0
                                              -                       dz  dz|z=0  +                    dz dz|z=0 ,
                                                        ∂z ∂Y∂z2                            ∂Y ∂z3
```

<!-- Page 26 -->

" z " z 2 & z

```text
                            ∂Ψ0                         ∂ Ψ0                              ∂u0
          c14 (Y, T ) =            dz dz|z=-1 -       2      2
                                                                dz dz|z=-1 -        D′ R2      dz dz dz dz|z=-1
                             ∂T                          ∂Y                               ∂Y
                          & z                                        & z
                                         ∂u2                                 ∂Θ2
                        -          D′ R0     dz dz dz dz|z=-1 + S                dz dz dz dz|z=-1
                                          ∂Y                                 ∂Y
                          & z                                        & z
                                   ∂Ψ0 ∂3 Ψ0                                 ∂Ψ0 ∂3 Ψ0
                        +                       dz dz dz dz| z=-1  -                   dz dz dz dz|z=-1
                                     ∂z ∂Y∂z2                                 ∂Y ∂z3
                          1             1
                        - c12 (Y, T ) + c13 (Y, T ) + c15 (Y, T ),
                          6             2
```

" z " z 2 & z

```text
                               ∂Ψ0                      ∂ Ψ0                               ∂u0
           c15 (Y, T ) = -          dz dz|z=0 +       2      2
                                                                dz dz| z=0 +         D′ R2      dz dz dz dz|z=0
                               ∂T                        ∂Y                                ∂Y
                              & z                                 & z
                                        ∂u2                                  ∂Θ2
                            +     D′ R0     dz dz dz dz|z=0 - S                  dz dz dz dz|z=0
                                        ∂Y                                   ∂Y
                              & z                                 &        z
                                  ∂Ψ0 ∂3 Ψ0                                  ∂Ψ0 ∂3 Ψ0
                            -                 dz dz dz  dz| z=0 +                      dz dz dz dz|z=0 .
                                   ∂z ∂Y∂z2                                   ∂Y ∂z3
```

From u4 Z z Z z 2 Z z

```text
                                   ∂u2               ∂ u2               ∂Ψ2 ′
               c16 (Y, T ) = -         dz|z=0 +         2
                                                          dz|z=0 +          U dz|z=0
                                   ∂T                ∂Y                  ∂Y
                               Z z                   Z z                    Z z
                                   ∂Ψ0 ∂u2                ∂Ψ2 ∂u0               ∂Ψ0 ∂u2
                             +              dz|z=0 -               dz|z=0 -             dz|z=0 - γ2 u0 .
                                    ∂Y ∂z                  ∂z ∂Y                 ∂z ∂Y
```

From Θ4 for case I Z z Z z 2 Z z

```text
                                   1 ∂Θ2               ∂ Θ2                ∂Ψ2 H ′
          c18 (Y, T ) = -                dz|z=0 +          2
                                                              dz|z=0 +               dz|z=0
                                   τ ∂T                 ∂Y                  ∂Y τ
                             Z z                        Z z
                                   1 ∂Ψ0 ∂Θ2                 1 ∂Ψ0 ∂Θ2
                           +                  dz|z=0 -                  dz|z=0
                                   τ ∂Y ∂z                   τ ∂z ∂Y
                                  " z                       " z 2                    " z
                                        1 ∂Θ2                     ∂ Θ2                      ∂Ψ2 H ′
                           -β1,0 (             dz dz|z=0 -              dz dz| z=0 -                  dz dz|z=0
                                        τ ∂T                       ∂Y 2                      ∂Y τ
                             " z                             " z
                                    1 ∂Ψ0 ∂Θ2                      1 ∂Ψ0 ∂Θ2
                           -                    dz dz|z=0 +                     dz dz|z=0 + c19 (Y, T )) - β1,2 Θ2 |z=0 ,
                                     τ ∂Y ∂z                       τ ∂z ∂Y
```

Z z Z z 2 Z z Z z

```text
                              1 ∂Θ2                 ∂ Θ2                    ∂Ψ2 H ′                   1 ∂Ψ0 ∂Θ2
       c19 (Y, T ) = (-              dz|z=-1 +           2
                                                           dz|z=-1  +                 dz| z=-1  +                   dz|z=-1
                              τ ∂T                   ∂Y                      ∂Y τ                      τ ∂Y ∂z
                       Z z                             Z z                     Z z 2                 Z z
                             1 ∂Ψ0 ∂Θ2                      1 ∂Θ2                  ∂ Θ2                   ∂Ψ2 H ′
                     -                   dz|z=-1 - (-               dz|z=0 +            2
                                                                                           dz| z=0 +                dz|z=0
                             τ ∂z ∂Y                        τ ∂T                    ∂Y                     ∂Y τ
                       Z z                         Z z
                             1 ∂Ψ0 ∂Θ2                  1 ∂Ψ0 ∂Θ2
                     +                   dz|z=0 -                    dz|z=0
                             τ ∂Y ∂z                    τ ∂z ∂Y
                             " z                       " z 2                      " z
                                   1 ∂Θ2                     ∂ Θ2                         ∂Ψ2 H ′
                     -β1,0 (              dz dz|z=0 -                dz dz|z=0  -                  dz dz|z=0
                                   τ ∂T                       ∂Y 2                         ∂Y τ
                       " z                              " z
                               1 ∂Ψ0 ∂Θ2                      1 ∂Ψ0 ∂Θ2
                     -                     dz dz|z=0 +                      dz dz|z=0 ) - β1,2 Θ2 |z=0 )(1 - β2,0 )
                               τ ∂Y ∂z                        τ ∂z ∂Y
                            " z                         " z 2                        " z
                                   1 ∂Θ2                       ∂ Θ2                         ∂Ψ2 H ′
                     -β2,0 (              dz dz|z=-1 -                dz dz| z=-1 -                   dz dz|z=-1
                                   τ ∂T                         ∂Y 2                         ∂Y τ
                       " z                               " z
                               1 ∂Ψ0 ∂Θ2                        1 ∂Ψ0 ∂Θ2
                     -                     dz dz|z=-1 +                       dz dz|z=-1 )
                               τ ∂Y ∂z                          τ ∂z ∂Y
                     -β2,2 Θ2 |z=-1 )/(-β1,0 + β2,0 β1,0 + β2,0 ).
```

<!-- Page 27 -->

From Θ4 for case II Z z Z z 2 Z z Z z

```text
                       1 ∂Θ2                   ∂ Θ2                ∂Ψ2 H ′               1 ∂Ψ2 ∂Θ0
  c18 (Y, T ) = -              dz|z=0 +           2
                                                     dz|z=0 +               dz|z=0 +                  dz|z=0
                       τ ∂T                    ∂Y                   ∂Y   τ                τ ∂Y ∂z
                  Z z                          Z z                        Z z                              " z
                       1 ∂Ψ0 ∂Θ2                    1 ∂Ψ2 ∂Θ0                  1 ∂Ψ0 ∂Θ2                        1 ∂Θ2
                +                   dz|z=0 -                    dz|z=0 -                   dz|z=0 - β1,0 (             dz dz|z=0
                       τ ∂Y ∂z                      τ ∂z ∂Y                    τ ∂z ∂Y                          τ ∂T
                  " z 2                       " z                        " z
                        ∂ Θ2                        ∂Ψ2 H ′                     1 ∂Ψ2 ∂Θ0
                -               dz dz|z=0  -                 dz dz|z=0 -                    dz dz|z=0
                          ∂Y 2                       ∂Y τ                       τ ∂Y ∂z
                  " z                               " z                            " z
                        1 ∂Ψ0 ∂Θ2                         1 ∂Ψ2 ∂Θ0                      1 ∂Ψ0 ∂Θ2
                -                     dz dz|z=0 +                      dz dz|z=0 +                    dz dz|z=0 + c19 (Y, T ))
                         τ ∂Y ∂z                          τ ∂z ∂Y                        τ ∂z ∂Y
                -β1,2 Θ2 |z=0 - β1,4 Θ0 |z=0 .
```

For case II, the second coupled partial differential equation in terms of c9 (Y, T ) and c11 (Y, T ) at O(l4 ) is Z z Z z 2 Z z Z z

```text
             1 ∂Θ2                 ∂ Θ2                     ∂Ψ2 H ′                    1 ∂Ψ2 ∂Θ0
                    dz|z=-1 -           2
                                           dz|z=-1 -                  dz|z=-1 -                     dz|z=-1
             τ ∂T                   ∂Y                       ∂Y τ                      τ ∂Y ∂z
          Z z                          Z z                            Z z
               1 ∂Ψ0 ∂Θ2                    1 ∂Ψ2 ∂Θ0                      1 ∂Ψ0 ∂Θ2
       -                   dz|z=-1 +                      dz|z=-1 +                      dz|z=-1
               τ ∂Y ∂z                      τ ∂z ∂Y                        τ ∂z ∂Y
          Z z                  Z z 2                   Z z                        z
                                                            ∂Ψ2 H ′
                                                                               Z
               1 ∂Θ2                ∂ Θ2                                            1 ∂Ψ2 ∂Θ0
       -              dz|z=0 +              dz| z=0 +                 dz|z=0 +                    dz|z=0
               τ ∂T                  ∂Y 2                    ∂Y τ                    τ ∂Y ∂z
          Z z                        Z z                           Z z
               1 ∂Ψ0 ∂Θ2                   1 ∂Ψ2 ∂Θ0                    1 ∂Ψ0 ∂Θ2
       +                   dz|z=0 -                      dz|z=0 -                     dz|z=0
               τ ∂Y ∂z                     τ ∂z ∂Y                      τ ∂z ∂Y
               " z                        " z 2                       " z                           " z
                     1 ∂Θ2                       ∂ Θ2                        ∂Ψ2 H ′                      1 ∂Ψ2 ∂Θ0
       -β1,0 (              dz dz|z=0 -                  dz dz|z=0 -                   dz dz|z=0 -                    dz dz|z=0
                     τ ∂T                         ∂Y 2                        ∂Y τ                        τ ∂Y ∂z
          " z                              " z                               " z
                1 ∂Ψ0 ∂Θ2                         1 ∂Ψ2 ∂Θ0                         1 ∂Ψ0 ∂Θ2
       -                     dz dz|z=0 +                        dz dz|z=0 +                       dz dz|z=0 )
                 τ ∂Y ∂z                          τ ∂z ∂Y                           τ ∂z ∂Y
               " z                         " z 2                         " z                            " z
                     1 ∂Θ2                         ∂ Θ2                        ∂Ψ2 H ′                         1 ∂Ψ2 ∂Θ0
       +β2,0 (              dz dz|z=-1 -               2
                                                          dz dz| z=-1 -                   dz dz|z=-1 -                    dz dz|z=-1
                     τ ∂T                           ∂Y                          ∂Y τ                           τ ∂Y ∂z
          " z                               " z                                 " z
                1 ∂Ψ0 ∂Θ2                           1 ∂Ψ2 ∂Θ0                          1 ∂Ψ0 ∂Θ2
       -                     dz dz|z=-1 +                         dz dz|z=-1 +                      dz dz|z=-1
                 τ ∂Y ∂z                            τ ∂z ∂Y                            τ ∂z ∂Y
          Z z                  Z z 2                   Z z                        z
                                                            ∂Ψ2 H ′
                                                                               Z
               1 ∂Θ2                ∂ Θ2                                            1 ∂Ψ2 ∂Θ0
       +              dz|z=0 -           2
                                            dz| z=0 -                 dz|z=0 -                    dz|z=0
               τ ∂T                  ∂Y                      ∂Y τ                    τ ∂Y ∂z
          Z z                        Z z                           Z z
               1 ∂Ψ0 ∂Θ2                   1 ∂Ψ2 ∂Θ0                    1 ∂Ψ0 ∂Θ2
       -                   dz|z=0 +                      dz|z=0 +                     dz|z=0
               τ ∂Y ∂z                     τ ∂z ∂Y                      τ ∂z ∂Y
               " z                        " z 2                       " z                           " z
                     1 ∂Θ2                       ∂ Θ2                        ∂Ψ2 H ′                      1 ∂Ψ2 ∂Θ0
       +β1,0 (              dz dz|z=0 -                  dz dz|z=0 -                   dz dz|z=0 -                    dz dz|z=0
                     τ ∂T                         ∂Y 2                        ∂Y τ                        τ ∂Y ∂z
          " z                              " z                               " z
                1 ∂Ψ0 ∂Θ2                         1 ∂Ψ2 ∂Θ0                         1 ∂Ψ0 ∂Θ2
       -                     dz dz|z=0 +                        dz dz|z=0 +                       dz dz|z=0 ))
                 τ ∂Y ∂z                          τ ∂z ∂Y                           τ ∂z ∂Y
       +(-β1,2 Θ2 |z=0 - β1,4 Θ0 |z=0 )(1 - β2,0 ) + β2,2 Θ2 |z=-1 + β2,4 Θ0 |z=-1 = 0.
```

### 7.4 A theorem for a class of nonlinear differential equations

The following Theorem A formalises a procedure outlined in Hildebrand (1956): Theorem A Provided that the L + 1 term Maclaurin series of the exact general solution, L X dl A xl

```text
                                                         A=               |
                                                                        x=0                                                    (7.4.1)
                                                              l=0
                                                                    dxl     l!
```

to an M th order ordinary differential equation dM A

```text
                                                                 =ξ                                                         (7.4.2)
                                                           dx M
```

exists and all the derivatives and integrals of A are defined at x = 0, it only solves the coefficients of xl , l ∈ {0, 1, . . . , L- M} in the residual of (7.4.2) provided ξ is expandable in a Maclaurin series as ∞ X dl ξ xl ξ= | l x=0

```text
                                                                                    ,                                          (7.4.3)
                                                              l=0
                                                                    dx           l!
```

<!-- Page 28 -->

M where all the derivatives and integrals of ξ are defined at x = 0 and the right hand side of (7.4.2) does not contain ddxMA . Proof of Theorem A Since the Maclaurin series of A and ξ exist and all their derivatives and integrals are defined at x = 0, we can integrate (7.4.2) M times and substitute the result into (7.4.1) to find L X d(l-M) ξ xl

```text
                                                      A=                    |
                                                                       (l-M) x=0
                                                                                      .                              (7.4.4)
                                                           l=0
                                                                  dx               l!
```

Substituting (7.4.4) into the residual r of (7.4.2) then gives L ∞ X d(l-M) ξ xl-M X dl ξ xl

```text
                                      r=                   |
                                                      (l-M) x=0
                                                                          -         | x=0    ,                       (7.4.5)
                                           l=0
                                                 dx               (l - M)! l=0 dxl        l!
```

provided ξ is expandable in a Maclaurin series as in (7.4.3). Equating like powers of x in (7.4.5) then yields ∞ X dl ξ xl

```text
                                                   r=-                  | x=0    ,                                   (7.4.6)
                                                           l=L-M+1
                                                                   dxl        l!
```

which shows that Theorem A is true. 

### 7.5 Another theorem for a class of nonlinear differential equations

The following Theorem B is of the essence of that given in various texts (see for example Muscalu & Schlag, 2013): Theorem B Provided that the 2L + 1 term complex Fourier series of the exact general solution L X

```text
                                           A=            P(A, einlx )einlx , 0 < l < ∞,                              (7.5.1)
                                                  n=-L
```

to an M th order ordinary differential equation dM A

```text
                                                                  = ξ,                                               (7.5.2)
                                                             dx M
```

exists, it only solves the coefficients of einlx for n ∈ [-L, L] in the residual of (7.5.2) if ξ is expandable as a complex Fourier series as X∞

```text
                                            ξ=         P(ξ, einlx )einlx , 0 < l < ∞.                                 (7.5.3)
                                                 n=-∞
```

Here A and ξ are periodic with period 2πl and all of their derivatives and integrals are continuous for all x. Moreover M the right hand side of (7.5.2) must not contain ddxMA and P(a, einlx ) denotes the projection of a onto einlx . Proof of Theorem B Since the complex Fourier series of A and ξ exist and because A and ξ are periodic with period 2πl and all their

```text
derivatives and integrals are continuous for all x, we can integrate (7.5.2) M times and substitute the result into (7.5.1)
```

to find L X d(-M) ξ

```text
                                             A=         P( (-M) , einlx )einlx ,                                    (7.5.4)
                                                   n=-L
                                                          dx
                   d  (-M)
                         ξ
```

where the notation dx th integral of ξ with respect to x. Substituting (7.5.4) into the residual r of (-M) denotes the M (7.5.2) then gives L ∞ d M X d(-M) ξ inlx inlx X

```text
                               r= M           P( (-M) , e )e -            P(ξ, einlx )einlx ,                    (7.5.5)
                                   dx n=-L dx                       n=-∞
```

provided ξ is expandable in complex Fourier series as in (7.5.3). Then equation (7.5.5) can be written as X

```text
                                            r=-            P(ξ, einlx )einlx ,                                       (7.5.6)
                                                          n<[-L,L]
```

which shows that Theorem B is true. 

<!-- Page 29 -->

Acknowledgement I thank Professor William Phillips for helping me with this paper.

## References

[1] Cox, S. M., & Leibovich, S. 1993. Langmuir circulations in a surface layer bounded by a strong thermocline. J. Phys. Ocean, 23, 1330-1345. [2] Craik, A.D.D., & Leibovich, S. 1976. A rational model for Langmuir circulations, J. Fluid. Mech., 73, 401-426. [3] Hardy, G. H. 1949. Divergent series. Oxford University Press. [4] Hayes, D. T., & Phillips, W. R. C. 2016. An asymptotic study of instability to Langmuir circulation in shallow layers. Geophy. Astrophys. Fluid Dyn, 110, 295-316. [5] Hayes, D. T., & Phillips, W. R. C. 2017. Nonlinear steady states to Langmuir circulation in shallow layers: an asymptotic study. Geophys. Astrophys. Fluid Dyn, 111, 65-90. [6] Hildebrand, F. B. 1956. Introduction to Numerical Analysis. McGraw-Hill. [7] Langmuir, I. 1938. Surface motion of water induced by wind. Science, 87, 119-123. [8] Leibovich, S. 1980. On wave-current interaction theories of Langmuir circulation. J. Fluid Mech., 99, 715-724. [9] Leibovich, S. 1983. The form and dynamics of Langmuir circulations. Ann. Rev. Fluid Mech., 15, 391-427. [10] Muscalu, C., & Schlag, W. 2013. Classical and Multilinear Harmonic Analysis. Vol 1. Cambridge University Press. [11] Plueddemann, A. J., Smith, J. A., Farmer, D. M., Weller, R. A., Crawford, W. R., Pinlel, R., Vagle, S. & Gnanade- silan, A. 1996. Structure and variability of Langmuir circulation during the Surface Waves Processes Program. J. Geophys. Res., 21, 85-102. [12] Smith, J. A. 1992. Observed growth of Langmuir circulation. J. Geophys. Res., 97, 5651-5664. [13] Thorpe, S. A. 2004. Langmuir circulation. Annu. Rev. Fluid Mech., 36, 55-79.
