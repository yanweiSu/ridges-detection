### MHRD code: './TF_analysis/MultiCurveExt.m'
- There are two functions in it, `CurveMultiExt_init_3curves` and `CurveMultiExt_init_2curves`.
- They are written in C, and need to be compiled using mex to generate `.mexw64` files beforehand.

### Numerical simulation code
- `SingleCurveExttest.m`, `FullyAdaptiveRD_Marcelo.m`, `RRP_RDtest.m` and `MultiCurveExttest.m`: Generate extracted curves and rmse of Single-RD, FM-RD, RRP-RD and MHRD
- `compare.m`: load the generated rmse and show the statistical result
