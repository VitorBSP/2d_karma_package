# Monte Carlo results, 1,000 replications per cell

Intermediate-size run of the Monte Carlo study of Section 4.1. This is **not**
the article's configuration of 10,000 replications per cell, which was not
run in the audit environment. Command:

```
Rscript analysis/simulations/run_monte_carlo.R --mode=full --replications=1000 \
  --cores=3 --output=analysis/output/simulations_r1000
Rscript analysis/simulations/summarize_monte_carlo.R \
  --input=analysis/output/simulations_r1000
```

Settings: the article's true values (scenario set `paper`), master seed
20260101, karma2d 0.1.0.9000, R 4.3.3. The run took 1,494 s of wall time on
3 cores (522 s for KARMA(1,1) 50 x 50). The full study should therefore take
about 4 h on similar hardware. All 9,000 fits converged; none failed.
Machine-readable results: `docs/monte-carlo-r1000-comparison.csv`.

Columns: mean of the estimates (`mean`) and the article's value (`paper`);
`z` = difference divided by its Monte Carlo standard error (sd / sqrt(1000));
MSE = mean squared error; coverage of nominal 95 % Wald intervals with
standard errors from the expected (`CR exp`) and observed (`CR obs`)
information, and the article's coverage (`CR paper`). Monte Carlo standard
error of a coverage near 0.95: 0.007.

## Summary

* **Means.** 51 of 54 are within +/-2 Monte Carlo standard errors of the
  article (KARMA(1,0) 10 x 10 alpha: z = 3.5; KARMA(1,1) 10 x 10 phi(1,1):
  -2.1; KARMA(0,1) 50 x 50 theta(0,1): 2.1). This is consistent with the
  article's Tables 1-3.
* **MSE.** Matches the article closely in every cell except the 10 x 10
  images. There the MSE is larger, markedly so for KARMA(1,1) alpha (0.133 vs
  0.017) and KARMA(1,0) alpha (0.054 vs 0.038). This points to a few
  extreme estimates on 10 x 10 images that do not appear in the article's
  numbers. Possible causes: different replications, exclusion of fits, or the
  MSE definition `bias^2 + sd^2` of the original scripts, which is
  numerically equivalent and therefore cannot explain a factor of 8.
* **Coverage.** Reproduces the article's pattern, including the
  under-coverage for 10 x 10 images and for phi(1,1) and theta(1,1) in the
  ARMA model. Expected-information intervals agree with the article to within
  about 0.02 for 25 x 25 and 50 x 50. The legacy mixed matrix gives the same
  coverage as the expected information for the gamma parameters, and almost
  the same for the shape.
* The KARMA(1,0) and KARMA(0,1) scenarios use the article's true values,
  which differ from those in the committed scripts (see the technical audit,
  Sec. 5.1).

## Full table

| Scenario | Size | Parameter | Mean | Paper | z | MSE | MSE paper | CR exp | CR obs | CR paper |
|---|---|---|---|---|---|---|---|---|---|---|
| karma01 | 10 | alpha | 1.1141 | 1.1134 | 0.3 | 0.00547 | 0.0054 | 0.923 | 0.926 | 0.924 |
| karma01 | 10 | shape | 12.6360 | 12.6818 | -1.2 | 1.62 | 1.62 | 0.936 | 0.938 | 0.941 |
| karma01 | 10 | theta(0,1) | 0.3406 | 0.3459 | -1.6 | 0.0111 | 0.0124 | 0.906 | 0.932 | 0.887 |
| karma01 | 10 | theta(1,0) | 0.2820 | 0.2844 | -0.7 | 0.0135 | 0.0131 | 0.875 | 0.904 | 0.885 |
| karma01 | 10 | theta(1,1) | 0.1856 | 0.1810 | 1.2 | 0.016 | 0.0163 | 0.885 | 0.909 | 0.882 |
| karma01 | 25 | alpha | 1.1105 | 1.1113 | -0.9 | 0.000804 | 0.0008 | 0.956 | 0.954 | 0.946 |
| karma01 | 25 | shape | 12.3623 | 12.3501 | 0.9 | 0.204 | 0.195 | 0.937 | 0.937 | 0.941 |
| karma01 | 25 | theta(0,1) | 0.3442 | 0.3460 | -1.7 | 0.00112 | 0.0011 | 0.933 | 0.938 | 0.943 |
| karma01 | 25 | theta(1,0) | 0.2827 | 0.2820 | 0.7 | 0.00123 | 0.0012 | 0.932 | 0.937 | 0.944 |
| karma01 | 25 | theta(1,1) | 0.1812 | 0.1826 | -1.2 | 0.00134 | 0.0013 | 0.939 | 0.944 | 0.948 |
| karma01 | 50 | alpha | 1.1109 | 1.1112 | -0.7 | 0.000191 | 0.0002 | 0.941 | 0.941 | 0.948 |
| karma01 | 50 | shape | 12.3519 | 12.3483 | 0.6 | 0.0443 | 0.0472 | 0.945 | 0.944 | 0.942 |
| karma01 | 50 | theta(0,1) | 0.3470 | 0.3459 | 2.1 | 0.000263 | 0.0002 | 0.944 | 0.947 | 0.948 |
| karma01 | 50 | theta(1,0) | 0.2822 | 0.2821 | 0.1 | 0.00026 | 0.0003 | 0.954 | 0.957 | 0.949 |
| karma01 | 50 | theta(1,1) | 0.1844 | 0.1843 | 0.2 | 0.000283 | 0.0003 | 0.951 | 0.954 | 0.948 |
| karma10 | 10 | alpha | 0.7317 | 0.7073 | 3.5 | 0.0544 | 0.0376 | 0.919 | 0.942 | 0.929 |
| karma10 | 10 | phi(0,1) | 0.2643 | 0.2686 | -1.4 | 0.0098 | 0.0093 | 0.922 | 0.937 | 0.922 |
| karma10 | 10 | phi(1,0) | 0.1749 | 0.1761 | -0.4 | 0.00965 | 0.0095 | 0.927 | 0.949 | 0.926 |
| karma10 | 10 | phi(1,1) | 0.1581 | 0.1633 | -1.7 | 0.00955 | 0.0102 | 0.941 | 0.958 | 0.927 |
| karma10 | 10 | shape | 20.5253 | 20.4864 | 0.6 | 4.67 | 4.83 | 0.944 | 0.949 | 0.936 |
| karma10 | 25 | alpha | 0.6714 | 0.6706 | 0.3 | 0.00619 | 0.0049 | 0.948 | 0.949 | 0.949 |
| karma10 | 25 | phi(0,1) | 0.2776 | 0.2796 | -1.8 | 0.00116 | 0.0011 | 0.939 | 0.942 | 0.948 |
| karma10 | 25 | phi(1,0) | 0.1876 | 0.1877 | -0.1 | 0.00119 | 0.0011 | 0.938 | 0.945 | 0.945 |
| karma10 | 25 | phi(1,1) | 0.1654 | 0.1639 | 1.4 | 0.00126 | 0.0012 | 0.945 | 0.942 | 0.945 |
| karma10 | 25 | shape | 19.7549 | 19.7205 | 1.5 | 0.536 | 0.529 | 0.958 | 0.958 | 0.946 |
| karma10 | 50 | alpha | 0.6629 | 0.6638 | -0.8 | 0.00138 | 0.0013 | 0.954 | 0.955 | 0.948 |
| karma10 | 50 | phi(0,1) | 0.2817 | 0.2810 | 1.3 | 0.000247 | 0.0003 | 0.955 | 0.958 | 0.950 |
| karma10 | 50 | phi(1,0) | 0.1890 | 0.1890 | -0.1 | 0.000284 | 0.0003 | 0.942 | 0.944 | 0.950 |
| karma10 | 50 | phi(1,1) | 0.1652 | 0.1652 | -0.0 | 0.000279 | 0.0003 | 0.951 | 0.950 | 0.950 |
| karma10 | 50 | shape | 19.6417 | 19.6295 | 1.1 | 0.123 | 0.123 | 0.946 | 0.946 | 0.948 |
| karma11 | 10 | alpha | 0.2895 | 0.2763 | 1.1 | 0.133 | 0.017 | 0.808 | 0.842 | 0.801 |
| karma11 | 10 | phi(0,1) | 0.2891 | 0.2796 | 1.3 | 0.0504 | 0.039 | 0.775 | 0.813 | 0.793 |
| karma11 | 10 | phi(1,0) | -0.1641 | -0.1724 | 1.2 | 0.0465 | 0.0488 | 0.784 | 0.812 | 0.780 |
| karma11 | 10 | phi(1,1) | 0.3103 | 0.3251 | -2.1 | 0.0567 | 0.0379 | 0.781 | 0.809 | 0.746 |
| karma11 | 10 | shape | 17.9409 | 18.0321 | -1.6 | 4.09 | 4.04 | 0.914 | 0.920 | 0.910 |
| karma11 | 10 | theta(0,1) | 0.1595 | 0.1599 | -0.1 | 0.0559 | 0.06 | 0.754 | 0.794 | 0.760 |
| karma11 | 10 | theta(1,0) | -0.0910 | -0.0767 | -1.5 | 0.0912 | 0.0823 | 0.740 | 0.779 | 0.733 |
| karma11 | 10 | theta(1,1) | 0.1667 | 0.1542 | 1.5 | 0.0734 | 0.0797 | 0.746 | 0.801 | 0.740 |
| karma11 | 25 | alpha | 0.2860 | 0.2856 | 0.3 | 0.00187 | 0.0018 | 0.914 | 0.916 | 0.914 |
| karma11 | 25 | phi(0,1) | 0.2926 | 0.2937 | -0.6 | 0.00352 | 0.0038 | 0.934 | 0.948 | 0.934 |
| karma11 | 25 | phi(1,0) | -0.1544 | -0.1547 | 0.1 | 0.00495 | 0.0049 | 0.924 | 0.932 | 0.934 |
| karma11 | 25 | phi(1,1) | 0.2712 | 0.2721 | -0.5 | 0.00462 | 0.0048 | 0.852 | 0.868 | 0.855 |
| karma11 | 25 | shape | 17.0323 | 17.0415 | -0.5 | 0.334 | 0.321 | 0.948 | 0.948 | 0.947 |
| karma11 | 25 | theta(0,1) | 0.1569 | 0.1561 | 0.4 | 0.00376 | 0.0041 | 0.942 | 0.951 | 0.930 |
| karma11 | 25 | theta(1,0) | -0.0991 | -0.0966 | -1.1 | 0.00541 | 0.0053 | 0.934 | 0.937 | 0.931 |
| karma11 | 25 | theta(1,1) | 0.2113 | 0.2121 | -0.4 | 0.00618 | 0.0063 | 0.892 | 0.903 | 0.894 |
| karma11 | 50 | alpha | 0.2920 | 0.2921 | -0.1 | 0.00048 | 0.0005 | 0.920 | 0.917 | 0.921 |
| karma11 | 50 | phi(0,1) | 0.2975 | 0.2963 | 1.3 | 0.00086 | 0.0009 | 0.946 | 0.950 | 0.942 |
| karma11 | 50 | phi(1,0) | -0.1517 | -0.1515 | -0.2 | 0.00117 | 0.0012 | 0.944 | 0.945 | 0.944 |
| karma11 | 50 | phi(1,1) | 0.2517 | 0.2523 | -0.7 | 0.00114 | 0.0012 | 0.879 | 0.878 | 0.858 |
| karma11 | 50 | shape | 16.9930 | 16.9797 | 1.6 | 0.0718 | 0.0752 | 0.957 | 0.959 | 0.944 |
| karma11 | 50 | theta(0,1) | 0.1518 | 0.1536 | -1.9 | 0.000869 | 0.0009 | 0.953 | 0.954 | 0.943 |
| karma11 | 50 | theta(1,0) | -0.0986 | -0.0987 | 0.1 | 0.00123 | 0.0013 | 0.932 | 0.938 | 0.942 |
| karma11 | 50 | theta(1,1) | 0.2308 | 0.2300 | 0.8 | 0.00142 | 0.0015 | 0.912 | 0.921 | 0.906 |
