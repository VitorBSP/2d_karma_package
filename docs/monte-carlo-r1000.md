# Monte Carlo results, 1,000 replications per cell

Intermediate-size run of the Monte Carlo study of Section 4.1. This is **not**
the article's configuration of 10,000 replications per cell, which was not
run in the audit environment. Command:

```
Rscript analysis/simulations/run_monte_carlo.R --mode=full --replications=1000 \
  --cores=3 --output=analysis/output/simulations_r1000b
Rscript analysis/simulations/summarize_monte_carlo.R \
  --input=analysis/output/simulations_r1000b
```

Settings: the article's true values (scenario set `paper`), the authors'
simulator (`burn_in = max(p, q)`), master seed 20260101, karma2d
0.1.0.9000 with the objective scaled by the pixel count, R 4.3.3. The run
took 655 s of wall time on 3 cores (238 s for KARMA(1,1) 50 x 50), so the
full study should take about 2 h on similar hardware. Of the 9,000 fits,
none failed and one (KARMA(1,1), 10 x 10) did not converge and is excluded.
Machine-readable results: `docs/monte-carlo-r1000-comparison.csv`.

This replaces an earlier table produced before the optimizer fix described
in `docs/validation-authors-code.md`, Sec. 4.5; that table showed inflated
MSEs on 10 x 10 images caused by one package fit, not by the article.

Columns: mean of the estimates (`mean`) and the article's value (`paper`);
`z` = difference divided by its Monte Carlo standard error (sd / sqrt(1000));
MSE = mean squared error; coverage of nominal 95 % Wald intervals with
standard errors from the expected (`CR exp`), observed (`CR obs`) and mixed
(`CR mix`, the matrix of the authors' code) information, and the article's
coverage (`CR paper`). Monte Carlo standard error of a coverage near 0.95:
0.007.

## Summary

* **Means.** 50 of 54 are within +/-2 Monte Carlo standard errors of the
  article and all 54 within +/-3.
* **MSE.** Between 0.90 and 1.27 times the article's value for every
  parameter whose published MSE is at least 0.001 (smaller ones are printed
  with one significant digit).
* **Coverage.** The article's pattern is reproduced, including the
  under-coverage on 10 x 10 images and of `phi(1,1)` and `theta(1,1)` in the
  ARMA model at every size. The largest difference from the article is
  0.033: `theta(1,1)` of KARMA(0,1) at 25 x 25, 0.915 against 0.948, about
  four Monte Carlo standard errors. A second stream of 500 replications
  gave 0.930, and the authors' optimizer settings give the same coverage as
  the package on the same images, so this looks like simulation noise, but
  only the 10,000-replication run can settle it. Coverage with the mixed
  matrix equals coverage with the expected information to three decimals in
  almost every cell.

## Full table

| Scenario | Size | Parameter | Mean | Paper | z | MSE | MSE paper | CR exp | CR obs | CR mix | CR paper |
|---|---|---|---|---|---|---|---|---|---|---|---|
| karma01 | 10 | alpha | 1.1117 | 1.1134 | -0.7 | 0.00531 | 0.0054 | 0.929 | 0.932 | 0.930 | 0.924 |
| karma01 | 10 | shape | 12.6100 | 12.6818 | -1.8 | 1.58 | 1.62 | 0.940 | 0.940 | 0.938 | 0.941 |
| karma01 | 10 | theta(0,1) | 0.3443 | 0.3459 | -0.5 | 0.0112 | 0.0124 | 0.896 | 0.921 | 0.896 | 0.887 |
| karma01 | 10 | theta(1,0) | 0.2782 | 0.2844 | -1.7 | 0.0127 | 0.0131 | 0.885 | 0.918 | 0.885 | 0.885 |
| karma01 | 10 | theta(1,1) | 0.1739 | 0.1810 | -1.7 | 0.0166 | 0.0163 | 0.891 | 0.916 | 0.891 | 0.882 |
| karma01 | 25 | alpha | 1.1103 | 1.1113 | -1.2 | 0.000809 | 0.0008 | 0.938 | 0.938 | 0.937 | 0.946 |
| karma01 | 25 | shape | 12.3536 | 12.3501 | 0.2 | 0.196 | 0.195 | 0.941 | 0.941 | 0.941 | 0.941 |
| karma01 | 25 | theta(0,1) | 0.3437 | 0.3460 | -2.1 | 0.00121 | 0.0011 | 0.927 | 0.932 | 0.927 | 0.943 |
| karma01 | 25 | theta(1,0) | 0.2817 | 0.2820 | -0.3 | 0.00119 | 0.0012 | 0.942 | 0.952 | 0.942 | 0.944 |
| karma01 | 25 | theta(1,1) | 0.1846 | 0.1826 | 1.6 | 0.0015 | 0.0013 | 0.915 | 0.916 | 0.915 | 0.948 |
| karma01 | 50 | alpha | 1.1107 | 1.1112 | -1.1 | 0.00019 | 0.0002 | 0.951 | 0.951 | 0.951 | 0.948 |
| karma01 | 50 | shape | 12.3524 | 12.3483 | 0.6 | 0.0428 | 0.0472 | 0.946 | 0.945 | 0.945 | 0.942 |
| karma01 | 50 | theta(0,1) | 0.3464 | 0.3459 | 1.1 | 0.000248 | 0.0002 | 0.953 | 0.954 | 0.953 | 0.948 |
| karma01 | 50 | theta(1,0) | 0.2820 | 0.2821 | -0.3 | 0.000268 | 0.0003 | 0.944 | 0.948 | 0.944 | 0.949 |
| karma01 | 50 | theta(1,1) | 0.1839 | 0.1843 | -0.9 | 0.00027 | 0.0003 | 0.952 | 0.951 | 0.952 | 0.948 |
| karma10 | 10 | alpha | 0.7140 | 0.7073 | 1.1 | 0.0407 | 0.0376 | 0.914 | 0.935 | 0.914 | 0.929 |
| karma10 | 10 | phi(0,1) | 0.2693 | 0.2686 | 0.2 | 0.00945 | 0.0093 | 0.911 | 0.921 | 0.911 | 0.922 |
| karma10 | 10 | phi(1,0) | 0.1755 | 0.1761 | -0.2 | 0.00934 | 0.0095 | 0.928 | 0.945 | 0.928 | 0.926 |
| karma10 | 10 | phi(1,1) | 0.1617 | 0.1633 | -0.5 | 0.0109 | 0.0102 | 0.926 | 0.937 | 0.926 | 0.927 |
| karma10 | 10 | shape | 20.5971 | 20.4864 | 1.7 | 5.19 | 4.83 | 0.931 | 0.935 | 0.936 | 0.936 |
| karma10 | 25 | alpha | 0.6709 | 0.6706 | 0.2 | 0.00477 | 0.0049 | 0.954 | 0.954 | 0.954 | 0.949 |
| karma10 | 25 | phi(0,1) | 0.2785 | 0.2796 | -1.0 | 0.00114 | 0.0011 | 0.936 | 0.944 | 0.936 | 0.948 |
| karma10 | 25 | phi(1,0) | 0.1885 | 0.1877 | 0.7 | 0.0011 | 0.0011 | 0.943 | 0.943 | 0.943 | 0.945 |
| karma10 | 25 | phi(1,1) | 0.1638 | 0.1639 | -0.1 | 0.0013 | 0.0012 | 0.945 | 0.943 | 0.945 | 0.945 |
| karma10 | 25 | shape | 19.7535 | 19.7205 | 1.5 | 0.535 | 0.529 | 0.942 | 0.943 | 0.943 | 0.946 |
| karma10 | 50 | alpha | 0.6638 | 0.6638 | 0.0 | 0.00129 | 0.0013 | 0.953 | 0.953 | 0.953 | 0.948 |
| karma10 | 50 | phi(0,1) | 0.2818 | 0.2810 | 1.6 | 0.000255 | 0.0003 | 0.946 | 0.946 | 0.946 | 0.950 |
| karma10 | 50 | phi(1,0) | 0.1890 | 0.1890 | -0.1 | 0.00025 | 0.0003 | 0.963 | 0.965 | 0.963 | 0.950 |
| karma10 | 50 | phi(1,1) | 0.1645 | 0.1652 | -1.3 | 0.000292 | 0.0003 | 0.942 | 0.942 | 0.942 | 0.950 |
| karma10 | 50 | shape | 19.6389 | 19.6295 | 0.9 | 0.124 | 0.123 | 0.949 | 0.949 | 0.950 | 0.948 |
| karma11 | 10 | alpha | 0.2848 | 0.2763 | 1.8 | 0.0217 | 0.017 | 0.805 | 0.839 | 0.805 | 0.801 |
| karma11 | 10 | phi(0,1) | 0.2766 | 0.2796 | -0.5 | 0.0431 | 0.039 | 0.778 | 0.808 | 0.778 | 0.793 |
| karma11 | 10 | phi(1,0) | -0.1910 | -0.1724 | -2.6 | 0.0511 | 0.0488 | 0.780 | 0.798 | 0.780 | 0.780 |
| karma11 | 10 | phi(1,1) | 0.3309 | 0.3251 | 1.1 | 0.0408 | 0.0379 | 0.728 | 0.769 | 0.728 | 0.746 |
| karma11 | 10 | shape | 17.9801 | 18.0321 | -1.0 | 3.89 | 4.04 | 0.906 | 0.914 | 0.910 | 0.910 |
| karma11 | 10 | theta(0,1) | 0.1657 | 0.1599 | 0.7 | 0.0621 | 0.06 | 0.747 | 0.785 | 0.747 | 0.760 |
| karma11 | 10 | theta(1,0) | -0.0511 | -0.0767 | 2.7 | 0.0907 | 0.0823 | 0.724 | 0.772 | 0.724 | 0.733 |
| karma11 | 10 | theta(1,1) | 0.1562 | 0.1542 | 0.2 | 0.0761 | 0.0797 | 0.732 | 0.778 | 0.732 | 0.740 |
| karma11 | 25 | alpha | 0.2864 | 0.2856 | 0.6 | 0.0019 | 0.0018 | 0.911 | 0.924 | 0.911 | 0.914 |
| karma11 | 25 | phi(0,1) | 0.2949 | 0.2937 | 0.6 | 0.00393 | 0.0038 | 0.920 | 0.936 | 0.920 | 0.934 |
| karma11 | 25 | phi(1,0) | -0.1574 | -0.1547 | -1.2 | 0.00545 | 0.0049 | 0.912 | 0.930 | 0.912 | 0.934 |
| karma11 | 25 | phi(1,1) | 0.2710 | 0.2721 | -0.7 | 0.00455 | 0.0048 | 0.865 | 0.875 | 0.865 | 0.855 |
| karma11 | 25 | shape | 17.0377 | 17.0415 | -0.2 | 0.318 | 0.321 | 0.939 | 0.943 | 0.943 | 0.947 |
| karma11 | 25 | theta(0,1) | 0.1545 | 0.1561 | -0.8 | 0.00402 | 0.0041 | 0.927 | 0.939 | 0.927 | 0.930 |
| karma11 | 25 | theta(1,0) | -0.0944 | -0.0966 | 0.9 | 0.00595 | 0.0053 | 0.918 | 0.927 | 0.918 | 0.931 |
| karma11 | 25 | theta(1,1) | 0.2130 | 0.2121 | 0.4 | 0.00679 | 0.0063 | 0.880 | 0.895 | 0.880 | 0.894 |
| karma11 | 50 | alpha | 0.2920 | 0.2921 | -0.2 | 0.000439 | 0.0005 | 0.929 | 0.927 | 0.929 | 0.921 |
| karma11 | 50 | phi(0,1) | 0.2970 | 0.2963 | 0.7 | 0.000886 | 0.0009 | 0.945 | 0.947 | 0.945 | 0.942 |
| karma11 | 50 | phi(1,0) | -0.1512 | -0.1515 | 0.3 | 0.00113 | 0.0012 | 0.952 | 0.955 | 0.952 | 0.944 |
| karma11 | 50 | phi(1,1) | 0.2520 | 0.2523 | -0.4 | 0.00121 | 0.0012 | 0.867 | 0.866 | 0.867 | 0.858 |
| karma11 | 50 | shape | 16.9980 | 16.9797 | 2.2 | 0.0686 | 0.0752 | 0.955 | 0.955 | 0.955 | 0.944 |
| karma11 | 50 | theta(0,1) | 0.1527 | 0.1536 | -1.0 | 0.000892 | 0.0009 | 0.953 | 0.952 | 0.953 | 0.943 |
| karma11 | 50 | theta(1,0) | -0.0990 | -0.0987 | -0.3 | 0.00117 | 0.0013 | 0.955 | 0.955 | 0.955 | 0.942 |
| karma11 | 50 | theta(1,1) | 0.2306 | 0.2300 | 0.6 | 0.00144 | 0.0015 | 0.912 | 0.916 | 0.912 | 0.906 |
