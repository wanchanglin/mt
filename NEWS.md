
# mt 2.0-1.1:

- `pca.outlier`: add prefix to ellipse to avoid conflict with package car.

# mt 2.0-1.2:

- `corrgram.circle` and `corrgram.ellipse`: move 'scales' into argument list.

# mt 2.0-1.3:  

- `cl.auc`: adjust auc: `if (auc < 0.5) auc <- 1 - auc`

# mt 2.0-1.4: (14-05-2012)   

- `accest`: Fix bugs and add more output.

# mt 2.0-1.5: (16-05-2012)   

- `cl.roc`: change output and fix bugs.

# mt 2.0-1.6: (19-05-2012)   

- `cl.roc`: add output of positive level.
- `cl.perf`: add confusion matrix, Kappa and positive outputs.

# mt 2.0-1.7: (21-05-2012)   

- `plsc`: minor changes for returned component pls.out so that R2 can be
  applied.
- `plslda`: minor changes for returned component pls.out so that R2 can be
	applied.

# mt 2.0-1.8: (22-05-2012)   

- `classifier`: bugs fix

# mt 2.0-1.9: (22-05-2013)

- `df.summ`: add dots argument and modify an user-defined function.

# mt 2.0-1.10: (09-01-2014)

- `stats.vec` and `stats.mat`: summarising two groups stats, such as
	p-values.
- `fs.cl.2`, `perf` and `perf.aam`: internal functions (21-01-2014).

# mt 2.0-1.11: (22-01-2014)

- re-write `shrink.list`
- add `proprec.auto` as hidden function and tidy up

# mt 2.0-1.12: (07-08-2014)

- fix bug in `pval.reject`: add na.rm = TRUE 
- change in `.grpplot`: remove xlab and ylab inside.

# mt 2.0-1.13: (10-06-2015)

- two functions `mdsplot` and `mds.plot.wrap` for MDS plotting.

# mt 2.0-1.14: (16-07-2015)

- Re-write two functions for ellipse: `panel.elli` and `panel.elli.1`
- Provide `panel.outl` for labelling data points outside ellipse. They can be
	considered as outliers.

# mt 2.0-1.15: (19-11-2015)

- `.path.package` is defunct so use `path.package` instead in `csv2xls`.
- Use `pls:::coef.mvr` replace `coef.mvr` since it hides in the new version
	of 'pls'.

# mt 2.0-1.16: (30-11-2015)

- Add log2 fold-change for `stats.mat` and `stats.vec`.
- Add adjusted p-values for `stats.mat` only.
- Provide documents for `pca.outlier` and `pca.outlier.1`. The former is
lattice version and the later is the basic graphics version.

# mt 2.0-1.17: (26-01-2016)

- Remove returned values of "direction" in `stats.vec` and `stats.mat`. So 
	the results are numeric, not character.

# mt 2.0-1.18: (06-11-2021)

- Replace `melt` and `ddply` with base R functions. So no need to load two
  packages `reshape` and `plyr`.
- Move out some un-documented functions. 
- Change CHANGELOG as NEWS.md
- Move packages from Depends to Imports in DESCRIPTION
- Spell check on R scripts using RStudio


# mt 2.0-1.19: (31-01-2022)

- no changes but test against new version of randomForest 4.7.1 as requested.