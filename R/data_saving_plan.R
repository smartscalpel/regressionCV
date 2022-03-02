plan <- drake_plan(
  trace = TRUE,
  files = target(
    load_dataset(res,mode,mz,diag),
    # Define an analysis target for each combination of
    # tuning_setting, mean_value, and model_function.
    transform = cross(
      res=!!c(2),
      mode=!!c(2),#(1,2),
      mz=!!c(2),#(2,1),
      diag=!!c(6),#(3,6),
      .id = c(diag,res,mode,mz),.tag_out=dataset
    )
  ),
  fm=target(prepare_feature_matrix(files,monoisotopic=FALSE),transform=map(files,.id = c(diag,res,mode,mz),.tag_out=dataset)),
  patDF = target(get_pat_df(fm),
                 transform = map(fm)),
  specDF = target(get_spec_df(fm),
                  transform = map(fm)),
  dimPat=target(get_dim(patDF),transform = map(patDF)),
  dimSpec=target(get_dim(specDF),transform = map(specDF)),
  patStat= target(bind_rows(dimPat),transform = combine(
    dimPat)),
  specStat= target(bind_rows(dimSpec),transform = combine(
    dimSpec)),
  wrtPatStat=write.csv(patStat,file=file_out('patStat.csv')),
  wrtSpecStat=write.csv(specStat,file=file_out('specStat.csv')),
  wrtSpec=target(write.csv(fm,
              file = file_out(sprintf('peak2019.reduced.diag_%d.res_%d.mode_%d.mz_%d.csv',
                                diag,res,mode,mz))),transform = map(fm))
)