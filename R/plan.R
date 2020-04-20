plan <- drake_plan(
  ms_setup=target(c(res=r,mode=m,mz=z),
                  transform = cross(r=!!c(2),
                                    m=!!c(1,2),
                                    z=!!c(1,2))),
  fm = target(
    load_dataset(ms_setup,diag,expt),
    # Define an analysis target for each combination of
    # tuning_setting, mean_value, and model_function.
    transform = cross(
      ms_setup,diag=!!c(3,6),expt=!!c(1,2),
      .id = c(diag,expt,ms_setup),.tag_out=dataset
    )
  ),
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
  wrtPatStat=write.csv(patStat,file='patStat.csv'),
  wrtSpecStat=write.csv(specStat,file='specStat.csv')
  # model = target(
  #   get_model(ms_setup,diag,expt,method,idx),
  #   # Define an analysis target for each combination of
  #   # tuning_setting, mean_value, and model_function.
  #   transform = cross(
  #     ms_setup,diag=!!c(3,6,15),expt=!!c(1,2),
  #     method = c("RF", "XGB"),
  #     idx = !!c("pat","spec"),
  #     .id = c(diag,expt,ms_setup,method,idx),.tag_out=dataset
  #   )
  # ),
  # predict=target(
  #   predict_dataset(model,ms_setup,ddiag,dexpt),
  #   transform = cross(
  #     model,ddiag=!!c(3,6,15),dexpt=!!c(1,2),
  #     .id = c(diag,expt,ms_setup,method,idx,ddiag,dexpt),.tag_out=preds
  #   )
  # ),
  # matrix=target(
  #   bind_rows(preds,.id=preds),
  #   transform = combine(preds,.by=c(method,idx,ms_setup))
  # )
)
