plan <- drake_plan(
  trace = TRUE,
  files = target(
    load_dataset(res,mode,mz,diag),
    # Define an analysis target for each combination of
    # tuning_setting, mean_value, and model_function.
    transform = cross(
      res=!!c(2),
      mode=!!c(1,2),
      mz=!!c(2,1),
      diag=!!c(3,6),
      .id = c(diag,res,mode,mz),.tag_out=dataset
    )
  ),
  fm=target(prepare_feature_matrix(files),transform=map(files,.id = c(diag,res,mode,mz),.tag_out=dataset)),
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
  smpl_splited_fm=target(smpl_split_fm(fm,split=0.75),transform = map(fm)),# split feature matrix into train/test parts by patientid
  normalized_fm=target(normalize(fm=smpl_splited_fm,normtype),transform = cross(smpl_splited_fm,normtype=!!c('None'))),#,'Autoscaling','Pareto'))),# scale feature matrix \cite{vandenBerg:2006hm}
  filter_fm=target(feature_filter(fm=normalized_fm,ftype),transform = cross(normalized_fm,ftype=!!c('None','ZVar','Corr'))),# reduce feature space 
  # splited_fm=target(transform = map(normalized_fm)),# split feature matrix into train/test parts by spectrumid with respect to percentage
  # pca=target(transform = cross(fm=normalized_fm,color=!!c('diagnosis','spectrumid','patientid'))),# make PCA plots for transformed feature matrices
  # umap=target(transform = cross(fm=normalized_fm,color=!!c('diagnosis','spectrumid','patientid'))),# make umap plots for transformed feature matrices
  rf_cv10=target(train_model(fm=filter_fm,modeltype='rf'),transform = map(filter_fm)),# train regression model with CV10
  test_rf=target(test_model(rf_cv10),transform = map(rf_cv10)),
  plot_rf_point=target(plot_test_point(fm=test_rf),transform = map(test_rf)),
  plot_rf_box=target(plot_test_box(fm=test_rf),transform = map(test_rf)),
  plot_rf_pointT=target(plot_train_point(fm=test_rf),transform = map(test_rf)),
  plot_rf_boxT=target(plot_train_box(fm=test_rf),transform = map(test_rf)),
  xgb_cv10=target(train_model(fm=filter_fm,modeltype='xgb'),transform = map(filter_fm)),# train regression model with CV10
  test_xgb=target(test_model(xgb_cv10),transform = map(xgb_cv10)),
  plot_xgb_point=target(plot_test_point(fm=test_xgb),transform = map(test_xgb)),
  plot_xgb_box=target(plot_test_box(fm=test_xgb),transform = map(test_xgb)),
  plot_xgb_pointT=target(plot_train_point(fm=test_xgb),transform = map(test_xgb)),
  plot_xgb_boxT=target(plot_train_box(fm=test_xgb),transform = map(test_xgb)),
  save_plot_xgb_point = target(ggsave(
    filename = file_out(sprintf('peak2019.diag_%d.res_%d.mode_%d.mz_%d.ftype_%s.norm_%s.cv10.%s.plot_point.pdf',
                                diag,res,mode,mz,normtype,ftype,"xgb")),
    plot = plot_xgb_point,
    width = 8,
    height = 8
  ),transform = map(plot_xgb_point)),
  save_plot_xgb_pointT = target(ggsave(
    filename = file_out(sprintf('peak2019.diag_%d.res_%d.mode_%d.mz_%d.ftype_%s.norm_%s.cv10.%s.plot_point.train.pdf',
                                diag,res,mode,mz,normtype,ftype,"xgb")),
    plot = plot_xgb_pointT,
    width = 8,
    height = 8
  ),transform = map(plot_xgb_pointT)),
  save_plot_rf_point = target(ggsave(
    filename = file_out(sprintf('peak2019.diag_%d.res_%d.mode_%d.mz_%d.ftype_%s.norm_%s.cv10.%s.plot_point.pdf',
                                diag,res,mode,mz,normtype,ftype,"rf")),
    plot = plot_rf_point,
    width = 8,
    height = 8
  ),transform = map(plot_rf_point)),
  save_plot_rf_pointT = target(ggsave(
    filename = file_out(sprintf('peak2019.diag_%d.res_%d.mode_%d.mz_%d.ftype_%s.norm_%s.cv10.%s.plot_point.train.pdf',
                                diag,res,mode,mz,normtype,ftype,"rf")),
    plot = plot_rf_pointT,
    width = 8,
    height = 8
  ),transform = map(plot_rf_pointT)),
  save_plot_xgb_box = target(ggsave(
    filename = file_out(sprintf('peak2019.diag_%d.res_%d.mode_%d.mz_%d.ftype_%s.norm_%s.cv10.%s.plot_box.pdf',
                                diag,res,mode,mz,normtype,ftype,"xgb")),
    plot = plot_xgb_box,
    width = 8,
    height = 8
  ),transform = map(plot_xgb_box)),
  save_plot_xgb_boxT = target(ggsave(
    filename = file_out(sprintf('peak2019.diag_%d.res_%d.mode_%d.mz_%d.ftype_%s.norm_%s.cv10.%s.plot_box.train.pdf',
                                diag,res,mode,mz,normtype,ftype,"xgb")),
    plot = plot_xgb_boxT,
    width = 8,
    height = 8
  ),transform = map(plot_xgb_boxT)),
  save_plot_rf_box = target(ggsave(
    filename = file_out(sprintf('peak2019.diag_%d.res_%d.mode_%d.mz_%d.ftype_%s.norm_%s.cv10.%s.plot_box.pdf',
                                diag,res,mode,mz,normtype,ftype,"rf")),
    plot = plot_rf_box,
    width = 8,
    height = 8
  ),transform = map(plot_rf_box)),
  save_plot_rf_boxT = target(ggsave(
    filename = file_out(sprintf('peak2019.diag_%d.res_%d.mode_%d.mz_%d.ftype_%s.norm_%s.cv10.%s.plot_box.train.pdf',
                                diag,res,mode,mz,normtype,ftype,"rf")),
    plot = plot_rf_boxT,
    width = 8,
    height = 8
  ),transform = map(plot_rf_boxT)),
  get_xgb_imp=target(xgb_importance(mod=xgb_cv10),transform = map(xgb_cv10)),
  plot_xgb_imp=target(xgb_plot_importance(imp=get_xgb_imp),transform = map(get_xgb_imp)),
  save_plot_xgb_imp = target(ggsave(
    filename = file_out(sprintf('importance.diag_%d.res_%d.mode_%d.mz_%d.ftype_%s.norm_%s.cv10.%s.plot_point.pdf',
                                diag,res,mode,mz,normtype,ftype,"xgb")),
    plot = plot_xgb_imp,
    width = 8,
    height = 24
  ),transform = map(plot_xgb_imp)),
  plot_shap=target(get_shap_plot(imp=get_xgb_imp),transform = map(get_xgb_imp)),
  save_plot_shap = target(ggsave(
    filename = file_out(sprintf('shap_summary.diag_%d.res_%d.mode_%d.mz_%d.ftype_%s.norm_%s.cv10.%s.plot_point.png',
                                diag,res,mode,mz,normtype,ftype,"xgb")),
    plot = plot_shap,
    width = 8,
    height = 8
  ),transform = map(plot_shap)),
  plot_shap_xgb=target(get_shap_dep_plot(imp=get_xgb_imp),transform = map(get_xgb_imp)),
  save_plot_shap_xgb = target(ggsave(
    filename = file_out(sprintf('shap_dep_plot.diag_%d.res_%d.mode_%d.mz_%d.ftype_%s.norm_%s.cv10.%s.plot_point.png',
                                diag,res,mode,mz,normtype,ftype,"xgb")),
    plot = plot_shap_xgb,
    width = 14,
    height = 14
  ),transform = map(plot_shap_xgb))  #,
  #report = knit(knitr_in("report.Rmd"), file_out("report.md"), quiet = TRUE)
  
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
