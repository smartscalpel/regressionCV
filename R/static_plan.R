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
  smpl_splited_fm=target(smpl_split_fm(fm,split=0.75),transform = map(fm)),# split feature matrix into train/test parts by patientid
  normalized_fm=target(normalize(fm=smpl_splited_fm,normtype),transform = cross(smpl_splited_fm,normtype=!!c('Pareto'))),#'None','Autoscaling','Pareto'))),# scale feature matrix \cite{vandenBerg:2006hm}
  filter_fm=target(feature_filter(fm=normalized_fm,ftype),transform = cross(normalized_fm,ftype=!!c('None'))),#'None','ZVar','Corr'))),# reduce feature space 
  xgb_cv10=target(train_model(fm=filter_fm,modeltype='xgb'),transform = map(filter_fm)),#,trigger = trigger(condition =FALSE)),# train regression model with CV10
  test_xgb=target(test_model(xgb_cv10),transform = map(xgb_cv10)),
  eval_xgb=target(eval_model(tst=test_xgb),transform = map(test_xgb)),
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
  plot_xgb_boxSmpl=target(plot_train_smpl_box(fm=test_xgb),transform = map(test_xgb)),
  save_plot_xgb_boxSmpl = target(ggsave(
    filename = file_out(sprintf('peak2019.diag_%d.res_%d.mode_%d.mz_%d.ftype_%s.norm_%s.cv10.%s.plot_box.test.smpl.pdf',
                                diag,res,mode,mz,normtype,ftype,"xgb")),
    plot = plot_xgb_boxSmpl,
    width = 8,
    height = 8
  ),transform = map(plot_xgb_boxSmpl)),
  #### TCP plots ####
  plot_xgb_tcp_boxSmpl=target(plot_tcp_smpl_box(fm=test_xgb,theme=smpl_box_theme,palette='Dark2'),transform = map(test_xgb)),
  save_plot_xgb_tcp_boxSmpl = target(ggsave(
    filename = file_out(sprintf('peak2019.diag_%d.res_%d.mode_%d.mz_%d.ftype_%s.norm_%s.cv10.%s.plot_box.tcp.smpl.pdf',
                                diag,res,mode,mz,normtype,ftype,"xgb")),
    plot = plot_xgb_tcp_boxSmpl,
    width = 8,
    height = 8
  ),transform = map(plot_xgb_tcp_boxSmpl)),
  plot_xgb_tcp_box=target(plot_tcp_box(fm=test_xgb),transform = map(test_xgb)),
  save_plot_xgb_tcp_box = target(ggsave(
    filename = file_out(sprintf('peak2019.diag_%d.res_%d.mode_%d.mz_%d.ftype_%s.norm_%s.cv10.%s.plot_tcp_box.pdf',
                                diag,res,mode,mz,normtype,ftype,"xgb")),
    plot = plot_xgb_tcp_box,
    width = 8,
    height = 8
  ),transform = map(plot_xgb_tcp_box)),
  plot_median_xgb_tcp_box=target(plot_median_spectrum_tcp_box(fm=test_xgb,theme=theme_light(base_size = 32)),transform = map(test_xgb)),
  save_median_plot_xgb_tcp_box = target(ggsave(
    filename = file_out(sprintf('peak2019.diag_%d.res_%d.mode_%d.mz_%d.ftype_%s.norm_%s.cv10.%s.plot_median_box.tcp.pdf',
                                diag,res,mode,mz,normtype,ftype,"xgb")),
    plot = plot_median_xgb_tcp_box,
    width = 8,
    height = 8
  ),transform = map(plot_median_xgb_tcp_box)),
  plot_mean_xgb_tcp_box=target(plot_mean_spectrum_tcp_box(fm=test_xgb,theme=theme_light(base_size = 32)),transform = map(test_xgb)),
  save_mean_plot_xgb_tcp_box = target(ggsave(
    filename = file_out(sprintf('peak2019.diag_%d.res_%d.mode_%d.mz_%d.ftype_%s.norm_%s.cv10.%s.plot_mean_box.tcp.pdf',
                                diag,res,mode,mz,normtype,ftype,"xgb")),
    plot = plot_mean_xgb_tcp_box,
    width = 8,
    height = 8
  ),transform = map(plot_mean_xgb_tcp_box)),
  
  #### Spectrum median plots ####
  plot_median_xgb_boxT=target(plot_median_spectrum_test_box(fm=test_xgb),transform = map(test_xgb)),
  plot_median_xgb_box=target(plot_median_spectrum_train_box(fm=test_xgb),transform = map(test_xgb)),
  save_median_plot_xgb_boxT = target(ggsave(
    filename = file_out(sprintf('peak2019.diag_%d.res_%d.mode_%d.mz_%d.ftype_%s.norm_%s.cv10.%s.plot_median_box.pdf',
                                diag,res,mode,mz,normtype,ftype,"xgb")),
    plot = plot_median_xgb_boxT,
    width = 8,
    height = 8
  ),transform = map(plot_median_xgb_boxT)),
  save_median_plot_xgb_box = target(ggsave(
    filename = file_out(sprintf('peak2019.diag_%d.res_%d.mode_%d.mz_%d.ftype_%s.norm_%s.cv10.%s.plot_median_box.train.pdf',
                                diag,res,mode,mz,normtype,ftype,"xgb")),
    plot = plot_median_xgb_box,
    width = 8,
    height = 8
  ),transform = map(plot_median_xgb_box)),
  #### Spectrum mean plots ####
  plot_mean_xgb_boxT=target(plot_mean_spectrum_test_box(fm=test_xgb),transform = map(test_xgb)),
  plot_mean_xgb_box=target(plot_mean_spectrum_train_box(fm=test_xgb),transform = map(test_xgb)),
  save_mean_plot_xgb_boxT = target(ggsave(
    filename = file_out(sprintf('peak2019.diag_%d.res_%d.mode_%d.mz_%d.ftype_%s.norm_%s.cv10.%s.plot_mean_box.pdf',
                                diag,res,mode,mz,normtype,ftype,"xgb")),
    plot = plot_mean_xgb_boxT,
    width = 8,
    height = 8
  ),transform = map(plot_mean_xgb_boxT)),
  save_mean_plot_xgb_box = target(ggsave(
    filename = file_out(sprintf('peak2019.diag_%d.res_%d.mode_%d.mz_%d.ftype_%s.norm_%s.cv10.%s.plot_mean_box.train.pdf',
                                diag,res,mode,mz,normtype,ftype,"xgb")),
    plot = plot_mean_xgb_box,
    width = 8,
    height = 8
  ),transform = map(plot_mean_xgb_box)),
  #### IML ####  
  get_xgb_imp=target(xgb_importance(mod=xgb_cv10),transform = map(xgb_cv10)),
  plot_xgb_imp=target(xgb_plot_importance(imp=get_xgb_imp),transform = map(get_xgb_imp)),
  save_plot_xgb_imp = target(ggsave(
    filename = file_out(sprintf('importance.diag_%d.res_%d.mode_%d.mz_%d.ftype_%s.norm_%s.cv10.%s.pdf',
                                diag,res,mode,mz,normtype,ftype,"xgb")),
    plot = plot_xgb_imp,
    width = 8,
    height = 24
  ),transform = map(plot_xgb_imp)),
  plot_shap=target(get_shap_plot(imp=get_xgb_imp),transform = map(get_xgb_imp)),
  save_plot_shap = target(ggsave(
    filename = file_out(sprintf('shap_summary.diag_%d.res_%d.mode_%d.mz_%d.ftype_%s.norm_%s.cv10.%s.png',
                                diag,res,mode,mz,normtype,ftype,"xgb")),
    plot = plot_shap,
    width = 8,
    height = 8
  ),transform = map(plot_shap)),
  get_xgb_shap=target(get_shap_values(imp=get_xgb_imp),transform = map(get_xgb_imp)),
  plot_shap_xgb=target(get_shap_dep_plot(shap_values=get_xgb_shap),transform = map(get_xgb_shap)),
  save_plot_shap_xgb = target(ggsave(
    filename = file_out(sprintf('shap_dep_plot.diag_%d.res_%d.mode_%d.mz_%d.ftype_%s.norm_%s.cv10.%s.png',
                                diag,res,mode,mz,normtype,ftype,"xgb")),
    plot = plot_shap_xgb,
    width = 14,
    height = 14
  ),transform = map(plot_shap_xgb)),
  shap_plot_data=target(get_shap_plot_data(shap_values=get_xgb_shap),transform = map(get_xgb_shap)),
  plot_shap_force=target(get_shap_force_plot(plot_data=shap_plot_data),transform = map(shap_plot_data)),
  save_plot_shap_force = target(ggsave(
    filename = file_out(sprintf('shap_force_plot.diag_%d.res_%d.mode_%d.mz_%d.ftype_%s.norm_%s.cv10.%s.png',
                                diag,res,mode,mz,normtype,ftype,"xgb")),
    plot = plot_shap_force,
    width = 14,
    height = 14
  ),transform = map(plot_shap_force)),
  plot_shap_fgroup=target(get_shap_force_group_plot(plot_data=shap_plot_data),transform = map(shap_plot_data)),
  save_plot_shap_fgroup = target(ggsave(
    filename = file_out(sprintf('shap_fgroup_plot.diag_%d.res_%d.mode_%d.mz_%d.ftype_%s.norm_%s.cv10.%s.png',
                                diag,res,mode,mz,normtype,ftype,"xgb")),
    plot = plot_shap_fgroup,
    width = 14,
    height = 14
  ),transform = map(plot_shap_fgroup)),
  #### Reduced model analysis ####
  reduced_fm=target(get_reduced_fm(imp=get_xgb_imp),transform = map(get_xgb_imp)),
  reduced_xgb_cv10=target(train_model(fm=reduced_fm,modeltype='xgb'),transform = map(reduced_fm)),#,trigger = trigger(condition =FALSE)),# train regression model with CV10
  reduced_test_xgb=target(test_model(reduced_xgb_cv10),transform = map(reduced_xgb_cv10)),
  plot_reduced_xgb_box=target(plot_test_box(fm=reduced_test_xgb),transform = map(reduced_test_xgb)),
  plot_reduced_xgb_boxT=target(plot_train_box(fm=reduced_test_xgb),transform = map(reduced_test_xgb)),
  save_plot_reduced_xgb_box = target(ggsave(
    filename = file_out(sprintf('peak2019.reduced.diag_%d.res_%d.mode_%d.mz_%d.ftype_%s.norm_%s.cv10.%s.plot_reduced_box.pdf',
                                diag,res,mode,mz,normtype,ftype,"xgb")),
    plot = plot_reduced_xgb_box,
    width = 8,
    height = 8
  ),transform = map(plot_reduced_xgb_box)),
  save_plot_reduced_xgb_boxT = target(ggsave(
    filename = file_out(sprintf('peak2019.reduced.diag_%d.res_%d.mode_%d.mz_%d.ftype_%s.norm_%s.cv10.%s.plot_reduced_box.train.pdf',
                                diag,res,mode,mz,normtype,ftype,"xgb")),
    plot = plot_reduced_xgb_boxT,
    width = 8,
    height = 8
  ),transform = map(plot_reduced_xgb_boxT)),
  plot_reduced_xgb_boxSmpl=target(plot_train_smpl_box(fm=reduced_test_xgb),transform = map(reduced_test_xgb)),
  save_plot_reduced_xgb_boxSmpl = target(ggsave(
    filename = file_out(sprintf('peak2019.reduced.diag_%d.res_%d.mode_%d.mz_%d.ftype_%s.norm_%s.cv10.%s.plot_reduced_box.test.smpl.pdf',
                                diag,res,mode,mz,normtype,ftype,"xgb")),
    plot = plot_reduced_xgb_boxSmpl,
    width = 8,
    height = 8
  ),transform = map(plot_reduced_xgb_boxSmpl)),
  plot_reduced_median_xgb_boxT=target(plot_median_spectrum_test_box(fm=reduced_test_xgb),transform = map(reduced_test_xgb)),
  plot_reduced_median_xgb_box=target(plot_median_spectrum_train_box(fm=reduced_test_xgb),transform = map(reduced_test_xgb)),
  save_median_plot_reduced_xgb_boxT = target(ggsave(
    filename = file_out(sprintf('peak2019.reduced.diag_%d.res_%d.mode_%d.mz_%d.ftype_%s.norm_%s.cv10.%s.plot_reduced_median_box.pdf',
                                diag,res,mode,mz,normtype,ftype,"xgb")),
    plot = plot_reduced_median_xgb_boxT,
    width = 8,
    height = 8
  ),transform = map(plot_reduced_median_xgb_boxT)),
  save_median_plot_reduced_xgb_box = target(ggsave(
    filename = file_out(sprintf('peak2019.reduced.diag_%d.res_%d.mode_%d.mz_%d.ftype_%s.norm_%s.cv10.%s.plot_reduced_median_box.train.pdf',
                                diag,res,mode,mz,normtype,ftype,"xgb")),
    plot = plot_reduced_median_xgb_box,
    width = 8,
    height = 8
  ),transform = map(plot_reduced_median_xgb_box)),
  plot_reduced_mean_xgb_boxT=target(plot_mean_spectrum_test_box(fm=reduced_test_xgb),transform = map(reduced_test_xgb)),
  plot_reduced_mean_xgb_box=target(plot_mean_spectrum_train_box(fm=reduced_test_xgb),transform = map(reduced_test_xgb)),
  save_mean_plot_reduced_xgb_boxT = target(ggsave(
    filename = file_out(sprintf('peak2019.reduced.diag_%d.res_%d.mode_%d.mz_%d.ftype_%s.norm_%s.cv10.%s.plot_reduced_mean_box.pdf',
                                diag,res,mode,mz,normtype,ftype,"xgb")),
    plot = plot_reduced_mean_xgb_boxT,
    width = 8,
    height = 8
  ),transform = map(plot_reduced_mean_xgb_boxT)),
  save_mean_plot_reduced_xgb_box = target(ggsave(
    filename = file_out(sprintf('peak2019.reduced.diag_%d.res_%d.mode_%d.mz_%d.ftype_%s.norm_%s.cv10.%s.plot_reduced_mean_box.train.pdf',
                                diag,res,mode,mz,normtype,ftype,"xgb")),
    plot = plot_reduced_mean_xgb_box,
    width = 8,
    height = 8
  ),transform = map(plot_reduced_mean_xgb_box))  
  )
