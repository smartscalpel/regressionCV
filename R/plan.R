peaks=get_peaks(dpath)#[1:12]
plan <- drake_plan(
  trace = TRUE,
  fm=target(prepare_feature_matrix(file_in(!!peak),norm_shift=10),transform = map(peak=!!peaks),fname=!!peak),# convert peaks to feature matrix
  diag_stat=target(fm %>% select(diagnosis) %>% unique %>% dim() %>% (function(.x)data.frame(x=.x[1],y=.x[2])),transform = map(fm),fname=!!peak),
  results = target(
    dplyr::bind_rows(diag_stat),
    transform = combine(diag_stat)
  ),
  smpl_splited_fm=target(smpl_split_fm(fm),transform = map(fm)),# split feature matrix into train/test parts by patientid
  normalized_fm=target(normalize(fm=smpl_splited_fm,normtype),transform = cross(smpl_splited_fm,normtype=!!normtypes)),# scale feature matrix \cite{vandenBerg:2006hm}
  # splited_fm=target(transform = map(normalized_fm)),# split feature matrix into train/test parts by spectrumid with respect to percentage
  # pca=target(transform = cross(fm=normalized_fm,color=!!c('diagnosis','spectrumid','patientid'))),# make PCA plots for transformed feature matrices
  # umap=target(transform = cross(fm=normalized_fm,color=!!c('diagnosis','spectrumid','patientid'))),# make umap plots for transformed feature matrices
  rf_cv10=target(train_model(fm=normalized_fm,modeltype='rf'),transform = map(normalized_fm)),# train regression model with CV10
  test_rf=target(test_model(fm=normalized_fm,model=rf),transform = map(rf=rf_cv10)),
  plot_rf=target(plot_test(fm=test_rf),transform = map(test_rf))
  #xgb_cv10=target(train_model(fm=normalized_fm,modeltype='xgb'),transform = cross(normalized_fm),trigger=trigger(condition=train_trigger(fm=normalized_fm)))#,# train regression model with CV10
  #test_xgb=target(test_model(fm=normalized_fm,model=xgb),transform = map(xgb=xgb_cv10)),
  #plot_xgb=target(plot_test(fm=test_xgb),transform = map(test_rf))
  # shap_data=target(transform = map(cv10)),# select small part of scans and calculate data shapley
  # shap_xgb=target(transform = map(shap_data)),# train XGB model to predict data shapley
  # select_shap=target(transform = map(fm=splited_fm,xgb=shap_xgb)), #calculate data shapley for all scans in the FM and keep 75% of highest value
)
