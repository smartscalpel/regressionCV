plan <- drake_plan(
  peaks=target(),# read peak files and organise them as dynamic target
  fm=target(prepare_feature_matrix(peaks),transform = map(peaks)),# convert peaks to feature matrix
  normalized_fm=target(transform = cross(fm,!!c('None','Autoscaling','Pareto'))),# scale feature matrix \cite{vandenBerg:2006hm}
  splited_fm=target(transform = map(normalized_fm)),# split feature matrix into train/test parts by spectrumid with respect to percentage
  pca=target(transform = cross(fm=normalized_fm,color=!!c('diagnosis','spectrumid','patientid'))),# make PCA plots for transformed feature matrices
  umap=target(transform = cross(fm=normalized_fm,color=!!c('diagnosis','spectrumid','patientid'))),# make umap plots for transformed feature matrices
)
