test_that("Testing the number of objects in my output", {
  #testthat::skip_on_cran()
  testthat::skip_if_not_installed("SuperLearner")
  testthat::skip_if_not_installed("survivalSL")
  testthat::skip_if_not_installed("ranger")
  dat = data(DATASET, package = "DRMLSurv")
  dat <- DATASET
  set.seed(123)
  folds <- caret::createFolds(dat$txgroup1L.sd, k = 3)

  test_idx <- folds[[1]]
  train_idx <- setdiff(seq_len(nrow(dat)), test_idx)

  train_data <- dat[train_idx, , drop = FALSE]
  test_data  <- dat[test_idx, , drop = FALSE]
  trainmod = Drmatch(
    data                       = train_data,
    id.var                     = 'patientid',
    eta2.var                   = 'eta2',
    Y1.var                     = 'OS_time.1L',
    Y2.var                     = 'OS_time.2L',
    delta.var                  = 'deathInd.raw',
    OY.var                     = 'OS_time',
    A1.var                     = 'txgroup1L.sd',
    A2.var                     = 'txgroup2L0.sd',
    exact1.vars                = NULL,
    exact2.vars                = NULL,
    names.var1                 = c("ageAt1L.sd","gender.sd","Albumin1st.sd","Lymphocyte1st.sd",
                                   "ECOG1st0.sd","ECOG1st1.sd", 'firstLineStartTime.sd'),
    names.var2                 =  c("ageAt1L.sd","OS_time.1L","Albumin2nd.sd",
                                    "gender.sd","Lymphocyte2nd.sd","ECOG2nd0.sd",
                                    "ECOG2nd1.sd", "txgroup1L.sd", 'firstLineStartTime.sd'),
    usecov                     = FALSE,
    cores                      = 1,
    tau                        = NULL,
    sl.seed                    = 1234,
    A.SL.library              = list(
      "SL.ranger", #Random forest with No screening version
      c("SL.glm", "screen.corRank"), #Logistic regression with screening
      "SL.glm",  #Logistic regression with No screening version
      "SL.glmnet"), #c("SL.glm", "SL.glmnet", "SL.ranger"),
    Y.SL.library              = c("LIB_COXlasso", "LIB_COXall",  "LIB_COXen"
                                  #"LIB_AFTweibull",
                                  #'LIB_AFTllogis',
                                  #"LIB_AFTggamma",
                                  #"LIB_RSF"
    ),
    A.method                   = "method.NNloglik",
    Y.method                   = "ibll",
    plotps                     = FALSE,
    ngrid                      = 5000,
    param.tune                 =  list(
      #LIB_COXlasso = list(lambda=NULL),
      LIB_COXlasso = list(lambda=seq(0.001, 0.25, length.out = 10)),
      #LIB_PHspline = NULL, #list(k=1:4)
      LIB_COXall = NULL,
      #LIB_AFTggamma = NULL,
      #LIB_AFTgamma =  NULL,
      #LIB_COXen = NULL,
      LIB_COXen = list(alpha=seq(.1, .9, length.out = 10),
                       lambda=seq(0.001, 0.1, length.out = 10)) #NULL,
      #LIB_AFTllogis = NULL,
      #LIB_AFTweibull = NULL,
      # LIB_PLANN = list(
      #   inter=1,
      #   size=c(2, 8),
      #   decay=c(0.001, 0.03),
      #   maxit=100,
      #   MaxNWts=10000)

      # LIB_RSF = list(
      #   mtry = unique(round(c(sqrt(length(names.var2)), seq(1, floor(length(names.var2)-1), length.out = 10)))),#round(nX/3), # Number of variables to consider at each split
      #   nodesize = 5,
      #   ntree = c(1000)
      # )
    ),
    param.weights.fix          = NULL,
    param.weights.init         = NULL,
    optim.method               = 'Nelder-Mead',
    maxit                      = 10000,
    usepenalty                 = FALSE,
    runseed                    = 2025,
    useds                      = TRUE,
    Taus                       = NULL,
    adjustdelta1               = FALSE,
    modeltype                  = "ranger",
    usecv                      = TRUE,
    doublepg                   = TRUE,
    model.pg                   = "cox",
    standardize                = FALSE,
    superLearn                 = TRUE,
    pslink                     = "logit",
    pglink                     = "lognormal",
    distance                   = 'mahalanobis',
    method                     = 'nearest',
    K                          = 3,
    replacement                = TRUE,
    cap_months                 = 24
  )

  testthat::expect_false(is.null(trainmod))

  pred <- predict.Drmatch(trainmod, newdata = test_data)

  testthat::expect_false(is.null(pred))
  #testthat::expect_equal(nrow(pred), nrow(test_data))
})

