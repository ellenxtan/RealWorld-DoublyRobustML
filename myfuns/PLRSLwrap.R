# # DML estimator, regress the residual of E(Y|X) on the residual of E(A|X) with cross-fitting. 
# # In this way, we can choose different models for E(Y|X) and E(A|X).
# # PLR-SL: not do variable selection. Just use all X to fit E(Y|X) and E(A|X) with SL then calculate ATE
# PLRSLwrap <- function(mydat, useCovar) {
# 
#   covars_all <- grep(paste0("^",useCovar), names(mydat), value=TRUE)
#   
#   if (length(unique(mydat$Y)) != 2) {  # continuous Y
#     sim_fam <- gaussian()
#   } else {  # binary Y
#     sim_fam <- binomial()
#   }
#   
# 
#   
#   # Load bonus data
#   df_bonus = fetch_bonus(return_type="data.table")
#   
#   # Specify the data and variables for the causal model
#   dml_data_bonus = DoubleMLData$new(df_bonus,
#                                     y_col = "inuidur1",
#                                     d_cols = "tg",
#                                     x_cols = c("female", "black", "othrace", "dep1", "dep2",
#                                                "q2", "q3", "q4", "q5", "q6", "agelt35", "agegt54",
#                                                "durable", "lusd", "husd"))
#   print(dml_data_bonus)
#     
#   
#   # surpress messages from mlr3 package during fitting
#   lgr::get_logger("mlr3")$set_threshold("warn")
#   lgr::get_logger("bbotk")$set_threshold("warn")
#   
#   # DML functions (paper replication code)
#   # file:///Users/xtan/Downloads/HESR-54-1273-s002.pdf
#   # https://onlinelibrary.wiley.com/doi/abs/10.1111/ectj.12097
#   
#   # https://github.com/mlr-org/mlr3extralearners/blob/main/R/learner_dbarts_regr_bart.R
#   # https://github.com/mlr-org/mlr3extralearners/blob/main/R/learner_mgcv_regr_gam.R
#   
#   # ps
#   learner = lrns("classif.bart", predict_type="prob")#, k=0.5
#   learner = lrns(c("classif.glmnet","classif.gam","classif.bart"), predict_type="prob")#, k=0.5
#   ml_m_bonus = learner$clone()
#   
#   # or
#   learner = lrns(c("regr.glm","regr.gam","regr.bart"), k=2)
#   ml_g_bonus = learner$clone()
#   
#   set.seed(3141)
#   obj_dml_plr_bonus = DoubleMLPLR$new(dml_data_bonus, ml_g=ml_g_bonus, ml_m=ml_m_bonus)
#   obj_dml_plr_bonus$fit()
#   print(obj_dml_plr_bonus)
#   
#   obj_dml_plr_bonus$coef
#   obj_dml_plr_bonus$summary()
#   obj_dml_plr_bonus$confint()
# }
