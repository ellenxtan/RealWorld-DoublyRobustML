LstToDF <- function(lst, df_name) {
  out_lst <- purrr::map(lst, df_name)
  # out_df <- purrr::map_dfr(out_lst, rbind)
  out_df <- c()
  for (r in 1:length(out_lst)) {
    out_df <- rbind(out_df, out_lst[[r]])
  }
  return(out_df)
}
