# ENM evaluation function
## ML Gaynor

# Load library
library(kuenm)

# proc function
# Function from https://jamiemkass.github.io/ENMeval/articles/ENMeval-2.0.0-vignette.html#running-enmeval
proc <- function(vars) {
  proc <- kuenm::kuenm_proc(vars$occs.val.pred, c(vars$bg.train.pred, vars$bg.val.pred))
  out <- data.frame(proc_auc_ratio = proc$pROC_summary[1], 
                    proc_pval = proc$pROC_summary[2], row.names = NULL)
  return(out)
}


