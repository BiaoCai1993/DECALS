compare_cross_and_within <- function(pseudo_bulk_sc, weighting = F, filtering = F, method = 'pearson',
                                     gaussian_matching = F){
  K <- length(pseudo_bulk_sc)
  p <- nrow(pseudo_bulk_sc[[1]])
  within_sum <- between_sum <- matrix(0, p, p)
  cor_p_val <- cov_all_cts <- array(NA, dim = c(p, p, K, K))
  ct_props_matched <- ct_props[match(names(pseudo_bulk_sc), names(ct_props))]
  print(ct_props_matched)
  if(gaussian_matching){
    for(k1 in 1:K){
      for(i in 1:p){
        tmp <- pseudo_bulk_sc[[k1]][i, ]
        tmp_qt <- rank(tmp) / length(tmp)
        tmp_qm <- qnorm(tmp_qt, mean(tmp), sd(tmp))
        tmp_qm[is.infinite(tmp_qm)] <- mean(tmp) + 3* sd(tmp)
        pseudo_bulk_sc[[k1]][i, ] <- tmp_qm
        
      }
    }
  }
  # gene 1
  for(i in 1:p){
    # gene 2
    for(j in 1:p){
      # cell type 1
      for(k1 in 1:K){
        # within-cell-type
        exp1 <- pseudo_bulk_sc[[k1]][i, ]
        exp2 <- pseudo_bulk_sc[[k1]][j, ]
        if(sd(exp1) > 0 & sd(exp2) > 0){
          if(weighting){
            cov_12 <- cov(exp1, exp2) * ct_props_matched[k1]^2
          }else{
            cov_12 <- cov(exp1, exp2)
          }
          pval <- cor.test(exp1, exp2, method = method)$p.value
          if(filtering){
            cov_12 <- ifelse(pval < 0.05, cov_12, 0)
          }
          cov_all_cts[i, j, k1, k1] <- cov_12
          within_sum[i, j] <- within_sum[i, j] + cov_12
          cor_p_val[i, j, k1, k1] <- pval
        }
        
        for(k2 in 1:K){
          # between-cell-type covariance
          if(k1 != k2){
            exp1 <- pseudo_bulk_sc[[k1]][i, ]
            exp2 <- pseudo_bulk_sc[[k2]][j, ]
            if(sd(exp1) > 0 & sd(exp2) > 0){
              if(weighting){
                cov_12 <- cov(exp1, exp2) * ct_props_matched[k1]  * ct_props_matched[k2]
              }else{
                cov_12 <- cov(exp1, exp2)
              }
              pval <- cor.test(exp1, exp2, method = method)$p.value
              if(filtering){
                cov_12 <- ifelse(pval < 0.05, cov_12, 0)
              }
              between_sum[i, j] <- between_sum[i, j] + cov_12
              cov_all_cts[i, j, k1, k2] <- cov_12
              cor_p_val[i, j, k1, k2] <- pval
            }
          }
        }
      }
    }
  }
  return(list(within = within_sum, between = between_sum, cov_rec = cov_all_cts,
              cor_p_val = cor_p_val))
}

extract_within_and_between <- function(within_cov){
  within_ct_ind <- array(F, dim = dim(within_cov))
  between_ct_ind <- array(T, dim = dim(within_cov))
  for(i in 1:dim(within_cov)[1]){
    for(j in 1:dim(within_cov)[2]){
      for(k1 in 1:dim(within_cov)[3]){
        for(k2 in 1:dim(within_cov)[4]){
          if(k1 == k2){
            within_ct_ind[i, j, k1, k2] <- T
            between_ct_ind[i, j, k1, k2] <- F
          }
        }
      }
    }
  }
  return(list(within = within_ct_ind, between = between_ct_ind))
}

eval_within_and_between <- function(res, inds, method = 'BH', p_cutoff = 0.05, thresholding = T, geneset = '',
                                    min_threshold = 1e-6, title = ''){
  if(method == 'BH'){
    within_p_adjusted <- p.adjust(res$cor_p_val, method = 'BH')[inds$within] 
    between_p_adjusted <- p.adjust(res$cor_p_val, method = 'BH')[inds$between]
  }else{
    within_p_adjusted <- res$cor_p_val[inds$within]
    between_p_adjusted <- res$cor_p_val[inds$between]
  }
  ngenes <- dim(inds$within)[1]
  p_cutoff <- ifelse(method == 'BH', p_cutoff, p_cutoff / (ngenes^2))
  if(!thresholding) p_cutoff = 1
  within <- res$cov_rec
  within[!inds$within] <- 0
  within[inds$within][within_p_adjusted > p_cutoff] <- 0
  within_sum <- apply(within, c(1,2), function(x) sum(x,na.rm=T))
  
  between <- res$cov_rec
  between[!inds$between] <- 0
  between[inds$between][between_p_adjusted > p_cutoff] <- 0
  between_sum <- apply(between, c(1,2), function(x) sum(x,na.rm=T))
  
  # compare for all gene pairs (j,j')
  cov_sum_df <- data.frame(within = within_sum[upper.tri(within_sum,diag=TRUE)] %>% abs, 
                           between = between_sum[upper.tri(between_sum,diag=TRUE)] %>% abs)
  g1 <- ggplot(cov_sum_df) +
    geom_jitter(aes(x = log10(between+1), y = log10(within+1)), alpha = 0.2, width = 0.1, height = 0.1) +
    geom_abline(slope = 1, intercept = 0, color = 'red') +
    theme_classic(base_size = 20) 
  
  print(sprintf('total=%i, between=0: %i, within=0:%i, both: %i',
                nrow(cov_sum_df),
                sum((cov_sum_df$within != 0 & cov_sum_df$between == 0)),
                sum((cov_sum_df$within == 0 & cov_sum_df$between != 0)),
                sum((cov_sum_df$within == 0 & cov_sum_df$between == 0))))
  
  # focus on gene pairs that have both within and between !=0
  
  g2 <- #ggplot(cov_sum_df[(abs(cov_sum_df$within) > min_threshold & abs(cov_sum_df$between) > min_threshold), ]) +
    ggplot(cov_sum_df[(cov_sum_df$within !=0 & cov_sum_df$between != 0), ]) +
    geom_histogram(aes(x = log10(within/between)), alpha = 0.2, bins = 30) +
    geom_vline(xintercept = 0, color = 'red') +
    theme_classic(base_size = 20) 
  
  
  sprintf('method: %s, prop of within equal to between: %.2f; & greater than between: %.2f (%.2f)', 
          method, mean((within_sum %>% abs) == (between_sum %>% abs), na.rm=T),
          mean((within_sum %>% abs) > (between_sum %>% abs), na.rm=T),
          sum((within_sum %>% abs) > (between_sum %>% abs), na.rm=T) / (sum((within_sum %>% abs) != (between_sum %>% abs), na.rm=T))) %>% print
  # return(cov_sum_df)
  
  # shift within and between by 1 to avoid log(0).
  g3 <- ggplot(cov_sum_df) +
    geom_histogram(aes(x = log10(within+1) - log10(between+1)), alpha = 0.2, bins = 30) +
    geom_vline(xintercept = 0, color = 'red') +
    labs(x = 'log10(in+1)-log10(bt+1)') +
    theme_classic(base_size = 20)   
  g <- grid.arrange(g1, g2, g3, nrow = 1, 
                    top = textGrob(title, gp=gpar(fontsize=25,font=8)))
  
}

options(repr.plot.width = 14, repr.plot.height = 4.5)
plot_within_versus_between <- function(gene_set, filtering, gene_set_name){
  if(filtering){
    add_inds <- top_inds[['10000']]
  }else{
    add_inds <- rep(T, length(features))
  }
  scaled_data_list <- lapply(major_cts, function(ct){
    data_list[[ct]]$scaled_data[features %in% gene_set & add_inds, ]
  })
  names(scaled_data_list) <- major_cts
  
  print(sprintf('Total genes: %i, after intersection with genes highly expressed in all four cell types: %i', 
                length(gene_set), nrow(scaled_data_list[[1]])))
  
  res <- compare_cross_and_within(scaled_data_list, weighting = T, filtering = F, method = 'pearson',
                                  gaussian_matching = F)
  inds <- extract_within_and_between(res$cov_rec)
  
  eval_within_and_between(res, inds, thresholding=F, title = sprintf('%s genes, no thresholding', gene_set_name))
  eval_within_and_between(res, inds, thresholding=T, title = sprintf('%s genes, BH p-values', gene_set_name))
}

#plot_within_versus_between(Biao_genes, T, 'DECALS')
