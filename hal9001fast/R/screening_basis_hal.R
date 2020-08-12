
#' @importFrom qlcMatrix corSparse
#' @export
# Screens a basis list by the pearson correlation value for each basis
# This is an unaggressive screening method designed to simply remove useless basis functions
screen_basis_by_cor = function(basis_list, X, y, pval_cutoff = 0.1, keep = c()){
  cutoff = pval_cutoff
  x_basis = make_design_matrix(X, basis_list)
  cors = qlcMatrix::corSparse(x_basis, Matrix::Matrix(y,ncol=1,sparse=T))
  cors[is.na(cors)] = 0

  t_stat = sapply(cors, function(c){c/sqrt((1-c^2)/(length(y)-2))})
  pvals = sapply(t_stat, function(t){2*pt(-abs(t), length(y)-2)})

  # pvals = (apply(x_basis,2,function(x,y){
  #   if(sd(x)==0){return(1)}
  #   return(cor.test(x,y)$p.value)
  # },y=y))

  keep = union(which(pvals <=cutoff), keep)
  return(list(basis_list[keep],x_basis[,keep] ))
}

# Discretizes/bins the covariates of a matrix by rounding each value to the nearest 1/bins quantile.
#' @export
quantizer = function(X, bins) {
  if (is.null(bins)) {
    return(X)
  }
  X = as.matrix(X)

  convertColumn = function(x) {
    quants = seq(0, 0.97, 1 / bins)
    q = quantile(x, quants)

    nearest <- findInterval(x, q)
    x <- q[nearest]
    return(x)
  }
  quantizer = function(X) {
    as.matrix(apply(X, MARGIN = 2, FUN = convertColumn))
  }
  return(quantizer(X))
}

# Aggressively screens a basis. First performs correlation screening and then uses hal to screen.
# Designed to be used for iteratively screening the basis when the basis list is constructed in a sequential way.
screen_basis = function(basis_list,
                        X,
                        y,
                        index_to_keep = NULL,
                        return_index = F,
                        lower.limits = -Inf,
                        upper.limits = Inf,
                        screen_at_which_lambda = NULL,
                        family = "gaussian", pval_cutoff = 0.05) {
  print(paste0("Current basis size is ", length(basis_list)))

  X = as.matrix(X)

  screen_by_cor_out = screen_basis_by_cor(basis_list,X,y, keep = index_to_keep, pval_cutoff = pval_cutoff)
  x_basis = screen_by_cor_out[[2]]
  basis_list = screen_by_cor_out[[1]]


  if (!is.null(screen_at_which_lambda)) {
    fit = glmnet::cv.glmnet(
      x = x_basis,
      y = y,
      family = family,
      lower.limits = lower.limits,
      upper.limits = upper.limits
    )
    lambda = fit$lambda.min

    #fit =  glmnet::glmnet(x = x_basis, y = y, family = family, lambda =  fit$lambda.min, lower.limits = lower.limits, upper.limits = upper.limits)

  }
  else{
    fit  = glmnet::glmnet(
      x = x_basis,
      y = y,
      family = family,
      lower.limits = lower.limits,
      upper.limits = upper.limits
    )
    lambda = min(fit$lambda)

  }


  if (is.null(screen_at_which_lambda)) {
    betas = fit$beta
    coefs = apply(as.matrix(betas), 2, function(v) {
      as.vector(which(v != 0))
    })
    keep = coefs[[length(coefs)]]

  }
  else if (screen_at_which_lambda == "lambda.min") {
    betas = stats::coef(fit, s = fit$lambda.min)
    keep = which(betas != 0)

  }
  else if (screen_at_which_lambda == "lambda.1se") {
    betas = stats::coef(fit, s = fit$lambda.1se)
    keep = which(betas != 0)
  }


  if (length(keep) == 1) {
    keep = union(keep, c(2, 1))
  }

  len  = length(basis_list)
  keep = setdiff(keep, index_to_keep)


  print(paste0("Amount of basis functions retained: ", length(keep)))
  if (return_index) {
    print(paste0(
      "Current basis size is ",
      length(index_to_keep) + length(keep)
    ))

    return(keep - max(index_to_keep))
  }
  basis_list = basis_list[keep]
  print(paste0(
    "Current basis size is ",
    length(index_to_keep) + length(basis_list)
  ))

  return(basis_list)
}

# Takes two basis lists and generates the interactions obtained by multiplyng basis functios
# whose knot points/cutoffs correspond to the same individual.
merge_basis = function(reduced_basis_list1,
                       reduced_basis_list2,
                       X) {
  if (length(reduced_basis_list1) == 0 | length(reduced_basis_list2)) {
    return(NULL)
  }
  if (is.null(reduced_basis_list1) == 0 |
      is.null(reduced_basis_list2)) {
    return(NULL)
  }
  X = as.matrix(X)
  len1 = length(reduced_basis_list1)
  len2 = length(reduced_basis_list2)
  mergeBasis = function(lst1, lst2) {
    lst = list()
    lst$cols <- c(lst1$cols, lst2$cols)
    order_them = order(lst$cols)
    lst$cols =  lst$cols[order_them]
    lst$cutoffs <- c(lst1$cutoffs, lst2$cutoffs)[order_them]
    lst$orders <- c(lst1$orders, lst2$orders)[order_them]

    if (length(unique(lst$cols)) != length(lst$cols)) {
      return(NULL)
    }
    return(lst)
  }

  map1 = lapply(reduced_basis_list1, function(basis) {
    as.vector(which(X[, basis$cols] == basis$cutoffs))
  })
  map2 = lapply(reduced_basis_list2, function(basis) {
    as.vector(which(X[, basis$cols] == basis$cutoffs))
  })

  #mat = combn(1:length(reduced_basis_list), 2)
  mat = tidyr::crossing(a = 1:len1, b = 1:len2)

  mat = mat[!duplicated(t(apply(mat, 1, sort))), ]

  merged = apply(mat, 1 , function(v) {
    if (any(map1[[v[1]]] %in% map2[[v[2]]])) {
      return(mergeBasis(reduced_basis_list1[[v[1]]], reduced_basis_list2[[v[2]]]))
    }
    else{
      return(NULL)
    }
  })
  return(merged[!sapply(merged, is.null)])

}

# From a basis list (containing main-term one way basis functions), generate higher order interactions corresponding with observed cutoffs/knot points.
get_higher_basis = function(reduced_basis_list,
                            max_dim,
                            X,
                            y,
                            screen_each_level = F,
                            max_num_two_way = Inf, pval_cutoff = 0.05) {
  if (length(reduced_basis_list) == 0) {
    return(NULL)
  }
  if (max_dim == 1) {
    return(reduced_basis_list)
  }

  basis_lists = get_higher_basis_up_to_three(reduced_basis_list, max_dim, X, y, screen_each_level, max_num_two_way, pval_cutoff = pval_cutoff)



  final = c(reduced_basis_list, basis_lists[[1]])

  merged = basis_lists[[2]]

  if (max_dim <= 3) {
    result = (c(final, basis_lists[[2]]))
    if (F) {
      return(final)
    }
    else{

      return(result)
    }
  }
  for (i in 1:(max_dim - 3)) {
    if (!is.null(merged) & screen_each_level) {
      merged = screen_basis(c(final, merged), X, y, index_to_keep = 1:length(final))
    }
    merged = merge_basis(reduced_basis_list, merged, X)


    if (!is.null(merged) & length(merged) != 0) {
      final = c(final, merged)
    }
    else{
      break
    }

  }
  return(final)
}


# Hard-coded but fast way to generate the second and third degree interactions from a main term basis list
get_higher_basis_up_to_three = function(reduced_basis_list,
                                        max_dim,
                                        X,
                                        y,
                                        screen_each_level, max_basis = 20000, pval_cutoff=0.05) {
  if (max_dim == 1 | ncol(X) == 1) {
    return(list(list(), list()))
  }

  getBasis = function(vec) {
    lst = list()
    lapply(vec, function(l) {
      basis =  reduced_basis_list[[l]]
      lst$cols <<- c(lst$cols, basis$cols)
      lst$cols <<-  lst$cols
      lst$cutoffs <<- c(lst$cutoffs, basis$cutoffs)
      lst$orders <<- c(lst$orders, basis$orders)
    })
    if (length(unique(lst$cols)) != length(lst$cols)) {
      return(NULL)
    }
    ordering = order(lst$cols)
    lst$cols = lst$cols[ordering]
    lst$cutoffs = lst$cutoffs[ordering]
    lst$orders = lst$orders[ordering]
    return(lst)
  }

  map = lapply(reduced_basis_list, function(basis) {
    as.vector(which(X[, basis$cols] == basis$cutoffs))
  })

  mat = combn(1:length(reduced_basis_list), 2)

  get_valid_two_way_ind = apply(mat, 2, function(v) {
    any(map[[v[1]]] %in% map[[v[2]]])
  })
  if (!any(get_valid_two_way_ind)) {
    return(list(list(), list()))
  }
  two_way_combos = mat[, get_valid_two_way_ind, drop = F]

  if (ncol(two_way_combos) == 0) {
    return(list(list(), list()))
  }
  way = apply(two_way_combos, 2, getBasis)
  throw = !sapply(way, is.null)
  way = way[throw]
  two_way_combos = two_way_combos[, throw, drop = F]
  t = proc.time()
  print(paste0("twoways (pre cor screening): ", length(way)))

  way = screen_basis_by_cor(way,X,y, pval_cutoff =pval_cutoff )[[1]]
  print(paste0("twoways (after cor screening): ", length(way)))

  if(length(way) > max_basis){
    #screen_each_level=T
    x_bas = make_design_matrix(X, way)
    cor_vals = apply(x_bas, 2, cor, y)
    way = way[which(rank(-abs(cor_vals))<= max_basis)]
  }
  if (screen_each_level) {


    keep = screen_basis(basis_list=c(reduced_basis_list, way),
                        X=X,
                        y=y,
                        index_to_keep=1:length(reduced_basis_list),
                        return_index=T, pval_cutoff = pval_cutoff)

    two_way_combos = two_way_combos[, keep, drop = F]
    way = way[keep]

  }

  if (max_dim == 2 | ncol(X) == 2) {
    return(list(way, list()))
  }

  res = lapply(unique(two_way_combos[1, ]), function(init) {
    matches = c(two_way_combos[2, two_way_combos[1, ] == init])
    triples = unlist(lapply(matches, function(n) {
      inter = c(intersect(two_way_combos[2, two_way_combos[1, ] == n], matches))
      if (length(inter) == 0) {
        return()
      }
      result = lapply(inter, function(s) {
        combo = c(init, n, s)
        if (length(intersect(map[[init]],
                             intersect(map[[n]],
                                       map[[s]])) != 0)) {
          return(combo)
        }
        return(NULL)
      })
      return(result[!sapply(result, is.null)])
    }), recursive = F)
    if (is.null(triples) | length(triples) == 0) {
      return()
    }
    triples =  triples[!sapply(triples, is.null)]
    if (length(triples) == 0 | is.null(triples)) {
      return()
    }

    as.matrix(do.call(cbind, triples))
  })

  res = res[!sapply(res, is.null)]


  if (length(res) == 0) {
    print(paste0("three: ",0))
    return(list(way, list()))
  }
  keeper = unlist(sapply(res, function(v) {
    nrow(v) != 0
  }))
  res = res[keeper]
  res = as.matrix(do.call(cbind, res))
  print(paste0("three: ", ncol(res)))
  new_basis = apply(res, 2, getBasis)


  new_basis = new_basis[!sapply(new_basis, is.null)]
  if (is.null(new_basis)) {
    new_basis = list()
  }
  up_to_three = list(way, new_basis)
  return(up_to_three)
}
