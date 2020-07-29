
#' @export
formula_hal9001 <- function(form, data, smoothness_orders = NULL, include_zero_order = F, bins = NULL){
  if(is.null(smoothness_orders)| !is.numeric(smoothness_orders)){
    smoothness_orders = round(rep(0,ncol(X)))
  }else{
    #recycle vector if needed.
    smoothness_orders = suppressWarnings(round(smoothness_orders) + rep(0,ncol(X)))
  }
  order_map = smoothness_orders
  form = stringr::str_replace_all(form, " ", "")

  reg = "([^\\s])~([idh]\\(([^\\s()+]+)\\)|\\.(\\^[0-9])?)(?:\\+[ihd]\\(([^\\s()]+)\\))*(\\+\\.(\\^[0-9])?)?$"
  assertthat::assert_that(stringr::str_detect(form,reg), msg = "Incorrect format for formula.")
  outcome = stringr::str_match(form, "([^\\s]+)~")
  outcome = outcome[,2]
  X = data[, -which(colnames(data)==outcome), drop = F]
  X_orig = X
  X = quantizer(X, bins)
  y = data[,outcome]
  names = colnames(X)

  interactions = stringr::str_match_all(form, "[ihd]\\(([^\\s()]+)\\)")[[1]]
  interactions = interactions[,2]

  typeMonotone = stringr::str_match_all(form, "([[a-z]])\\([^\\s()]+\\)")[[1]]
  typeMonotone = typeMonotone[,2]

  if(stringr::str_detect(form, "\\.")){
    degree_rest = as.numeric(stringr::str_match(form, "\\.\\^([0-9]+)")[,2])
    if(is.na(degree_rest)) degree_rest = 1
  }
  else{
    degree_rest = NULL
  }

  getCombos = function(deg){
    x = combn(1:length(names), deg)
    allCombos = lapply(seq_len(ncol(x)), function(i) x[,i])

  }
  if(is.null(degree_rest)){
    allCombos = list()
  }
  else{
    allCombos = unlist(lapply(1:degree_rest,getCombos), recursive  =F)

  }
  getCov = function(term){
    cols = unlist(stringr::str_extract_all(term, "[^,]"))

    ind = unlist(lapply(as.vector(cols), function(x){(which(names ==x))}))
    ind = ind[unlist(lapply(ind, function(v){length(v)!=0}))]
    ind = sort(unique(ind))

    return( ind)
  }
  interactions_index = (lapply(interactions, getCov))
  interactions_index = interactions_index[unlist(lapply(interactions_index, function(v){length(v)!=0}))]

  not_dupes_index = which(!duplicated(interactions_index))


  interactions_index = interactions_index[not_dupes_index]
  typeMonotone = typeMonotone[not_dupes_index]


  outcome_index = getCov(outcome)



  allCombosLeft = setdiff(allCombos,interactions_index)

  lower.limits = c()
  upper.limits = c()
  basis_list = list()

  for(i in 1:length(interactions_index)){
    if(length(interactions_index)==0) break
    new_basis = basis_list_cols(interactions_index[[i]], X, order_map, include_zero_order)
    if(typeMonotone[i]== "i"){
      lower.limits = c(lower.limits, rep(0, length(new_basis)))
      upper.limits = c(upper.limits, rep(Inf, length(new_basis)))
    }
    else if (typeMonotone[i]== "d"){
      lower.limits = c(lower.limits, rep(-Inf, length(new_basis)))
      upper.limits = c(upper.limits, rep(0, length(new_basis)))
    }
    else{
      lower.limits = c(lower.limits, rep(-Inf, length(new_basis)))
      upper.limits = c(upper.limits, rep(Inf, length(new_basis)))
    }
    basis_list = c(basis_list, new_basis)
  }

  basis_listrest = unlist(lapply(allCombosLeft, basis_list_cols, X, order_map, include_zero_order ), recursive = F)
  len = length(basis_listrest)
  upper.limits = c(upper.limits, rep(Inf, len))
  lower.limits = c(lower.limits, rep(-Inf, len))
  basis_list = c(basis_list, basis_listrest)

  form_obj = list()
  form_obj$formula = form
  form_obj$basis_list = basis_list
  form_obj$upper.limits = upper.limits
  form_obj$lower.limits = lower.limits
  form_obj$X = as.matrix(X_orig)
  form_obj$y = as.vector(y)
  class(form_obj) <- "formula_hal9001"
  return(form_obj)
}




