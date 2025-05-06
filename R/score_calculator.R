
#devtools::use_data(diamonds, overwrite = TRUE)
internal_env <- new.env()
df_name <- load("data/all_para_tables.RData", internal_env)
#df_name <- load("all_para_tables.RData", internal_env)
#df_name_2 <- load("par_bmi_tables.RData", internal_env)
#get(x, envir=internal_env)
#rm(internal_env)

# initialize
sex_cats <- c("boys", "girls")
var_names <- c("bmi", "dbp", "glu", "hdl", "height", "homa", "insu", "MetS_shifted", "sbp", "trg", "waist")
par_cats <- c("mu", "sigma", "nu", "tau")

# approxfun for bivariate functions
approxfun2 <- function(x, y, Z, method = "linear") {
  force(x)
  force(y)
  force(Z)
  force(method)

  function(age, height) {
    #handle NAs
    none_na <- !(is.na(age) | is.na(height))
    rs <- rep(NA,length(none_na))
    if (sum(none_na)>0)
      rs[none_na] <- pracma::interp2(x, y, Z, age[none_na], height[none_na], method)
    return(rs)
  }
}

# vector for the distributions per variable
var_distribution <- list()

# fitsplines for each parameter, variable and sex
approx_param_functions <- list()
for (sex in sex_cats) {
  for (vname in var_names) {
    # get table
    curr_tab <- get(paste("par", vname, sex, sep="_"), envir= internal_env)

    # remember distribution type
    if (is.null(var_distribution[[sex]])) var_distribution[[sex]] <- list()
    var_distribution[[sex]][[vname]] <- as.character(unique(curr_tab[["dist"]])[1])

    age_values    <- curr_tab[["age"]]
    height_values <- curr_tab[["height"]]
    curr_cats <- intersect(par_cats, colnames(curr_tab))
    for (param in curr_cats){
      parameters <- curr_tab[[param]]
      #  force(param)
      if (!is.null(height_values)) {
        parameter_matrix <- tapply(curr_tab[[param]], list(curr_tab$height, curr_tab$age), FUN = identity)
        approx_param_functions[[paste(vname,sex,param, sep="_")]] <- approxfun2(unique(age_values), unique(height_values), parameter_matrix)
        rm(parameter_matrix)

        ## vvv hard to vectorize...
        #parameters <- curr_tab[[param]]
        #unique_ages <- unique(age_values)
        #for (age in unique_ages) {
        #  age_indices <- age_values == age
        #  approx_param_functions[[paste(vname,sex,param,age, sep="_")]] <- approxfun(height_values[age_indices], parameters[age_indices], ties=mean)
        #}
        #approx_param_functions[[paste(vname,sex,param, sep="_")]] <- function(age, height) {
        #  lowerAge <- floor(age*10)/10
        #  upperAge <- ceiling(age*10)/10
        #
        #  lowerAgeFun <- approx_param_functions[[paste(vname,sex,param,lowerAge, sep="_")]]
        #  upperAgeFun <- approx_param_functions[[paste(vname,sex,param,upperAge, sep="_")]]
        #
        #  if (is.null(lowerAgeFun) || is.null(upperAgeFun)) {
        #    warning("age out of range")
        #  }
        #  i_weight <- (age-lowerAge)/(upperAge-lowerAger)
        #  return(i_weight*upperAgeFun(height) +  (1-i_weight)*lowerAgeFun(height))
        #}

      } else {
        approx_param_functions[[paste(vname,sex,param, sep="_")]] <- approxfun(age_values, parameters, ties=mean)
      }
    }
  }
}

mkfun_sex_age <- function(vname, param) {
  function(sex, age) {
    ifelse(is.na(sex) | is.na(age),
           NA,
          ifelse(sex=="m",
                 approx_param_functions[[paste(vname,"boys", param, sep="_")]](age),
                 approx_param_functions[[paste(vname,"girls",param, sep="_")]](age))
           )

  }
}

mkfun_sex_age_height <- function(vname, param) {
  function(sex, age, height) {
    ifelse(is.na(sex) | is.na(age) | is.na(height),
           NA,
          ifelse(sex=="m",
                 approx_param_functions[[paste(vname,"boys", param, sep="_")]](age, height),
                 approx_param_functions[[paste(vname,"girls",param, sep="_")]](age, height))
    )
  }
}

for (vname in var_names) {
  all_cols <- lapply(sex_cats, function(x) colnames(get(paste("par", vname, x, sep="_"), envir= internal_env)))
  cols_for_all_sexes <- Reduce(intersect, all_cols)
  curr_cats <- intersect(par_cats, cols_for_all_sexes)
  for (param in curr_cats) {

    fname <- paste(vname, param, sep = "_")

    # Use a local scope to capture vname and param
    local({
      vname_ <- vname
      param_ <- param

      if (!"height" %in% cols_for_all_sexes) {
        approx_param_functions[[fname]] <<- mkfun_sex_age(vname_, param_)
      } else {
        approx_param_functions[[fname]] <<- mkfun_sex_age_height(vname_, param_)
      }
    })
  }
}

#approx_param_functions$dbp_boys_sigma(5, 120)
#approx_param_functions$dbp_mu(c("f","m","f","m"),c(5,5,5,5), c(120,120,121,121))
#approx_param_functions$bmi_sigma(c("f","m"),c(5,5))

#' Calculate percentiles and z-scores
#'
#' The function calculates the percentiles and z-scores for given value of
#' several clinical parameters ('waist', 'sbp', 'dbp', 'trg', 'hdl' and 'homa',
#' 'bmi', 'glu', 'height', 'insu')
#' using the respective sex-specifc (and for some parameters height-specific)
#' reference table.
#
# Arguments:
## data_input: data set with study data including the clinical parameters
## sex: sex of the subjects in the data set ('m' or 'f')
## p: clinical parameter that shall be investigated (use specific abbreviation)
## tablepath: the path of the file 'all_para_tables.RData'
#'
#' @param variable Either 'waist', 'sbp', 'dbp', 'trg', 'hdl' and 'homa',
#' 'bmi', 'glu', 'height' or 'insu'.
#' @param sex A vector of characters 'f' and 'm'.
#' @param age A numeric vector.
#' @param height A numeric vector.
#' @param values A numeric vector with values for the parameter given in `variable`.
#' @param return_values A character vector with the possible values 'percentile' and 'z.score'.
#' @return A A list with the values specified in `return_values`
#' @examples
#' get_scores(variable="dbp", sex=c("f","m"), age=c(5,5), height=c(120,110), values=c(70,60))
get_scores <- function(variable="waist", sex=c("f","m"), age=6:5, height=NULL, values=c(20,21), return_values=c("percentile","z.score"))  {
  if (length(variable)>1) warning("variable must be of length 1. Only the first value of variable is used.")

  sex_map <- c(f = "girls", m = "boys")
  dist <- var_distribution[[sex_map[sex][1] ]][[variable]]

  all_cols <- lapply(sex_cats, function(x) colnames(get(paste("par", variable, x, sep="_"), envir= internal_env)))
  cols_for_all_sexes <- Reduce(intersect, all_cols)
  curr_cats <- intersect(par_cats, cols_for_all_sexes)

  # assign parameters mu, sigma, nu, tau if used by current model
  if ("height" %in% cols_for_all_sexes)
    for (param in curr_cats) assign(param, approx_param_functions[[paste(variable,param,   sep="_")]](sex, age, height) )
  else
    for (param in curr_cats) assign(param, approx_param_functions[[paste(variable,param,   sep="_")]](sex, age) )
  #mu    <- approx_param_functions[[paste(variable,"mu",   sep="_")]](sex, age)
  #sigma <- approx_param_functions[[paste(variable,"sigma",sep="_")]](sex, age)
  #nu    <- approx_param_functions[[paste(variable,"nu",   sep="_")]](sex, age)
  #tau   <- approx_param_functions[[paste(variable,"tau",  sep="_")]](sex, age)

  if (dist == "BCCG") {
    #handle NAs
    none_na <- !(is.na(values) | is.na(mu) | is.na(sigma) | is.na(nu))

    percentiles <- rep(NA,length(none_na))
    percentiles[none_na] <- gamlss.dist::pBCCG(q = values[none_na], mu[none_na], sigma[none_na], nu[none_na])
  } else if (dist == "BCT") {

    #handle NAs
    none_na <- !(is.na(values) | is.na(mu) | is.na(sigma) | is.na(nu) | is.na(tau))

    percentiles <- rep(NA,length(none_na))
    percentiles[none_na] <- gamlss.dist::pBCT(q = values[none_na], mu[none_na], sigma[none_na], nu[none_na], tau[none_na])

  } else if (dist == "BCPE") {

    #handle NAs
    none_na <- !(is.na(values) | is.na(mu) | is.na(sigma) | is.na(nu) | is.na(tau))

    percentiles <- rep(NA,length(none_na))
    percentiles[none_na] <- gamlss.dist::pBCPE(q = values[none_na], mu[none_na], sigma[none_na], nu[none_na], tau[none_na])

  } else if (dist == "LO") {

    #handle NAs
    none_na <- !(is.na(values) | is.na(mu) | is.na(sigma))

    percentiles <- rep(NA,length(none_na))
    percentiles[none_na] <- gamlss.dist::pLO(q = values[none_na], mu[none_na], sigma[none_na])

  }

  # z.score requested?
  if ("z.score" %in% return_values) {
    zscores <- rep(NA,length(none_na))
    zscores[none_na] <- gamlss.dist::qNO(percentiles[none_na])
  }

  # build up return values as list
  rs <- list()
  if ("percentile" %in% return_values) rs$percentile <- percentiles
  if ("z.score" %in% return_values) rs$z.score <- zscores
  return(rs)
}

#' Calculate the MetS-score for a given data set
#'
#' The function takes the z-scores of the clinical parameters 'waist', 'sbp',
#' dbp', 'trg', 'hdl' and 'homa' to calculate the individual MetS-score
#' @param df data set with z-scores of of the clinical parameters, named:
#' hdl_z.score height_z.score homa_z.score trg_z.score waist_z.score sbp_z.score dbp_z.score
#' @return A numeric vector with the values of the MetS scores.
#' @examples
#' z.score_data <- data.frame(
#' hdl_z.score = c(-1.22,0.10,0.74,-0.39),
#' height_z.score = c(0.88,0.47,-0.16,1.09),
#' homa_z.score = c(-0.15,1.38,0.79,0.60),
#' trg_z.score = c(0.28,0.16,0.16,0.97),
#' waist_z.score = c(1.53,2.12,0.39,1.57),
#' sbp_z.score = c(0.85,-0.98,-0.53,-0.11),
#' dbp_z.score = c(0.64,-1.85,-0.61,0.41)
#' )
#' MetSScore(z.score_data)
MetSScore <- function(df) {
  if( is.null(df$waist_z.score) || is.null(df$homa_z.score) ||
      is.null(df$sbp_z.score) || is.null(df$dbp_z.score) || is.null(df$trg_z.score) || is.null(df$hdl_z.score) )
    stop("MatS score calculation requires the following columns not to be NULL: waist_z.score, homa_z.score, sbp_z.score, dbp_z.score, trg_z.score, hdl_z.score")
  MetS <- df$waist_z.score + df$homa_z.score +
    0.5*(df$sbp_z.score + df$dbp_z.score + df$trg_z.score - df$hdl_z.score)

  ## for distribution purposes: shifted value to calculate percentile and z-score
  #return_df$MetS_shifted <- return_df$MetS + 100
  return(MetS)
}


#' Derive the monitoring/action level status for component 'adiposity'
#'
#' The function decides for the 'adiposity' component whether a certain limit of
#' the percentile of 'waist' has been exceeded to classify a child (none,
#' monitoring, action)
#' If filter is NULL than all those action levels are determined, for which the necessary z.scores are ub df. If a vector with a filter string is defined, but the corresponding required z.scores are unavailable, NAs will be returned.
#' @param df data set with percentile values for adiposity component
#' @param lvl_name name of the classification level ('monit' or 'action')
#' @param perc_level the appropriate cut-off percentile for the the corresponding
#' @param append logical, if TRUE then the input data df are also returned in the rrsulting data.frame
#' @param filter NULL or a vector with any of the values "adiposity", "blood_pressure", "blood_lipids", "blood_glu_insu", "overall"
#' classification (0.9 or 0.95)
#' @return A ...
#' @examples
#' z.perc_data <- data.frame(
#' hdl_percentile = c(0.1,0.5,0.7,0.3),
#' homa_percentile = c(0.4,0.9,0.8,0.7),
#' trg_percentile = c(0.6,0.5,0.5,0.8),
#' waist_percentile = c(0.9,0.99,0.6,0.95),
#' sbp_percentile = c(0.8,0.1,0.3,0.5),
#' dbp_percentile = c(0.7,0.01,0.2,0.6)
#' )
#'
#' # determine available action levels
#' action_levels(z.perc_data)
#'
#' # get only overall adiposity and blood_pressure action level
#' action_levels(z.perc_data, filter=c("adiposity","blood_pressure"))
#'
#' # get only overall action level
#' action_levels(z.perc_data, filter="overall")
action_levels <- function(df, lvl_name=c("none","monit","action"), perc_level=c(0.9, 0.95), append=FALSE, filter=NULL) {
  n <- nrow(df)

  rs <- list()

  # helper function to determine action levels from a single percentile
  perc_to_actlev <- function(perc) {
    if (length(perc) == 0) # catches NULL as well as 1-NULL
      rep(NA, n)
    else
      cut(perc, c(-Inf,perc_level,Inf), lvl_name, ordered_result = TRUE)
  }

  if ("adiposity" %in% filter || "overall" %in% filter || (is.null(filter) && !is.null(df$waist_percentile)))
    rs$adiposity.action <- perc_to_actlev(df$waist_percentile)

  if ("blood_pressure" %in% filter || "overall" %in% filter || (is.null(filter) && !(is.null(df$dbp_percentile) && is.null(df$sbp_percentile) ))) {
    dbp.action <- perc_to_actlev(df$dbp_percentile)
    sbp.action <- perc_to_actlev(df$sbp_percentile)
    rs$blood_pressure.action <- pmax(dbp.action,sbp.action, na.rm=T)
  }

  if ("blood_lipids" %in% filter || "overall" %in% filter || (is.null(filter) && !(is.null(df$trg_percentile) && is.null(df$hdl_percentile) ))) {
    trg.action <- perc_to_actlev(df$trg_percentile)
    hdl.action <- perc_to_actlev(1-df$hdl_percentile)
    rs$blood_lipids.action <- pmax(trg.action,hdl.action, na.rm=T)
  }

  if ("blood_glu_insu" %in% filter || "overall" %in% filter || (is.null(filter) && !(is.null(df$homa_percentile) && is.null(df$glu_percentile) ))) {
    homa.action <- perc_to_actlev(df$homa_percentile)
    glu.action <- perc_to_actlev(df$glu_percentile)
    rs$blood_glu_insu.action <- pmax(homa.action,glu.action, na.rm=T)
  }

  if (is.null(filter) || "overall" %in% filter) {
    rs$overall.action <-
      apply(cbind(rs$adiposity.action,rs$blood_pressure.action,rs$blood_lipids.action,rs$blood_glu_insu.action), # converts also factors to integer...
            1, # apply rowwise
            function(x) {
              # is level 1 / 2 / 3 reached or exceeded?
              compare_with_123 <- sapply(x, function(lev) lev >= 1:(length(perc_level)+1))
              #rownames(compare_with_123) <- lvl_name
              # which of the levels 1, 2, 3 are exceeded/reached at least 3 times
              levelcheck_123 <- apply(compare_with_123, 1, function(x) sum(x, na.rm=T)>=3)
              # return maximum level that is reached or exceeded at least three times
              max((1:(length(perc_level)+1))[levelcheck_123]) #max(ordered(1:(length(perc_level)+1), 1:(length(perc_level)+1), lvl_name)[levelcheck_123])
            })
    rs$overall.action <- ordered(rs$overall.action, 1:(length(perc_level)+1), lvl_name)
  }

  if (!is.null(filter)) rs <- rs[names(rs) %in% paste0(filter,".action")]

  if (append) rs <- cbind(df,rs)

  return(as.data.frame(rs))
}

#' Calculate the metabolic syndrome score for a whole data set
#'
#' The function calculates the metabolic syndrome score for a given data set.
#'
#' @param df data set including the study data with clinical parameters (sex, age, height, waist, sbp, dbp, trg, hdl, homa).
#' @param return_values A character vector with the possible values 'percentile' and 'z.score'.
#' @return A list with the values specified in `return_values`
#' @examples
#' my_data <- data.frame(
#' sex= c("m","f","f","m"),
#' age= c(6.5,5.8,5.2,5.5),
#' height= c(126,118,112,119),
#' waist= c(59.5,60.4,53.0,57.5),
#' sbp= c(109,92,94,99),
#' dbp= c(67.5,52.5,58.0,66.0),
#' trg= c(48,50,45,78),
#' hdl= c(38,52,61,45),
#' homa= c(0.65,1.59,1.00,0.99)
#' ) # insu, bmi, glu ?
#'
#' # return alls scores appended to input data
#' ScoreCalc(my_data, return_input=TRUE)
#'
#' # only return z-scores and MetS score
#' ScoreCalc(my_data, return_values = c("z.score","MetS"))
#'
#' # only return percentiles
#' ScoreCalc(my_data, return_values = "percentile")
#'
#' # only return MetS scores
#' ScoreCalc(my_data, return_values = "MetS")
#'
ScoreCalc <- function(df, return_input = F, return_values=c("percentile","z.score", "MetS", "action")) {

  # names of the variables to wich get_scores will be applied
  vars <- c("bmi", "glu", "hdl", "height", "homa", "insu", "trg", "waist", "sbp", "dbp")

  # define vector with parameters necessary to calculate the
  # MetS-score
  necparas <- c("waist", "homa", "sbp", "dbp", "trg", "hdl")

  # calculation of percentiles and z-scores

  # MetS score computation requires z-scores
  return_values_temp <- return_values
  if ("MetS" %in% return_values_temp && !"z.score" %in% return_values_temp)
    return_values_temp <- c(return_values_temp, "z.score")

  # Generate score columns with appropriate prefixes
  score_results <- lapply(vars, function(var_name) {
    if (var_name %in% names(df)) {
      scores <- as.data.frame(get_scores(var_name, df$sex, df$age, df$height, df[[var_name]], return_values_temp))
      # Prefix the column names
      colnames(scores) <- paste0(var_name, "_", colnames(scores))
      return(scores)
    } else {
      NULL
    }
  })

  # remove NULL
  score_results <- score_results[!sapply(score_results, is.null)]

  # bind columns
  score_results <- do.call(cbind, score_results)

  if ("MetS" %in% return_values) {
    tryCatch({
      MetS <- MetSScore(score_results)
      score_results <- cbind(score_results, MetS)

      # compute percentiles and/or z.scores for "+100"-shifted MetS score
      MetS_score_results <- as.data.frame(get_scores("MetS_shifted", df$sex, df$age, df$height, MetS+100, return_values))
      # Prefix the column names
      if (length(colnames(MetS_score_results)) > 0) colnames(MetS_score_results) <- paste0("MetS", "_", colnames(MetS_score_results))
      if (length(MetS_score_results)>0) score_results <- cbind(score_results, MetS_score_results)

      # remove z-scores if not requested...
      if (!"z.score" %in% return_values)
        score_results <- score_results[!endsWith(names(score_results), "_z.score")]
    },
    error = function(e) warning("MetS score could not be computed. Required variables: waist, homa, sbp, dbp, trg, hdl")
    )
  }

  if ("action" %in% return_values) {
    score_results <- action_levels(score_results, append = TRUE)
  }

  # if input requested
  if (return_input) score_results <- cbind(df, score_results)

  return(score_results)
}
