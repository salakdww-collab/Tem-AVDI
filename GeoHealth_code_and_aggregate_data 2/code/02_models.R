# Statistical models and aggregate estimates.
# Usage: Rscript code/02_models.R /path/to/private_inputs /path/to/private_work
# Use the same paths as step 01. The 03 script can run independently from the public CSVs.
# Model dependencies include survival, splines, dlnm (2.4.7), mvmeta, mixmeta, Matrix and MASS.

# Preserve the adopted DLNM basis, covariates, lag range and Breslow fitting method.
fit_primary_model <- function() {
  suppressPackageStartupMessages({
    library(dlnm)
    library(splines)
    library(survival)
  })

  maximum_lag <- primary_max_lag

  d <- as.data.table(read_rds(file.path(OUTPUT_ROOT, "data/derived/case_crossover.rds")))

  setorder(d, case_id, ref_date)

  loaded <- load(resolve_input(file.path(RAW_ROOT, "weather/municipal_daily.RData")))

  if (!identical(loaded, "data_meteorology_BRA_list")) stop("Unexpected weather object")

  all_weather_temperature <- unlist(lapply(data_meteorology_BRA_list, `[[`, "Tmean"), use.names = FALSE)

  weather_q <- quantile(all_weather_temperature,
    probs = c(0, 0.1, 0.75, 0.9, 1), na.rm = TRUE, names = FALSE,
    type = 7
  )

  rm(all_weather_temperature, data_meteorology_BRA_list)

  invisible(gc())

  exposure_bounds <- weather_q[c(1L, 5L)]
  exposure_knots <- weather_q[2:4]
  lag_knots <- logknots(maximum_lag, 3L)

  history <- as.matrix(d[, paste0("tmean_lag", 0:maximum_lag), with = FALSE])
  basis <- crossbasis(history, lag = c(0L, maximum_lag), argvar = list(
    fun = "ns", knots = exposure_knots,
    Boundary.knots = exposure_bounds
  ), arglag = list(fun = "ns", knots = lag_knots, Boundary.knots = c(
    0L,
    maximum_lag
  )))

  rm(history)

  basis_columns <- paste0("cb", seq_len(ncol(basis)))
  fit_data <- copy(d[, .(case, rh_lag03, holiday, case_id)])

  fit_data[, `:=`((basis_columns), as.data.table(unclass(basis)))]

  model_formula <- as.formula(paste("case ~", paste(c(
    basis_columns, "ns(rh_lag03, df = 3)", "holiday",
    "strata(case_id)"
  ), collapse = " + ")))

  started <- Sys.time()
  fit <- clogit(model_formula, data = fit_data, method = "breslow")
  basis_coef <- coef(fit)[basis_columns]
  basis_vcov <- vcov(fit)[basis_columns, basis_columns, drop = FALSE]
  model_information <- data.table(
    Scope = "Brazil", Cases = uniqueN(d$case_id), Rows = nrow(d), Municipalities = uniqueN(d$Code),
    BasisColumns = ncol(basis), Parameters = attr(logLik(fit), "df"), LogLikelihood = as.numeric(logLik(fit)),
    ElapsedSeconds = as.numeric(difftime(Sys.time(), started, units = "secs")), Converged = is.null(fit$fail)
  )

  retained <- list(
    coefficients = basis_coef, vcov = basis_vcov, information = model_information, exposure_knots = exposure_knots,
    exposure_bounds = exposure_bounds, lag_knots = lag_knots, maximum_lag = maximum_lag, formula = deparse(model_formula)
  )

  save_rds(retained, "results/models/national_association.rds")
}

# Individual subgroups: sex, age, education and race.
fit_individual_subgroups <- function() {
  suppressPackageStartupMessages({
    library(dlnm)
    library(splines)
    library(survival)
  })

  maximum_lag <- primary_max_lag
  report_lag <- primary_cumulative_lag
  temp_columns <- paste0("tmean_lag", 0:maximum_lag)

  matched <- as.data.table(read_rds(file.path(OUTPUT_ROOT, "data/derived/case_crossover.rds")))

  setorder(matched, case_id, ref_date)

  matched[, age_binary := fcase(age_group %in% c("<=49", "50-64"), "<65", age_group == ">=65", ">=65",
    default = NA_character_
  )]

  matched[, education_binary := fcase(education == ">=8 years", ">=8 years", education %in% c(
    "1-7 years",
    "Never educated"
  ), "<8 years", default = NA_character_)]

  matched[, race_binary := fcase(race == "White", "White", race %in% c("Mixed", "Black", "Asian", "Indigenous"),
    "non-White",
    default = NA_character_
  )]

  retained <- read_rds(file.path(OUTPUT_ROOT, "results/models/national_association.rds"))

  specifications <- data.table(
    subgroup = rep(c("Sex", "Age (years)", "Education (years)", "Race"), each = 2),
    variable = rep(c("sex", "age_binary", "education_binary", "race_binary"), each = 2),
    level = c("Male", "Female", "<65", ">=65", ">=8 years", "<8 years", "White", "non-White")
  )

  build_basis <- function(d) {
    temperature_history <- as.matrix(d[, ..temp_columns])
    crossbasis(temperature_history, lag = c(0L, maximum_lag), argvar = list(
      fun = "ns", knots = retained$exposure_knots,
      Boundary.knots = retained$exposure_bounds
    ), arglag = list(
      fun = "ns", knots = retained$lag_knots,
      Boundary.knots = c(0L, maximum_lag)
    ))
  }

  fit_one <- function(variable, level) {
    d <- matched[get(variable) == level & !is.na(get(variable))]
    setorder(d, case_id, ref_date)
    cases <- uniqueN(d[case == 1L, case_id])
    rows <- nrow(d)
    municipalities <- uniqueN(d[case == 1L, Code])
    basis <- build_basis(d)
    basis_columns <- paste0("cb_sub", seq_len(ncol(basis)))
    fit_data <- copy(d[, .(case, rh_lag03, holiday, case_id)])
    fit_data[, `:=`((basis_columns), as.data.table(unclass(basis)))]
    model_formula <- as.formula(paste("case ~", paste(c(
      basis_columns, "ns(rh_lag03, df = 3)", "holiday",
      "strata(case_id)"
    ), collapse = " + ")))
    fit <- clogit(model_formula, data = fit_data, method = "breslow")
    basis_coef <- coef(fit)[basis_columns]
    basis_vcov <- vcov(fit)[basis_columns, basis_columns, drop = FALSE]
    list(coefficients = basis_coef, vcov = basis_vcov, cases = cases, rows = rows,
         municipalities = municipalities, converged = is.null(fit$fail))
  }

  all_models <- vector("list", nrow(specifications))

  for (i in seq_len(nrow(specifications))) {
    all_models[[i]] <- fit_one(specifications$variable[i], specifications$level[i])
  }

  names(all_models) <- paste(specifications$subgroup, specifications$level, sep = ": ")

  save_rds(all_models, "results/models/individual_subgroups.rds")
}

# Contextual subgroups: baseline population-weighted median cutoffs.
fit_contextual_subgroups <- function() {
  suppressPackageStartupMessages({
    library(dlnm)
    library(splines)
    library(survival)
  })

  maximum_lag <- primary_max_lag
  report_lag <- primary_cumulative_lag
  temp_columns <- paste0("tmean_lag", 0:maximum_lag)
  municipal <- municipal_indicators()
  spec <- read_csv(file.path(OUTPUT_ROOT, "results/estimates/supplement/TableS1.csv"))

  matched <- as.data.table(read_rds(file.path(OUTPUT_ROOT, "data/derived/case_crossover.rds")))

  setorder(matched, case_id, ref_date)

  matched[, Code := sprintf("%06d", as.integer(Code))]

  matched[municipal, `:=`(
    income_per_capita = i.income_per_capita, hospital_beds = i.hospital_beds, air_conditioning = i.air_conditioning,
    automobile = i.automobile
  ), on = "Code"]

  retained <- read_rds(file.path(OUTPUT_ROOT, "results/models/national_association.rds"))

  build_basis <- function(d) {
    crossbasis(as.matrix(d[, ..temp_columns]), lag = c(0L, maximum_lag), argvar = list(
      fun = "ns",
      knots = retained$exposure_knots, Boundary.knots = retained$exposure_bounds
    ), arglag = list(
      fun = "ns",
      knots = retained$lag_knots, Boundary.knots = c(0L, maximum_lag)
    ))
  }

  fit_group <- function(d, sp, group_name) {
    setorder(d, case_id, ref_date)
    basis <- build_basis(d)
    basis_cols <- paste0("cb_context_", seq_len(ncol(basis)))
    fit_data <- copy(d[, .(case, rh_lag03, holiday, case_id)])
    fit_data[, `:=`((basis_cols), as.data.table(unclass(basis)))]
    model_formula <- as.formula(paste("case ~", paste(c(
      basis_cols, "ns(rh_lag03, df = 3)", "holiday",
      "strata(case_id)"
    ), collapse = " + ")))
    fit <- clogit(model_formula, data = fit_data, method = "breslow")
    basis_coef <- coef(fit)[basis_cols]
    basis_vcov <- vcov(fit)[basis_cols, basis_cols, drop = FALSE]
    list(
      coefficients = basis_coef, vcov = basis_vcov, exposure_knots = retained$exposure_knots,
      exposure_bounds = retained$exposure_bounds, lag_knots = retained$lag_knots, maximum_lag = maximum_lag,
      report_lag = report_lag, cutoff_method = sp$cutoff_method, cutoff = sp$cutoff, group = group_name,
      cases = uniqueN(d$case_id)
    )
  }

  all_models <- list()

  for (i in seq_len(nrow(spec))) {
    sp <- spec[i]
    value_column <- sp$value_column
    d <- copy(matched[!is.na(get(value_column))])
    d[, group := fifelse(get(value_column) <= sp$cutoff, "Low", "High")]
    fits <- list()
    for (g in c("Low", "High")) {
      fits[[g]] <- fit_group(copy(d[group == g]), sp, g)
    }
    all_models[[sp$modifier]] <- fits
  }

  save_rds(all_models, "results/models/contextual_subgroups.rds")
}

# Manuscript Table1: characteristics of the complete model sample.
make_table1 <- function() {
  suppressPackageStartupMessages({
  })

  cases <- analysis_deaths()

  cases[, Year := as.integer(format(death_date, "%Y"))]

  cases[, sex_table := fcase(SEXO == 1, "Male", SEXO == 2, "Female", default = "Unknown")]

  cases[, age_table := fcase(age_years <= 49, "<=49 years", age_years >= 50 & age_years <= 64, "50-64 years",
    age_years >= 65, ">=65 years",
    default = "Unknown"
  )]

  cases[, education_table := fcase(ESC == 1, "No formal education", ESC %in% c(2, 3), "1-7 years", ESC %in%
    c(4, 5), ">=8 years", default = "Unknown")]

  cases[, race_table := fcase(RACACOR == 1, "White", RACACOR == 2, "Black", RACACOR == 3, "Asian", RACACOR ==
    4, "Mixed", RACACOR == 5, "Indigenous", default = "Unknown")]

  cases[, marital_table := fcase(ESTCIV == 1, "Single", ESTCIV %in% c(2, 5), "Married/consensual union",
    ESTCIV == 3, "Widowed", ESTCIV == 4, "Separated/divorced",
    default = "Unknown"
  )]

  cases[, place_table := fcase(LOCOCOR == 1, "Hospital", LOCOCOR == 2, "Other healthcare facility",
    LOCOCOR == 3, "Home", LOCOCOR == 4, "Public road", LOCOCOR == 5, "Other",
    default = "Unknown"
  )]

  cases[, region_table := fcase(region == "CenterWest", "Center-West", !is.na(region), region, default = "Unknown")]

  n_total <- nrow(cases)
  distribution_rows <- function(variable, categories, characteristic) {
    counts <- cases[, .(Deaths = .N), by = .(Category = get(variable))]
    result <- counts[data.table(Category = categories), on = "Category"]
    result[is.na(Deaths), Deaths := 0L]
    result[, `:=`(Characteristic = characteristic, Percent = 100 * Deaths / n_total)]
    result[, Display := sprintf("%s (%.1f)", format(Deaths, big.mark = ",", trim = TRUE), Percent)]
    result[, category_order := match(Category, categories)]
    setorder(result, category_order)
    result[, .(Characteristic, Category, Deaths, Percent, Display)]
  }

  table1 <- rbindlist(list(
    data.table(
      Characteristic = "Overall", Category = "All included deaths", Deaths = n_total,
      Percent = 100, Display = sprintf("%s (100.0)", format(n_total, big.mark = ","))
    ), distribution_rows(
      "sex_table",
      c("Male", "Female", "Unknown"), "Sex"
    ), distribution_rows("age_table", c(
      "<=49 years", "50-64 years",
      ">=65 years", "Unknown"
    ), "Age"), distribution_rows("education_table", c(
      "No formal education", "1-7 years",
      ">=8 years", "Unknown"
    ), "Education"), distribution_rows("race_table", c(
      "White", "Black", "Mixed",
      "Asian", "Indigenous", "Unknown"
    ), "Race/skin colour"), distribution_rows("marital_table", c(
      "Single",
      "Married/consensual union", "Widowed", "Separated/divorced", "Unknown"
    ), "Marital status"), distribution_rows(
      "place_table",
      c("Hospital", "Other healthcare facility", "Home", "Public road", "Other", "Unknown"), "Place of death"
    ),
    distribution_rows(
      "region_table", c("North", "Northeast", "Center-West", "Southeast", "South", "Unknown"),
      "Region"
    )
  ), use.names = TRUE)

  # Match the manuscript display: omit Unknown rows, retain all deaths in the denominator.
  write_csv(table1[Category != "Unknown"], "results/estimates/main/Table1.csv")
}

# Table S2: clinical context recorded on death certificates.
make_table_s2 <- function() {

  input_file <- file.path(OUTPUT_ROOT, "data/derived/death_certificates.rds")

  deaths <- as.data.table(read_rds(input_file))[SOURCE_YEAR >= 2000L & SOURCE_YEAR <= 2019L & trimws(as.character(CAUSABAS)) ==
    "K550"]

  registered <- as.data.table(read_rds(file.path(OUTPUT_ROOT, "data/derived/deaths.rds")))
  included <- analysis_deaths()
  cert <- copy(deaths)
  death_link <- copy(registered)
  for (v in c("CODMUNRES", "CODMUNOCOR")) {
    death_link[, (v) := substr(as.character(get(v)), 1L, 6L)]
    cert[, (v) := substr(as.character(get(v)), 1L, 6L)]
  }
  death_link[, OCUP := sub("^0+", "", as.character(OCUP))]
  cert[, OCUP := sub("^0+", "", as.character(OCUP))]
  match_vars <- c("DTOBITO", "CODMUNRES", "CODMUNOCOR", "DTNASC", "IDADE", "SEXO", "RACACOR", "ESTCIV", "ESC", "OCUP", "LOCOCOR", "CAUSABAS")
  make_key <- function(z) do.call(paste, c(lapply(z[, ..match_vars], function(v) fifelse(is.na(v), "<NA>", trimws(as.character(v)))), sep = "|"))
  death_link[, record_key := make_key(death_link)]
  cert[, record_key := make_key(cert)]
  stopifnot(!anyDuplicated(death_link$record_key), !anyDuplicated(cert$record_key), setequal(death_link$record_key, cert$record_key))
  cert_included <- cert[record_key %in% death_link[case_id %in% included$case_id, record_key]]
  stopifnot(nrow(cert_included) == nrow(included))
  deaths <- cert_included
  n_deaths <- nrow(deaths)
  parse_icd3 <- function(values) {
    values <- toupper(gsub("[[:space:]]+", "", ifelse(is.na(values), "", as.character(values))))
    lapply(strsplit(values, "*", fixed = TRUE), function(z) {
      z <- gsub("[^A-Z0-9]", "", z)
      unique(substr(z[nchar(z) >= 3L], 1L, 3L))
    })
  }

  count_regex <- function(code_lists, regex) {
    sum(vapply(code_lists, function(z) any(grepl(regex, z, perl = TRUE)), logical(1L)))
  }

  line_a_codes <- parse_icd3(deaths$LINHAA)
  part_ii_codes <- parse_icd3(deaths$LINHAII)
  part_ii_non_k55 <- lapply(part_ii_codes, function(z) z[z != "K55"])
  immediate_spec <- data.table(condition = c(
    "Sepsis", "Vascular disorders of intestine", "Other general symptoms and signs",
    "Shock", "Respiratory failure", "Other circulatory or respiratory symptoms and signs", "Complications of procedures",
    "Peritonitis", "Abdominal and pelvic pain", "Fluid, electrolyte and acid-base disorders"
  ), icd10_definition = c(
    "A40-A41",
    "K55", "R68", "R57", "J96", "R09", "T81", "K65", "R10", "E87"
  ), regex = c(
    "^A4[01]$", "^K55$", "^R68$",
    "^R57$", "^J96$", "^R09$", "^T81$", "^K65$", "^R10$", "^E87$"
  ))

  immediate_spec[, deaths := vapply(regex, count_regex, integer(1L), code_lists = line_a_codes)]

  immediate_spec[, percent := round(100 * deaths / n_deaths, 1L)]

  setorder(immediate_spec, -deaths, condition)

  immediate_spec[, `:=`(
    panel = "A. Immediate causes of death", row_type = "condition", parent_group = NA_character_,
    display_order = seq_len(.N)
  )]

  chapter_spec <- data.table(chapter = c("I", "E", "N", "J", "K"), condition = c(
    "Diseases of the circulatory system",
    "Endocrine, nutritional and metabolic diseases", "Diseases of the genitourinary system", "Diseases of the respiratory system",
    "Diseases of the digestive system"
  ), icd10_definition = c(
    "I00-I99", "E00-E90", "N00-N99", "J00-J99",
    "K00-K93"
  ))

  chapter_spec[, regex := paste0("^", chapter)]

  chapter_spec[, deaths := vapply(regex, count_regex, integer(1L), code_lists = part_ii_non_k55)]

  chapter_spec[, percent := round(100 * deaths / n_deaths, 1L)]

  setorder(chapter_spec, -deaths)

  child_spec <- data.table(
    chapter = c(
      "I", "I", "I", "E", "E", "E", "N", "N", "N", "J", "J", "J", "K",
      "K", "K"
    ), condition = c(
      "Hypertensive diseases", "Ischaemic heart diseases", "Heart failure", "Diabetes mellitus",
      "Obesity", "Unspecified protein-energy malnutrition", "Acute kidney failure", "Chronic kidney disease",
      "Unspecified kidney failure", "Other chronic obstructive pulmonary disease", "Pneumonia, organism unspecified",
      "Respiratory failure", "Fibrosis and cirrhosis of liver", "Other diseases of liver", "Other diseases of digestive system"
    ),
    icd10_definition = c(
      "I10-I15", "I20-I25", "I50", "E10-E14", "E66", "E46", "N17", "N18", "N19", "J44",
      "J18", "J96", "K74", "K76", "K92"
    ), regex = c(
      "^I1[0-5]$", "^I2[0-5]$", "^I50$", "^E1[0-4]$",
      "^E66$", "^E46$", "^N17$", "^N18$", "^N19$", "^J44$", "^J18$", "^J96$", "^K74$", "^K76$", "^K92$"
    )
  )

  child_spec[, deaths := vapply(regex, count_regex, integer(1L), code_lists = part_ii_non_k55)]

  child_spec[, percent := round(100 * deaths / n_deaths, 1L)]

  setorder(child_spec, chapter, -deaths, condition)

  child_spec[, child_rank := seq_len(.N), by = chapter]

  stopifnot(all(child_spec$child_rank <= 3L))

  contributing_rows <- rbindlist(lapply(seq_len(nrow(chapter_spec)), function(i) {
    parent <- chapter_spec[i]
    children <- child_spec[chapter == parent$chapter][order(child_rank)]
    rbindlist(
      list(data.table(
        panel = "B. Conditions contributing to death", row_type = "parent", parent_group = NA_character_,
        condition = parent$condition, icd10_definition = parent$icd10_definition, deaths = parent$deaths,
        percent = parent$percent, display_order = (i - 1L) * 4L + 1L
      ), data.table(
        panel = "B. Conditions contributing to death",
        row_type = "child", parent_group = parent$condition, condition = children$condition, icd10_definition = children$icd10_definition,
        deaths = children$deaths, percent = children$percent, display_order = (i - 1L) * 4L + 1L + children$child_rank
      )),
      use.names = TRUE
    )
  }))

  table_output <- rbindlist(list(immediate_spec[, .(
    panel, row_type, parent_group, condition, icd10_definition,
    deaths, percent, display_order
  )], contributing_rows), use.names = TRUE)

  table_output[, panel_order := fifelse(panel == "A. Immediate causes of death", 1L, 2L)]

  setorder(table_output, panel_order, display_order)

  setcolorder(table_output, c(
    "panel", "row_type", "parent_group", "condition", "icd10_definition", "deaths",
    "percent", "panel_order", "display_order"
  ))

  write_csv(table_output, "results/estimates/supplement/TableS2.csv", na = "")
}

# Table S3: temperature and relative humidity distributions.
make_table_s3 <- function() {

  deaths <- as.data.table(read_rds(file.path(OUTPUT_ROOT, "data/derived/deaths.rds")))
  included_codes <- unique(as.character(deaths$Code))
  loaded <- load(resolve_input(file.path(RAW_ROOT, "weather/municipal_daily.RData")))

  if (!identical(loaded, "data_meteorology_BRA_list")) {
    stop("Unexpected weather source object")
  }

  weather <- data_meteorology_BRA_list

  rm(data_meteorology_BRA_list)

  included_codes <- intersect(included_codes, names(weather))
  summarize_values <- function(codes, variable) {
    values <- unlist(lapply(weather[codes], function(city) city[[variable]]), use.names = FALSE)
    result <- c(
      Mean = mean(values, na.rm = TRUE), SD = sd(values, na.rm = TRUE), Min = min(values, na.rm = TRUE),
      P1 = unname(quantile(values, 0.01, na.rm = TRUE)), Median = median(values, na.rm = TRUE), P99 = unname(quantile(values,
        0.99,
        na.rm = TRUE
      )), Max = max(values, na.rm = TRUE)
    )
    rm(values)
    invisible(gc())
    result
  }

  regions <- c("Overall", "North", "Northeast", "CenterWest", "Southeast", "South")

  table_s1 <- rbindlist(lapply(c("Tmean", "RH"), function(variable) {
    rbindlist(lapply(regions, function(region) {
      codes <- if (region == "Overall") {
        included_codes
      } else {
        included_codes[unname(state_to_region[substr(included_codes, 1L, 2L)]) == region]
      }
      values <- summarize_values(codes, variable)
      data.table(
        Meteorology = ifelse(variable == "Tmean", "Temperature (C)", "Relative humidity (%)"),
        Region = region, Mean = values["Mean"], SD = values["SD"], Min = values["Min"], P1 = values["P1"],
        Median = values["Median"], P99 = values["P99"], Max = values["Max"]
      )
    }))
  }))

  write_csv(table_s1, "results/estimates/supplement/TableS3.csv")
}

# Reproduce the adopted occurrence-assigned model on the same 52,565 strata.
fit_occurrence_sensitivity <- function() {
  suppressPackageStartupMessages({library(dlnm);library(splines);library(survival)})
  cc <- as.data.table(read_rds(file.path(OUTPUT_ROOT, "data/derived/case_crossover.rds")))
  weather_env <- new.env()
  load(resolve_input(file.path(RAW_ROOT, "weather/municipal_daily.RData")), envir = weather_env)
  weather <- weather_env$data_meteorology_BRA_list
  rm(weather_env)
  fill_exposure <- function(d, code_vector) {
    n <- nrow(d)
    th <- matrix(NA_real_, n, 22L, dimnames = list(NULL, paste0("tmean_lag", 0:21)))
    rh <- matrix(NA_real_, n, 4L, dimnames = list(NULL, paste0("rh_lag", 0:3)))
    groups <- split(seq_len(n), code_vector)
    for (code in names(groups)) {
      rows <- groups[[code]]
      city <- weather[[code]]
      if (is.null(city) || !nrow(city)) next
      city_dates <- as.Date(city$date)
      target <- as.Date(d$ref_date[rows])
      for (lag in 0:21) th[rows, lag + 1L] <- city$Tmean[match(target - lag, city_dates)]
      for (lag in 0:3) rh[rows, lag + 1L] <- city$RH[match(target - lag, city_dates)]
    }
    ans <- as.data.table(th)
    ans[, rh_lag03 := rowMeans(rh, na.rm = FALSE)]
    ans
  }
  
  # Occurrence-municipality exposure assignment on the identical 52,565 strata.
  occ_code <- substr(as.character(cc$CODMUNOCOR), 1L, 6L)
  occ_exp <- fill_exposure(cc, occ_code)
  occ_complete <- complete.cases(occ_exp)
  if (!all(occ_complete)) stop("Occurrence assignment is not complete for the fixed analytic rows")
  
  occ <- copy(cc)
  occ[, (paste0("tmean_lag", 0:21)) := occ_exp[, paste0("tmean_lag", 0:21), with = FALSE]]
  occ[, rh_lag03 := occ_exp$rh_lag03]
  occ[, occurrence_code := occ_code]
  setorder(occ, case_id, ref_date)
  
  main_model <- read_rds(file.path(OUTPUT_ROOT, "results/models/national_association.rds"))
  main_result <- read_csv(file.path(OUTPUT_ROOT, "results/estimates/main/national_association.csv"))
  comparison_temperatures <- main_result[match(c("P1 cold", "P99 heat"), Exposure), TemperatureC]
  primary_mmt <- unique(main_result$ReferenceC)
  history <- as.matrix(occ[, paste0("tmean_lag", 0:21), with = FALSE])
  basis <- crossbasis(history, lag = c(0L, 21L), argvar = list(fun = "ns", knots = main_model$exposure_knots, Boundary.knots = main_model$exposure_bounds), arglag = list(fun = "ns", knots = main_model$lag_knots, Boundary.knots = c(0L, 21L)))
  basis_names <- paste0("occ_cb", seq_len(ncol(basis)))
  fit_data <- copy(occ[, .(case, rh_lag03, holiday, case_id)])
  fit_data[, (basis_names) := as.data.table(unclass(basis))]
  form <- as.formula(paste("case ~", paste(c(basis_names, "ns(rh_lag03, df = 3)", "holiday", "strata(case_id)"), collapse = " + ")))
  fit <- clogit(form, data = fit_data, method = "breslow")
  beta <- coef(fit)[basis_names]
  sigma <- vcov(fit)[basis_names, basis_names, drop = FALSE]
  pred <- crosspred(basis, coef = beta, vcov = sigma, at = comparison_temperatures, cen = primary_mmt, lag = c(0L, 5L), model.link = "log")
  case_temp <- occ[case == 1L, tmean_lag0]
  q <- quantile(case_temp, c(.01, .99), type = 7)
  mmt_grid <- sort(unique(c(seq(round(q[1], 2), round(q[2], 2), by = .02), primary_mmt)))
  prov <- crosspred(basis, coef = beta, vcov = sigma, at = mmt_grid, cen = primary_mmt, lag = c(0L, 5L), model.link = "log")
  own_mmt <- prov$predvar[which.min(prov$allRRfit)]
  occ_result <- data.table(
    Exposure = c("P1 cold", "P99 heat"), TemperatureC = comparison_temperatures,
    ReferenceC = primary_mmt, ModelSpecificMMT = own_mmt,
    OR = as.numeric(pred$allRRfit), Lower95 = as.numeric(pred$allRRlow), Upper95 = as.numeric(pred$allRRhigh),
    Cases = uniqueN(occ$case_id), Rows = nrow(occ), Municipalities = uniqueN(occ$occurrence_code),
    AIC = AIC(fit), LogLikelihood = as.numeric(logLik(fit))
  )
  occ_result[main_result[, .(Exposure, PrimaryOR = OR, PrimaryLower95 = Lower95, PrimaryUpper95 = Upper95)], on = "Exposure", `:=`(PrimaryOR = i.PrimaryOR, PrimaryLower95 = i.PrimaryLower95, PrimaryUpper95 = i.PrimaryUpper95)]
  occ_result[, PercentDifferenceFromPrimaryOR := 100 * (OR / PrimaryOR - 1)]
  write_csv(occ_result, "results/estimates/sensitivity/place_of_death.csv")
  save_rds(list(coefficients = beta, vcov = sigma, cases = uniqueN(occ$case_id), rows = nrow(occ), formula = deparse(form)), "results/models/place_of_death.rds")
}

# Exposure assignment and temperature/lag spline specifications.
fit_exposure_sensitivity <- function() {
  suppressPackageStartupMessages({
    library(dlnm)
    library(splines)
    library(survival)
  })

  maximum_lag <- primary_max_lag
  report_lag <- primary_cumulative_lag

  `%||%` <- function(x, y) if (is.null(x)) y else x

  main_data <- as.data.table(read_rds(file.path(OUTPUT_ROOT, "data/derived/case_crossover.rds")))

  setorder(main_data, case_id, ref_date)

  primary_summary <- read_csv(file.path(OUTPUT_ROOT, "results/estimates/main/national_association.csv"))

  primary_mmt <- unique(primary_summary$ReferenceC)
  primary_p1 <- primary_summary[Exposure == "P1 cold", TemperatureC]
  primary_p99 <- primary_summary[Exposure == "P99 heat", TemperatureC]

  stopifnot(length(primary_mmt) == 1L, length(primary_p1) == 1L, length(primary_p99) == 1L)

  comparison_temperatures <- c(primary_p1, primary_p99)

  loaded <- load(resolve_input(file.path(RAW_ROOT, "weather/municipal_daily.RData")))

  if (!identical(loaded, "data_meteorology_BRA_list")) {
    stop("Unexpected weather object")
  }

  all_temperature <- unlist(lapply(data_meteorology_BRA_list, `[[`, "Tmean"), use.names = FALSE)

  temperature_q <- quantile(all_temperature,
    probs = c(0, 0.1, 0.5, 0.75, 0.9, 1), na.rm = TRUE, names = FALSE,
    type = 7
  )

  rm(all_temperature, data_meteorology_BRA_list)

  invisible(gc())

  exposure_bounds <- temperature_q[c(1L, 6L)]
  exposure_knots_primary <- temperature_q[c(2L, 4L, 5L)]
  exposure_knots_2 <- temperature_q[c(2L, 5L)]
  exposure_knots_mid50 <- temperature_q[c(2L, 3L, 5L)]
  exposure_knots_4 <- temperature_q[c(2L, 3L, 4L, 5L)]
  deaths <- as.data.table(read_rds(file.path(OUTPUT_ROOT, "data/derived/deaths.rds")))

  deaths[, `:=`(residence_code = substr(as.character(CODMUNRES), 1L, 6L), occurrence_code = substr(
    as.character(CODMUNOCOR),
    1L, 6L
  ), analysis_code = as.character(Code))]

  fallback <- deaths[analysis_code != residence_code & analysis_code == occurrence_code]

  fallback_ids <- intersect(fallback$case_id, unique(main_data$case_id))

  stopifnot(nrow(fallback) == 65L, length(fallback_ids) == 65L)

  residence_only_data <- main_data[!case_id %in% fallback_ids]

  stopifnot(uniqueN(residence_only_data$case_id) == 52500L)

  fit_specification <- function(
    d, model_id, family, specification, exposure_knots = exposure_knots_primary,
    lag_knot_count = 3L, humidity_adjustment = TRUE
  ) {
    setorder(d, case_id, ref_date)
    history <- as.matrix(d[, paste0("tmean_lag", 0:maximum_lag), with = FALSE])
    lag_knots <- logknots(maximum_lag, lag_knot_count)
    basis <- crossbasis(history, lag = c(0L, maximum_lag), argvar = list(
      fun = "ns", knots = exposure_knots,
      Boundary.knots = exposure_bounds
    ), arglag = list(fun = "ns", knots = lag_knots, Boundary.knots = c(
      0L,
      maximum_lag
    )))
    rm(history)
    basis_names <- paste0("sens_cb", seq_len(ncol(basis)))
    fit_data <- copy(d[, .(case, rh_lag03, holiday, case_id)])
    fit_data[, `:=`((basis_names), as.data.table(unclass(basis)))]
    model_terms <- c(
      basis_names, if (humidity_adjustment) "ns(rh_lag03, df = 3)" else character(), "holiday",
      "strata(case_id)"
    )
    model_formula <- as.formula(paste("case ~", paste(model_terms, collapse = " + ")))
    message(sprintf("Fitting %s: %s cases", specification, format(uniqueN(d$case_id), big.mark = ",")))
    started <- Sys.time()
    fit <- clogit(model_formula, data = fit_data, method = "breslow")
    beta <- coef(fit)[basis_names]
    sigma <- vcov(fit)[basis_names, basis_names, drop = FALSE]
    case_temperature <- d[case == 1L, tmean_lag0]
    q <- quantile(case_temperature, c(0.01, 0.99), type = 7, na.rm = TRUE)
    mmt_grid <- seq(round(q[1L], 2), round(q[2L], 2), by = 0.02)
    mmt_grid <- sort(unique(c(mmt_grid, primary_mmt)))
    provisional <- crosspred(basis, coef = beta, vcov = sigma, at = mmt_grid, cen = primary_mmt, lag = c(
      0L,
      report_lag
    ), model.link = "log")
    own_mmt <- provisional$predvar[which.min(provisional$allRRfit)]
    pred <- crosspred(basis,
      coef = beta, vcov = sigma, at = comparison_temperatures, cen = primary_mmt,
      lag = c(0L, report_lag), model.link = "log"
    )
    model_row <- data.table(
      ModelID = model_id, Family = family, Specification = specification, Cases = uniqueN(d$case_id),
      Rows = nrow(d), ExcludedCases = uniqueN(main_data$case_id) - uniqueN(d$case_id), MaximumLag = maximum_lag,
      ReportedCumulativeLag = paste0("0-", report_lag), ExposureInternalKnots = length(exposure_knots),
      ExposureKnotPercentiles = if (identical(exposure_knots, exposure_knots_2)) {
        "P10;P90"
      } else if (identical(exposure_knots, exposure_knots_mid50)) {
        "P10;P50;P90"
      } else if (identical(exposure_knots, exposure_knots_4)) {
        "P10;P50;P75;P90"
      } else {
        "P10;P75;P90"
      }, ExposureKnotValuesC = paste(round(exposure_knots, 4), collapse = ";"), LagInternalKnots = lag_knot_count,
      LagKnotValuesDays = paste(round(lag_knots, 4), collapse = ";"), HumidityAdjustment = if (humidity_adjustment) {
        "Natural spline of lag 0-3 mean RH (3 df)"
      } else {
        "None"
      }, ReferenceTemperatureC = primary_mmt, ModelSpecificMMT_C = own_mmt, Parameters = attr(
        logLik(fit),
        "df"
      ), LogLikelihood = as.numeric(logLik(fit)), AIC = AIC(fit), ElapsedSeconds = as.numeric(difftime(Sys.time(),
        started,
        units = "secs"
      )), Converged = is.null(fit$fail)
    )
    result <- data.table(
      ModelID = model_id, Family = family, Specification = specification, Exposure = c(
        "P1 cold",
        "P99 heat"
      ), TemperatureC = comparison_temperatures, ReferenceTemperatureC = primary_mmt, CumulativeLag = paste0(
        "0-",
        report_lag
      ), OR = as.numeric(pred$allRRfit), Lower95 = as.numeric(pred$allRRlow), Upper95 = as.numeric(pred$allRRhigh),
      Cases = uniqueN(d$case_id), Rows = nrow(d), ExcludedCases = uniqueN(main_data$case_id) - uniqueN(d$case_id),
      Converged = is.null(fit$fail)
    )
    list(summary = model_row, results = result, model = list(
      coefficients = beta, vcov = sigma, formula = deparse(model_formula),
      exposure_knots = exposure_knots, exposure_bounds = exposure_bounds, lag_knots = lag_knots, maximum_lag = maximum_lag,
      cases = uniqueN(d$case_id), rows = nrow(d)
    ))
  }

  specifications <- list(
    list(data = main_data, id = "primary", family = "Primary model", label = "Primary model"),
    list(data = residence_only_data, id = "residence_only", family = "Exposure-address assignment", label = "Valid residence municipality only"),
    list(
      data = main_data, id = "exposure_2_knots", family = "Temperature exposure basis", label = "Exposure knots: P10 and P90",
      exposure_knots = exposure_knots_2
    ), list(
      data = main_data, id = "exposure_4_knots", family = "Temperature exposure basis",
      label = "Exposure knots: P10, P50, P75, and P90", exposure_knots = exposure_knots_4
    ), list(
      data = main_data,
      id = "exposure_mid50_knots", family = "Temperature exposure basis", label = "Exposure knots: P10, P50, and P90",
      exposure_knots = exposure_knots_mid50
    ), list(
      data = main_data, id = "lag_2_knots", family = "Lag basis",
      label = "Lag basis: 2 log-spaced internal knots", lag_knot_count = 2L
    ), list(
      data = main_data,
      id = "lag_4_knots", family = "Lag basis", label = "Lag basis: 4 log-spaced internal knots", lag_knot_count = 4L
    ),
    list(
      data = main_data, id = "no_humidity", family = "Humidity adjustment", label = "No relative-humidity adjustment",
      humidity_adjustment = FALSE
    )
  )

  specifications <- Filter(function(x) {
    x$id %in% c(
      "residence_only", "exposure_2_knots", "exposure_4_knots",
      "lag_2_knots", "lag_4_knots"
    )
  }, specifications)

  fits <- lapply(specifications, function(s) {
    do.call(fit_specification, list(
      d = copy(s$data), model_id = s$id, family = s$family, specification = s$label,
      exposure_knots = s$exposure_knots %||% exposure_knots_primary, lag_knot_count = s$lag_knot_count %||%
        3L, humidity_adjustment = s$humidity_adjustment %||% TRUE
    ))
  })

  model_summary <- rbindlist(lapply(fits, `[[`, "summary"), fill = TRUE)
  results <- rbindlist(lapply(fits, `[[`, "results"), fill = TRUE)
  models <- setNames(lapply(fits, `[[`, "model"), model_summary$ModelID)

  if (!all(model_summary$Converged)) {
    stop("At least one selected sensitivity model did not converge")
  }

  primary_effect <- read_csv(file.path(OUTPUT_ROOT, "results/estimates/main/national_association.csv"))[
    ,
    .(Exposure, PrimaryOR = OR, PrimaryLower95 = Lower95, PrimaryUpper95 = Upper95)
  ]

  results[primary_effect, `:=`(PrimaryOR = i.PrimaryOR, PrimaryLower95 = i.PrimaryLower95, PrimaryUpper95 = i.PrimaryUpper95),
    on = "Exposure"
  ]

  results[, PercentDifferenceFromPrimaryOR := 100 * (OR / PrimaryOR - 1)]

  write_csv(results, "results/estimates/sensitivity/exposure_and_splines.csv")

  save_rds(models, "results/models/exposure_and_splines.rds")
}

# Relative humidity adjustment with 4 df.
fit_humidity_sensitivity <- function() {
  suppressPackageStartupMessages({
    library(dlnm)
    library(splines)
    library(survival)
  })

  main_data <- as.data.table(read_rds(file.path(OUTPUT_ROOT, "data/derived/case_crossover.rds")))

  setorder(main_data, case_id, ref_date)

  primary_summary <- read_csv(file.path(OUTPUT_ROOT, "results/estimates/main/national_association.csv"))

  comparison_temperatures <- primary_summary[, TemperatureC]
  loaded <- load(resolve_input(file.path(RAW_ROOT, "weather/municipal_daily.RData")))

  if (!identical(loaded, "data_meteorology_BRA_list")) stop("Unexpected weather object")

  all_temperature <- unlist(lapply(data_meteorology_BRA_list, `[[`, "Tmean"), use.names = FALSE)

  temperature_probability_grid <- quantile(all_temperature,
    probs = c(0, 0.1, 0.5, 0.75, 0.9, 1), na.rm = TRUE,
    type = 7
  )

  rm(all_temperature)

  exposure_bounds <- temperature_probability_grid[c(1L, 6L)]

  exposure_sets <- list(`2 internal knots (P10/P90; df 3)` = temperature_probability_grid[c(2L, 5L)], `Primary 3 knots (P10/P75/P90; df 4)` = temperature_probability_grid[c(
    2L,
    4L, 5L
  )], `4 internal knots (P10/P50/P75/P90; df 5)` = temperature_probability_grid[c(
    2L, 3L, 4L,
    5L
  )])

  fit_specification <- function(
    d, family, specification, max_lag = 21L, exposure_knots = exposure_sets[[2L]],
    lag_knot_count = 3L, humidity = "ns3"
  ) {
    setorder(d, case_id, ref_date)
    history <- as.matrix(d[, paste0("tmean_lag", 0:max_lag), with = FALSE])
    basis <- crossbasis(history, lag = c(0L, max_lag), argvar = list(
      fun = "ns", knots = exposure_knots,
      Boundary.knots = exposure_bounds
    ), arglag = list(
      fun = "ns", knots = logknots(max_lag, lag_knot_count),
      Boundary.knots = c(0L, max_lag)
    ))
    basis_names <- paste0("sens_cb", seq_len(ncol(basis)))
    fit_data <- copy(d[, .(case, rh_lag03, holiday, case_id)])
    fit_data[, `:=`((basis_names), as.data.table(unclass(basis)))]
    humidity_term <- switch(humidity,
      none = character(),
      linear = "rh_lag03",
      ns3 = "ns(rh_lag03, df = 3)",
      ns4 = "ns(rh_lag03, df = 4)"
    )
    formula <- as.formula(paste("case ~", paste(c(basis_names, humidity_term, "holiday", "strata(case_id)"),
      collapse = " + "
    )))
    fit <- clogit(formula, data = fit_data, method = "breslow")
    beta <- coef(fit)[basis_names]
    sigma <- vcov(fit)[basis_names, basis_names, drop = FALSE]
    case_temperature <- d[case == 1L, tmean_lag0]
    q <- quantile(case_temperature, c(0.01, 0.5, 0.99), type = 7)
    grid <- seq(round(q[1L], 2), round(q[3L], 2), by = 0.03)
    preliminary <- crosspred(basis, coef = beta, vcov = sigma, at = grid, cen = round(q[2L], 2), lag = c(
      0L,
      min(5L, max_lag)
    ), model.link = "log")
    mmt <- preliminary$predvar[which.min(preliminary$allRRfit)]
    pred <- crosspred(basis, coef = beta, vcov = sigma, at = comparison_temperatures, cen = mmt, lag = c(
      0L,
      min(5L, max_lag)
    ), model.link = "log")
    data.table(
      Family = family, Specification = specification, MaximumLag = max_lag, ReportedCumulativeLag = "0-5",
      Cases = uniqueN(d$case_id), Rows = nrow(d), ExposureBasisColumns = length(exposure_knots) + 1L,
      LagInternalKnots = lag_knot_count, HumidityAdjustment = humidity, MMT_C = mmt, MMT_Percentile = mean(case_temperature <=
        mmt) * 100, ColdTemperatureC = comparison_temperatures[1L], ColdOR = as.numeric(pred$allRRfit[1L]),
      ColdLower95 = as.numeric(pred$allRRlow[1L]), ColdUpper95 = as.numeric(pred$allRRhigh[1L]), HeatTemperatureC = comparison_temperatures[2L],
      HeatOR = as.numeric(pred$allRRfit[2L]), HeatLower95 = as.numeric(pred$allRRlow[2L]), HeatUpper95 = as.numeric(pred$allRRhigh[2L]),
      Converged = is.null(fit$fail)
    )
  }

  sensitivity <- fit_specification(copy(main_data), "Humidity adjustment", "ns(df=4) lag 0-3 RH", humidity = "ns4")

  write_csv(sensitivity, "results/estimates/sensitivity/humidity_adjustment.csv")
}

# Cumulative lag windows 0--7, 0--14 and 0--21.
compute_window_sensitivity <- function() {

  suppressPackageStartupMessages(library(dlnm))

  maximum_lag <- 21L
  window_ends <- c(7L, 14L, 21L)

  d <- as.data.table(read_rds(file.path(OUTPUT_ROOT, "data/derived/case_crossover.rds")))

  setorder(d, case_id, ref_date)

  retained <- read_rds(file.path(OUTPUT_ROOT, "results/models/national_association.rds"))

  temperature_history <- as.matrix(d[, paste0("tmean_lag", 0:maximum_lag), with = FALSE])

  basis <- crossbasis(temperature_history, lag = c(0L, maximum_lag), argvar = list(
    fun = "ns", knots = retained$exposure_knots,
    Boundary.knots = retained$exposure_bounds
  ), arglag = list(
    fun = "ns", knots = retained$lag_knots,
    Boundary.knots = c(0L, maximum_lag)
  ))

  rm(temperature_history)

  invisible(gc())

  case_temperature <- d[case == 1L, tmean_lag0]
  case_quantiles <- quantile(case_temperature,
    probs = c(0.01, 0.5, 0.99), na.rm = TRUE, names = FALSE,
    type = 7
  )

  names(case_quantiles) <- c("P1", "P50", "P99")

  curve_grid <- seq(round(case_quantiles[["P1"]], 2), round(case_quantiles[["P99"]], 2), by = 0.02)

  curve_grid <- sort(unique(c(curve_grid, round(case_quantiles, 2))))
  calculate_window <- function(window_end) {
    preliminary <- crosspred(basis,
      coef = retained$coefficients, vcov = retained$vcov, at = curve_grid,
      cen = round(case_quantiles[["P50"]], 2), lag = c(0L, window_end), model.link = "log"
    )
    mmt <- preliminary$predvar[which.min(preliminary$allRRfit)]
    mmt_percentile <- mean(case_temperature <= mmt, na.rm = TRUE) * 100
    prediction_temperatures <- sort(unique(c(curve_grid, round(case_quantiles[c("P1", "P99")], 2), mmt)))
    pred <- crosspred(basis,
      coef = retained$coefficients, vcov = retained$vcov, at = prediction_temperatures,
      cen = mmt, lag = c(0L, window_end), model.link = "log"
    )
    window_label <- sprintf("Lag 0-%d", window_end)
    curve <- data.table(
      WindowEnd = window_end, Window = window_label, TemperatureC = pred$predvar, OR = as.numeric(pred$allRRfit),
      Lower95 = as.numeric(pred$allRRlow), Upper95 = as.numeric(pred$allRRhigh), MMT_C = mmt, MMT_Percentile = mmt_percentile
    )
    contrast_values <- round(case_quantiles[c("P1", "P99")], 2)
    summary <- rbindlist(lapply(names(contrast_values), function(label) {
      temperature <- unname(contrast_values[[label]])
      row <- curve[which.min(abs(TemperatureC - temperature))]
      data.table(
        WindowEnd = window_end, Window = window_label, Exposure = if (label == "P1") {
          "P1 cold"
        } else {
          "P99 heat"
        }, TemperatureC = temperature, MMT_C = mmt, MMT_Percentile = mmt_percentile, OR = row$OR,
        Lower95 = row$Lower95, Upper95 = row$Upper95, PercentChange = (row$OR - 1) * 100, PercentChangeLower95 = (row$Lower95 -
          1) * 100, PercentChangeUpper95 = (row$Upper95 - 1) * 100, Cases = uniqueN(d$case_id)
      )
    }))
    summary
  }

  summaries <- rbindlist(lapply(window_ends, calculate_window))
  window_order <- sprintf("Lag 0-%d", window_ends)

  summaries[, Window := factor(Window, levels = window_order)]

  write_csv(summaries, "results/estimates/sensitivity/cumulative_windows.csv")
}

# State two-stage model: at least 100 cases, unstructured REML.
fit_state_first_stage <- function() {
  suppressPackageStartupMessages({library(dlnm);library(mvmeta);library(splines);library(survival)})
  captured<-function(expr){warnings<-character();value<-withCallingHandlers(
    tryCatch(expr,error=function(e)e),warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')})
    list(value=value,warnings=unique(warnings))}
  state_abbrev<-c('11'='RO','12'='AC','13'='AM','14'='RR','15'='PA','16'='AP','17'='TO',
   '21'='MA','22'='PI','23'='CE','24'='RN','25'='PB','26'='PE','27'='AL','28'='SE','29'='BA',
   '31'='MG','32'='ES','33'='RJ','35'='SP','41'='PR','42'='SC','43'='RS',
   '50'='MS','51'='MT','52'='GO','53'='DF')
  d<-as.data.table(read_rds(file.path(OUTPUT_ROOT, "data/derived/case_crossover.rds")))
  setorder(d,case_id,ref_date)
  primary<-read_rds(file.path(OUTPUT_ROOT, "results/models/national_association.rds"))
  stopifnot(uniqueN(d$case_id)==52565L,nrow(d)==231152L,primary$maximum_lag==21L)
  states<-sort(unique(d$state_code));stopifnot(length(states)==27L,all(states%in%names(state_abbrev)))
  cols<-paste0('tmean_lag',0:21)
  make_basis<-function(h)crossbasis(h,lag=c(0,21),
    argvar=list(fun='ns',knots=primary$exposure_knots,Boundary.knots=primary$exposure_bounds),
    arglag=list(fun='ns',knots=primary$lag_knots,Boundary.knots=c(0,21)))
  targets<-c(10.74,29.10);reference<-23.14
  all_diagnostics <- retained <- reductions <- list()
  for(s in states){
    z<-d[state_code==s];abbr<-unname(state_abbrev[s])
    h<-as.matrix(z[,..cols]);b<-make_basis(h);ncb<-ncol(b)
    groups<-match(z$case_id,unique(z$case_id));ct<-tabulate(groups)
    within_basis<-unclass(b)-rowsum(unclass(b),groups,reorder=FALSE)[groups,,drop=FALSE]/ct[groups]
    rank_within<-qr(within_basis,tol=1e-9)$rank
    cn<-paste0('cb',seq_len(ncb));fit_data<-copy(z[,.(case,case_id,rh_lag03,holiday)])
    fit_data[,(cn):=as.data.table(unclass(b))]
    informative_holiday_strata<-z[,.(vary=diff(range(holiday))>0),by=case_id][,sum(vary)]
    holiday_term<-if(informative_holiday_strata>0)'holiday' else character()
    form<-as.formula(paste('case ~',paste(c(cn,'ns(rh_lag03,df=3)',holiday_term,'strata(case_id)'),collapse=' + ')))
    cap<-captured(clogit(form,data=fit_data,method='breslow'));f<-cap$value
    info<-data.table(StateCode=s,State=abbr,Region=z$region[1],Cases=sum(z$case),Rows=nrow(z),
      Municipalities=uniqueN(z$Code),BasisColumns=ncb,WithinStratumBasisRank=rank_within,
      MinHistoryC=min(h),MaxHistoryC=max(h),ColdTargetWithinHistoryRange=min(h)<=targets[1]&&max(h)>=targets[1],
      HeatTargetWithinHistoryRange=min(h)<=targets[2]&&max(h)>=targets[2],
      HolidayInformativeStrata=informative_holiday_strata,HolidayTermIncluded=length(holiday_term)>0,
      NuisanceNote=if(length(holiday_term))'' else 'Holiday constant within every stratum; redundant term omitted',
      Converged=FALSE,FiniteBasis=FALSE,FiniteAllCoefficients=FALSE,MinCovarianceEigenvalue=NA_real_,
      CovarianceConditionNumber=NA_real_,IncludedInMeta=FALSE,Reason='',Warning=paste(cap$warnings,collapse=' | '))
    if(inherits(f,'error')){info[,Reason:=conditionMessage(f)];all_diagnostics[[s]]<-info;next}
    beta<-coef(f)[cn];V<-vcov(f)[cn,cn,drop=FALSE]
    finite_basis<-all(is.finite(beta))&&all(is.finite(V))
    finite_all<-all(is.finite(coef(f)))
    eig<-if(finite_basis)eigen(V,symmetric=TRUE,only.values=TRUE)$values else NA_real_
    info[,`:=`(Converged=is.null(f$fail),FiniteBasis=finite_basis,FiniteAllCoefficients=finite_all,
      MinCovarianceEigenvalue=if(finite_basis)min(eig) else NA_real_,
      CovarianceConditionNumber=if(finite_basis&&min(eig)>0)max(eig)/min(eig) else NA_real_)]
    reason<-character()
    if(!is.null(f$fail))reason<-c(reason,paste('Fit failure:',f$fail))
    if(rank_within<ncb)reason<-c(reason,'Common cross-basis not fully identifiable within matched strata')
    if(!finite_basis)reason<-c(reason,'Non-finite temperature coefficients or covariance')
    if(!finite_all)reason<-c(reason,'Non-finite fitted coefficients, including nuisance terms')
    if(finite_basis&&min(eig)<=0)reason<-c(reason,'Temperature covariance not positive definite')
    if(any(grepl('infinite|converg|iteration|singular',cap$warnings,ignore.case=TRUE)))
      reason<-c(reason,'Convergence/infinite-coefficient warning')
    retained[[s]]<-list(state=abbr,coefficients=beta,vcov=V,formula=deparse(form),
      exposure_knots=primary$exposure_knots,exposure_bounds=primary$exposure_bounds,lag_knots=primary$lag_knots,
      cases=sum(z$case),rows=nrow(z),warnings=cap$warnings)
    if(!length(reason)){
      cap_red<-captured(crossreduce(b,coef=beta,vcov=V,type='overall',lag=c(0,5),cen=reference,model.link='logit'))
      red<-cap_red$value
      if(inherits(red,'error'))reason<-c(reason,paste('Reduction failed:',conditionMessage(red))) else {
        v_red<-vcov(red);b_red<-coef(red)
        if(any(!is.finite(b_red))||any(!is.finite(v_red))||min(eigen(v_red,symmetric=TRUE,only.values=TRUE)$values)<=0)
          reason<-c(reason,'Reduced coefficients/covariance failed estimability check')
      }
    }
    if(length(reason)){info[,Reason:=paste(unique(reason),collapse=' | ')]}else{
      info[,IncludedInMeta:=TRUE];reductions[[s]]<-list(state=abbr,coefficients=b_red,vcov=v_red,cases=sum(z$case),rows=nrow(z))
    }
    all_diagnostics[[s]]<-info
  }
  diagnostics<-rbindlist(all_diagnostics,fill=TRUE)
  save_rds(list(models=retained,reduced=reductions,diagnostics=diagnostics),'results/models/state_first_stage.rds')
}

fit_state_sensitivity <- function() {

  stage_path <- file.path(OUTPUT_ROOT, "results/models/state_first_stage.rds")

  suppressPackageStartupMessages({
    library(dlnm)
    library(mvmeta)
    library(survival)
    library(splines)
  })

  captured <- function(expr) {
    w <- character()
    v <- withCallingHandlers(tryCatch(expr, error = function(e) e), warning = function(e) {
      w <<- c(w, conditionMessage(e))
      invokeRestart("muffleWarning")
    })
    list(value = v, warnings = unique(w))
  }

  data_path <- file.path(OUTPUT_ROOT, "data/derived/case_crossover.rds")

  primary_path <- file.path(OUTPUT_ROOT, "results/models/national_association.rds")

  x <- read_rds(stage_path)
  primary <- read_rds(primary_path)
  d <- as.data.table(read_rds(data_path))

  setorder(d, case_id, ref_date)

  cols <- paste0("tmean_lag", 0:21)
  targets <- c(10.74, 29.1)
  reference <- 23.14
  make_basis <- function(h) {
    crossbasis(h, lag = c(0, 21), argvar = list(
      fun = "ns", knots = primary$exposure_knots,
      Boundary.knots = primary$exposure_bounds
    ), arglag = list(fun = "ns", knots = primary$lag_knots, Boundary.knots = c(
      0,
      21
    )))
  }

  grid <- sort(unique(c(seq(targets[1], targets[2], by = 0.02), targets, reference)))

  one <- onebasis(grid, fun = "ns", knots = primary$exposure_knots, Boundary.knots = primary$exposure_bounds)

  summaries <- list()
  diagnostics <- list()
  predict_red <- function(name, beta, V, k, cases, rows) {
    p <- crosspred(one, coef = beta, vcov = V, at = grid, cen = reference, model.link = "logit")
    i <- match(targets, p$predvar)
    list(summary = data.table(
      Threshold = k, Model = name, Exposure = c("Cold", "Heat"), TemperatureC = targets,
      ReferenceC = reference, OwnMMT = p$predvar[which.min(p$allRRfit)], OR = as.numeric(p$allRRfit[i]),
      Lower95 = as.numeric(p$allRRlow[i]), Upper95 = as.numeric(p$allRRhigh[i]), Cases = cases, Rows = rows,
      MaximumLag = 21L, CumulativeWindow = "0-5"
    ), curve = data.table(
      Threshold = k, Model = name,
      TemperatureC = p$predvar, OR = as.numeric(p$allRRfit), Lower95 = as.numeric(p$allRRlow), Upper95 = as.numeric(p$allRRhigh)
    ))
  }

  for (k in 100L) {
    codes <- as.character(x$diagnostics[IncludedInMeta == TRUE & Cases >= k, StateCode])
    r <- x$reduced[codes]
    stopifnot(length(r) == length(codes), !any(vapply(r, is.null, logical(1))))
    z <- d[state_code %in% codes]
    cases <- sum(z$case)
    rows <- nrow(z)
    stopifnot(cases == sum(vapply(r, `[[`, numeric(1), "cases")), rows == sum(vapply(
      r, `[[`, integer(1),
      "rows"
    )), uniqueN(z$case_id) == cases, z[, all(sum(case) == 1 & .N >= 2), by = case_id][, all(V1)])
    y <- do.call(rbind, lapply(r, `[[`, "coefficients"))
    rownames(y) <- vapply(r, `[[`, character(1), "state")
    S <- lapply(r, `[[`, "vcov")
    metas <- list(REML_unstructured = captured(mvmeta(y, S = S, method = "reml", bscov = "unstr", control = list(
      maxiter = 1000,
      hessian = TRUE
    ))))
    for (n in names(metas)) {
      cap <- metas[[n]]
      f <- cap$value
      if (inherits(f, "error")) {
        stop(n, ": ", conditionMessage(f))
      }
      stopifnot(all(is.finite(coef(f))), all(is.finite(vcov(f))), is.null(f$converged) || isTRUE(f$converged))
      q <- qtest(f)
      e <- if (!is.null(f$Psi)) {
        eigen(f$Psi, symmetric = TRUE, only.values = TRUE)$values
      } else {
        numeric()
      }
      diagnostics[[paste(k, n)]] <- data.table(
        Threshold = k, Model = n, Converged = if (is.null(f$converged)) {
          NA
        } else {
          f$converged
        }, States = length(codes), Q = unname(q$Q[1]), Qdf = unname(q$df[1]), Qp = unname(q$pvalue[1]),
        I2 = max(0, (q$Q[1] - q$df[1]) / q$Q[1]) * 100, BetweenCovarianceRank = if (length(e)) {
          sum(e > max(e) * sqrt(.Machine$double.eps))
        } else {
          NA_integer_
        }, LargestEigenvalue = if (length(e)) {
          max(e)
        } else {
          NA_real_
        }, SmallestEigenvalue = if (length(e)) {
          min(e)
        } else {
          NA_real_
        }, Warning = paste(cap$warnings, collapse = " | ")
      )
      p <- predict_red(n, coef(f), vcov(f), k, cases, rows)
      summaries[[paste(k, n)]] <- p$summary
    }
    save_rds(list(
      meta = metas, included_states = codes, threshold = k, cases = cases, rows = rows, exposure_knots = primary$exposure_knots,
      exposure_bounds = primary$exposure_bounds, lag_knots = primary$lag_knots, maximum_lag = 21L,
      cumulative_window = c(0, 5), reference = reference
    ), "results/models/state_two_stage.rds")
  }

  ss <- rbindlist(summaries)

  write_csv(ss, "results/estimates/sensitivity/state_two_stage.csv")

  write_csv(rbindlist(diagnostics), "results/estimates/sensitivity/state_heterogeneity.csv")
}

# Municipality-specific temperature percentile model.
fit_percentile_sensitivity <- function() {

  suppressPackageStartupMessages({
    library(dlnm)
    library(mvmeta)
    library(splines)
    library(survival)
  })

  captured <- function(expr) {
    warnings <- character()
    value <- withCallingHandlers(tryCatch(expr, error = function(e) e), warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
    list(value = value, warnings = unique(warnings))
  }

  d <- as.data.table(read_rds(file.path(OUTPUT_ROOT, "data/derived/case_crossover.rds")))

  setorder(d, case_id, ref_date)

  stopifnot(!anyNA(d$rh_lag03), !anyNA(d$holiday))

  primary <- read_rds(file.path(OUTPUT_ROOT, "results/models/national_association.rds"))

  history_cols <- paste0("tmean_lag", 0:21)
  history <- as.matrix(d[, ..history_cols])

  stopifnot(all(is.finite(history)))

  make_basis <- function(h, m) {
    crossbasis(h, lag = c(0, 21), argvar = list(
      fun = "ns", knots = m$exposure_knots,
      Boundary.knots = m$exposure_bounds
    ), arglag = list(fun = "ns", knots = m$lag_knots, Boundary.knots = c(
      0,
      21
    )))
  }

  curves <- list()
  contrasts <- list()

  predict_reduced <- function(name, beta, V, scale, grid, one, at, fixed_ref = NULL) {
    provisional <- crosspred(one, coef = beta, vcov = V, at = grid, cen = median(grid), model.link = "logit")
    own_mmt <- provisional$predvar[which.min(provisional$allRRfit)]
    reference <- if (is.null(fixed_ref)) {
      own_mmt
    } else {
      fixed_ref
    }
    p <- crosspred(one,
      coef = beta, vcov = V, at = sort(unique(c(grid, at, reference))), cen = reference,
      model.link = "logit"
    )
    curve <- data.table(
      Model = name, Scale = scale, ExposureValue = p$predvar, OR = as.numeric(p$allRRfit),
      Lower95 = as.numeric(p$allRRlow), Upper95 = as.numeric(p$allRRhigh), Reference = reference, OwnMMT = own_mmt,
      ReferenceRule = if (is.null(fixed_ref)) {
        "Own MMT"
      } else {
        "Fixed primary MMT"
      }
    )
    ix <- match(at, p$predvar)
    stopifnot(!anyNA(ix))
    contrast <- data.table(
      Model = name, Scale = scale, Exposure = c("Cold", "Heat"), ExposureValue = at,
      Reference = reference, OwnMMT = own_mmt, ReferenceRule = if (is.null(fixed_ref)) {
        "Own MMT"
      } else {
        "Fixed primary MMT"
      }, OR = as.numeric(p$allRRfit[ix]), Lower95 = as.numeric(p$allRRlow[ix]),
      Upper95 = as.numeric(p$allRRhigh[ix]), Cases = 52565L, Rows = 231152L, MaximumLag = 21L, CumulativeWindow = "0-5",
      MMTAtBoundary = own_mmt == min(grid) || own_mmt == max(grid)
    )
    list(curve = curve, contrast = contrast)
  }

  weather_env <- new.env()

  loaded <- load(resolve_input(file.path(RAW_ROOT, "weather/municipal_daily.RData")), envir = weather_env)

  stopifnot("data_meteorology_BRA_list" %in% loaded)

  weather <- weather_env$data_meteorology_BRA_list

  percentile_history <- matrix(NA_real_, nrow(history), ncol(history), dimnames = list(NULL, paste0(
    "percentile_lag",
    0:21
  )))

  row_groups <- split(seq_len(nrow(d)), d$Code)
  alignment_error <- 0

  for (j in seq_along(row_groups)) {
    code <- names(row_groups)[j]
    ix <- row_groups[[j]]
    w <- weather[[code]]
    if (is.null(w)) {
      stop(paste("Missing municipality weather", code))
    }
    dates <- as.Date(w$date)
    take <- dates >= as.Date("2000-01-01") & dates <= as.Date("2019-12-31") & is.finite(w$Tmean)
    climate <- w$Tmean[take]
    stopifnot(length(climate) > 100)
    F <- ecdf(climate)
    percentile_history[ix, ] <- matrix(100 * F(as.vector(history[ix, , drop = FALSE])), length(ix), 22)
    observed <- w$Tmean[match(as.Date(d$ref_date[ix]), dates)]
    stopifnot(all(is.finite(observed)))
    alignment_error <- max(alignment_error, max(abs(observed - history[ix, 1])))
  }

  stopifnot(all(is.finite(percentile_history)), min(percentile_history) >= 0, max(percentile_history) <=
    100, alignment_error < 1e-05)

  rm(weather, weather_env, history, row_groups)

  invisible(gc())

  message("Fitting municipality-percentile national conditional logistic DLNM")

  basis_p <- crossbasis(percentile_history, lag = c(0, 21), argvar = list(fun = "ns", knots = c(
    10, 75,
    90
  ), Boundary.knots = c(0, 100)), arglag = list(fun = "ns", knots = primary$lag_knots, Boundary.knots = c(
    0,
    21
  )))

  cn <- paste0("pct_cb", seq_len(ncol(basis_p)))
  fit_data <- copy(d[, .(case, case_id, rh_lag03, holiday)])

  fit_data[, `:=`((cn), as.data.table(unclass(basis_p)))]

  formula_p <- as.formula(paste("case ~", paste(c(cn, "ns(rh_lag03,df=3)", "holiday", "strata(case_id)"),
    collapse = " + "
  )))

  fit_result <- captured(clogit(formula_p, data = fit_data, method = "breslow"))
  fit <- fit_result$value

  if (inherits(fit, "error")) stop(conditionMessage(fit))

  beta <- coef(fit)[cn]
  V <- vcov(fit)[cn, cn, drop = FALSE]

  stopifnot(is.null(fit$fail), all(is.finite(beta)), all(is.finite(V)), fit$n == nrow(d), fit$nevent ==
    52565L)

  red_p <- crossreduce(basis_p, coef = beta, vcov = V, type = "overall", lag = c(0, 5), cen = 50, model.link = "logit")

  grid_p <- seq(1, 99, by = 0.1)

  one_p <- onebasis(grid_p, fun = "ns", knots = c(10, 75, 90), Boundary.knots = c(0, 100))

  z <- predict_reduced(
    "Municipal percentile national", coef(red_p), vcov(red_p), "Municipal climate percentile",
    grid_p, one_p, c(1, 99)
  )

  curves[[length(curves) + 1L]] <- z$curve

  contrasts[[length(contrasts) + 1L]] <- z$contrast

  write_csv(rbindlist(contrasts), "results/estimates/sensitivity/municipal_percentiles.csv")

  save_rds(list(
    coefficients = beta, vcov = V, reduced_coefficients = coef(red_p), reduced_vcov = vcov(red_p),
    formula = deparse(formula_p), exposure_knots = c(10, 75, 90), exposure_bounds = c(0, 100), lag_knots = primary$lag_knots,
    maximum_lag = 21L, cumulative_lag = 5L, reference = z$contrast$Reference[1], mapping = "100 * empirical CDF within municipality, all finite 2000-2019 Tmean days",
    cases = fit$nevent, rows = fit$n, warnings = fit_result$warnings
  ), "results/models/municipal_percentiles.rds")
}

# Format the adopted results for supplementary Table S4.
make_table_s4 <- function() {
  num <- file.path(OUTPUT_ROOT, "results/estimates")
  fmt <- function(or, lo, hi) sprintf("%.3f (%.3f, %.3f)", or, lo, hi)
  primary <- read_csv(file.path(num, "main/national_association.csv"))

  rows <- list(data.table(`Model specification` = "", `Cold: 1st percentile vs MMT` = "OR (95% CI)", `Heat: 99th percentile vs MMT` = "OR (95% CI)"))

  add <- function(label, or, lo, hi) {
    rows[[length(rows) + 1L]] <<- data.table(
      `Model specification` = label,
      `Cold: 1st percentile vs MMT` = fmt(or[1], lo[1], hi[1]), `Heat: 99th percentile vs MMT` = fmt(
        or[2],
        lo[2], hi[2]
      )
    )
  }

  add("Primary model", primary$OR, primary$Lower95, primary$Upper95)

  basis <- read_csv(file.path(num, "sensitivity/exposure_and_splines.csv"))
  occurrence <- read_csv(file.path(num, "sensitivity/place_of_death.csv"))
  occurrence <- occurrence[match(c("P1 cold", "P99 heat"), Exposure)]
  add("Exposure assigned according to municipality of death occurrence for all deaths", occurrence$OR, occurrence$Lower95, occurrence$Upper95)

  labels <- c(
    exposure_2_knots = "Exposure knots at P10 and P90",
    exposure_4_knots = "Exposure knots at P10, P50, P75, and P90", lag_2_knots = "Two log-spaced internal lag knots",
    lag_4_knots = "Four log-spaced internal lag knots"
  )

  for (id in names(labels)) {
    r <- basis[ModelID == id]
    r <- r[order(Exposure)]
    add(labels[[id]], r$OR, r$Lower95, r$Upper95)
  }

  windows <- read_csv(file.path(num, "sensitivity/cumulative_windows.csv"))

  for (k in c(7L, 14L, 21L)) {
    r <- windows[WindowEnd == k][order(Exposure)]
    add(paste0("Cumulative lag 0–", k, " days"), r$OR, r$Lower95, r$Upper95)
  }

  rh <- read_csv(file.path(num, "sensitivity/humidity_adjustment.csv"))

  stopifnot(nrow(rh) == 1L)

  add(
    "Relative-humidity natural spline with 4 df", c(rh$ColdOR, rh$HeatOR), c(rh$ColdLower95, rh$HeatLower95),
    c(rh$ColdUpper95, rh$HeatUpper95)
  )

  # The current supplement contains the primary model and nine sensitivity rows.
  write_csv(rbindlist(rows), "results/estimates/supplement/TableS4.csv")

  for (entry in list(c("state_two_stage.csv", "State-level two-stage model"), c(
    "municipal_percentiles.csv",
    "Municipality-specific temperature percentile model"
  ))) {
    r <- read_csv(file.path(num, "sensitivity", entry[1]))
    r <- r[order(Exposure)]
    add(entry[2], r$OR, r$Lower95, r$Upper95)
  }

  write_csv(rbindlist(rows), "results/estimates/sensitivity/all_models.csv")
}

# Prepare municipal mortality and the monthly series using the adopted IBGE2024 denominator.
prepare_figure1_data <- function() {
  pop <- read_csv(file.path(OUTPUT_ROOT, "data/derived/municipal_population.csv"), colClasses = c(code6 = "character"))

  d <- analysis_deaths()

  d[, `:=`(Code = as.character(Code), year = as.integer(format(death_date, "%Y")), month = as.integer(format(
    death_date,
    "%m"
  )))]

  stopifnot(nrow(pop) == 5570L * 20L, nrow(d) == 52565L, !anyDuplicated(pop[, .(code6, year)]))

  map <- read_csv(file.path(OUTPUT_ROOT,"data/derived/municipal_exposure_municipal_means.csv"), colClasses = c(Code = "character"))[
    ,
    .(Code, MeanTemperatureC, region)
  ]

  den <- pop[, .(PopulationPersonYears = sum(population), PopulationYears = .N), by = .(Code = code6)]

  counts <- d[, .(K550Deaths = .N), by = Code]

  map <- merge(merge(map, den, by = "Code", all = TRUE), counts, by = "Code", all.x = TRUE)

  map[is.na(K550Deaths), K550Deaths := 0L]

  map[, K550CrudeRatePer100kPY := 1e+05 * K550Deaths / PopulationPersonYears]

  setcolorder(map, c(
    "Code", "K550Deaths", "PopulationPersonYears", "PopulationYears", "K550CrudeRatePer100kPY",
    "MeanTemperatureC", "region"
  ))

  setorder(map, Code)

  stopifnot(nrow(map) == 5570L, sum(map$K550Deaths) == 52565L, all(map$PopulationYears == 20L))

  write_csv(map, "results/estimates/figures/Figure1_municipal.csv")

  monthly <- read_csv(file.path(OUTPUT_ROOT,"data/derived/municipal_exposure_monthly_temperature.csv"))

  monthly[, date := as.Date(date)]
  counts_monthly <- d[, .(deaths = .N), by = .(year, month)]
  monthly[, deaths := 0L]
  monthly[counts_monthly, deaths := i.deaths, on = .(year, month)]

  annual <- pop[, .(population = sum(population)), by = year]
  monthly <- merge(monthly, annual, by = "year", all.x = TRUE, sort = FALSE)

  setorder(monthly, date)

  days <- as.integer(as.Date(sprintf("%04d-%02d-01", monthly$year + as.integer(monthly$month == 12), monthly$month %% 12 +
    1)) - monthly$date)

  monthly[, annualized_monthly_rate_per_100k := 1e+05 * deaths / (as.numeric(population) * days / 365.25)]
  monthly[, displayed_rate_per_100k := annualized_monthly_rate_per_100k]
  # Preserve the adopted October 27 display: incomplete initial lag-history
  # month retained in numeric data, omitted from the displayed mortality line.
  monthly[date == as.Date("2000-01-01"), displayed_rate_per_100k := NA_real_]

  monthly[, population := NULL]

  setcolorder(monthly, c(
    "year", "month", "mean_temperature_c", "population_coverage", "represented_municipalities",
    "deaths", "annualized_monthly_rate_per_100k", "date"
  ))

  stopifnot(nrow(monthly) == 240L, sum(monthly$deaths) == 52565L)

  write_csv(monthly, "results/estimates/figures/Figure1_monthly.csv")
  list(municipal = map, monthly = monthly)
}

# Prepare cumulative estimates, lag-specific estimates and case-day temperature frequencies.
prepare_figure2_data <- function() {
  suppressPackageStartupMessages({
    library(dlnm)
  })

  maximum_lag <- primary_max_lag
  report_lag <- primary_cumulative_lag

  d <- as.data.table(read_rds(file.path(OUTPUT_ROOT, "data/derived/case_crossover.rds")))

  setorder(d, case_id, ref_date)

  retained <- read_rds(file.path(OUTPUT_ROOT, "results/models/national_association.rds"))
  history <- as.matrix(d[, paste0("tmean_lag", 0:maximum_lag), with = FALSE])
  basis <- crossbasis(history, lag = c(0L, maximum_lag), argvar = list(
    fun = "ns", knots = retained$exposure_knots,
    Boundary.knots = retained$exposure_bounds
  ), arglag = list(
    fun = "ns", knots = retained$lag_knots,
    Boundary.knots = c(0L, maximum_lag)
  ))

  rm(history)

  invisible(gc())

  case_temperature <- d[case == 1L, tmean_lag0]
  case_quantiles <- quantile(case_temperature,
    probs = c(0.01, 0.5, 0.99), na.rm = TRUE, names = FALSE,
    type = 7
  )

  names(case_quantiles) <- c("P1", "P50", "P99")

  curve_grid <- seq(round(case_quantiles[["P1"]], 2), round(case_quantiles[["P99"]], 2), by = 0.01)

  curve_grid <- sort(unique(c(curve_grid, round(case_quantiles, 2))))
  provisional <- crosspred(basis,
    coef = retained$coefficients, vcov = retained$vcov, at = curve_grid,
    cen = round(case_quantiles[["P50"]], 2), lag = c(0L, report_lag), model.link = "log"
  )

  mmt <- provisional$predvar[which.min(provisional$allRRfit)]
  mmt_percentile <- mean(case_temperature <= mmt, na.rm = TRUE) * 100

  boundary_mmt <- abs(mmt - case_quantiles[["P1"]]) <= 0.02 || abs(mmt - case_quantiles[["P99"]]) <= 0.02

  prediction_temperatures <- sort(unique(c(curve_grid, round(case_quantiles[c("P1", "P99")], 2), mmt)))

  pred <- crosspred(basis,
    coef = retained$coefficients, vcov = retained$vcov, at = prediction_temperatures,
    cen = mmt, lag = c(0L, report_lag), model.link = "log"
  )

  curve <- data.table(
    TemperatureC = pred$predvar, OR = as.numeric(pred$allRRfit), Lower95 = as.numeric(pred$allRRlow),
    Upper95 = as.numeric(pred$allRRhigh), ReferenceC = mmt, ReferencePercentile = mmt_percentile, CumulativeLag = sprintf(
      "0-%d",
      report_lag
    )
  )

  comparison_values <- round(case_quantiles[c("P1", "P99")], 2)
  summary_rows <- rbindlist(lapply(names(comparison_values), function(label) {
    temperature <- unname(comparison_values[[label]])
    row <- curve[which.min(abs(TemperatureC - temperature))]
    data.table(
      Exposure = if (label == "P1") {
        "P1 cold"
      } else {
        "P99 heat"
      }, TemperatureC = temperature, Reference = "Lag 0-5 cumulative MMT", ReferenceC = mmt,
      ReferencePercentile = mmt_percentile, MMTAtBoundary = boundary_mmt, CumulativeLag = "0-5", OR = row$OR,
      Lower95 = row$Lower95, Upper95 = row$Upper95, PercentChange = 100 * (row$OR - 1), PercentChangeLower95 = 100 *
        (row$Lower95 - 1), PercentChangeUpper95 = 100 * (row$Upper95 - 1), Cases = uniqueN(d$case_id)
    )
  }))

  write_csv(summary_rows, "results/estimates/main/national_association.csv")

  lag_prediction <- crosspred(basis, coef = retained$coefficients, vcov = retained$vcov, at = c(
    comparison_values[["P1"]],
    comparison_values[["P99"]]
  ), cen = mmt, cumul = TRUE, model.link = "log")

  lag_curves <- rbindlist(lapply(seq_along(lag_prediction$predvar), function(i) {
    data.table(
      contrast = c("P1 cold", "P99 heat")[i], temperature_c = lag_prediction$predvar[i], reference_c = mmt,
      lag = 0:maximum_lag, or = as.numeric(lag_prediction$matRRfit[i, ]), lower95 = as.numeric(lag_prediction$matRRlow[i, ]), upper95 = as.numeric(lag_prediction$matRRhigh[i, ]), cumulative_or = as.numeric(lag_prediction$cumRRfit[i, ]), cumulative_lower95 = as.numeric(lag_prediction$cumRRlow[i, ]), cumulative_upper95 = as.numeric(lag_prediction$cumRRhigh[i, ])
    )
  }))

  write_csv(lag_curves, "results/estimates/figures/Figure2_lags.csv")

  h <- hist(case_temperature[between(case_temperature, min(curve$TemperatureC), max(curve$TemperatureC))],
    breaks = seq(10, 30, 0.5), plot = FALSE
  )

  histogram <- data.table(xmin = head(h$breaks, -1L), xmax = tail(h$breaks, -1L), proportion = h$counts / sum(h$counts))

  write_csv(histogram, "results/estimates/figures/Figure2_histogram.csv")

  curve[, exposure := ifelse(TemperatureC < mmt, "Cold", "Heat")]
  write_csv(curve, "results/estimates/figures/Figure2_curve.csv")
  list(curve = curve, lag = lag_curves, histogram = histogram)
}

prepare_figure3_data <- function() {
  suppressPackageStartupMessages(library(dlnm))
  primary <- read_rds(file.path(OUTPUT_ROOT, "results/models/national_association.rds"))
  individual <- read_rds(file.path(OUTPUT_ROOT, "results/models/individual_subgroups.rds"))
  contextual <- read_rds(file.path(OUTPUT_ROOT, "results/models/contextual_subgroups.rds"))
  summary <- read_csv(file.path(OUTPUT_ROOT, "results/estimates/main/national_association.csv"))
  targets <- summary$TemperatureC
  reference <- unique(summary$ReferenceC)
  history <- matrix(reference, nrow = 20L, ncol = 22L)
  basis <- crossbasis(history,
    lag = c(0L, 21L),
    argvar = list(fun = "ns", knots = primary$exposure_knots, Boundary.knots = primary$exposure_bounds),
    arglag = list(fun = "ns", knots = primary$lag_knots, Boundary.knots = c(0L, 21L))
  )

  definitions <- list(
    list(panel = "A", group = "Sex", levels = c("Male", "Female"), label = "Sex"),
    list(panel = "A", group = "Age (years)", levels = c("<65", ">=65"), label = "Age"),
    list(panel = "A", group = "Education (years)", levels = c(">=8 years", "<8 years"), label = "Education"),
    list(panel = "A", group = "Race", levels = c("White", "non-White"), label = "Race/colour"),
    list(panel = "B", group = "Air-conditioning prevalence", levels = c("Low", "High"), key = "air_conditioning", unit = "%"),
    list(panel = "B", group = "Per-capita income", levels = c("Low", "High"), key = "economic", unit = " BRL"),
    list(panel = "B", group = "Hospital-bed density", levels = c("Low", "High"), key = "hospital_beds", unit = " per 1000"),
    list(panel = "B", group = "Automobile prevalence", levels = c("Low", "High"), key = "transport", unit = "%")
  )
  rows <- list()
  for (definition in definitions) {
    for (level in definition$levels) {
      model <- if (definition$panel == "A") individual[[paste(definition$group, level, sep = ": ")]] else contextual[[definition$key]][[level]]
      prediction <- crosspred(basis, coef = model$coefficients, vcov = model$vcov, at = targets, cen = reference, lag = c(0L, 5L), model.link = "log")
      n_text <- if (model$cases < 10000L) as.character(model$cases) else format(model$cases, big.mark = " ", scientific = FALSE, trim = TRUE)
      if (definition$panel == "A") {
        display_level <- sub("^>=", "≥", level)
        row_label <- sprintf("%s: %s  (n=%s)", definition$label, display_level, n_text)
      } else {
        operator <- if (level == "Low") "≤" else ">"
        row_label <- sprintf("%s: %s (%s%.2f%s; n=%s)", definition$group, level, operator, model$cutoff, definition$unit, n_text)
      }
      rows[[length(rows) + 1L]] <- data.table(
        panel = definition$panel, subgroup = definition$group, level = level, row_label = row_label,
        exposure = c("Cold", "Heat"), temperature_c = targets, reference_c = reference, cumulative_lag = "0-5",
        log_or = as.numeric(prediction$allfit), se_log_or = as.numeric(prediction$allse), cases = model$cases,
        comparison = if (level == definition$levels[1]) "Reference" else if (definition$panel == "A") paste("vs", definition$levels[1]) else "High vs low"
      )
    }
  }
  data <- rbindlist(rows)
  data[, `:=`(or = exp(log_or), lower95 = exp(log_or - 1.96 * se_log_or), upper95 = exp(log_or + 1.96 * se_log_or), comparison_p = NA_real_, contextual_comparison_p = NA_real_, contextual_comparison_p_bh = NA_real_)]
  for (group in unique(data$subgroup)) {
    for (exposure_name in c("Cold", "Heat")) {
      pair <- data[subgroup == group & exposure == exposure_name]
      p_value <- 2 * pnorm(-abs(diff(pair$log_or) / sqrt(sum(pair$se_log_or^2))))
      if (pair$panel[1] == "A") {
        data[subgroup == group & exposure == exposure_name & comparison != "Reference", comparison_p := p_value]
      } else {
        data[subgroup == group & exposure == exposure_name, contextual_comparison_p := p_value]
      }
    }
  }
  tests <- unique(data[panel == "B", .(subgroup, exposure, contextual_comparison_p)])
  tests[, adjusted := p.adjust(contextual_comparison_p, method = "BH")]
  for (i in seq_len(nrow(tests))) data[panel == "B" & subgroup == tests$subgroup[i] & exposure == tests$exposure[i], contextual_comparison_p_bh := tests$adjusted[i]]
  data[, c("log_or", "se_log_or") := NULL]
  setcolorder(data, c("panel", "subgroup", "level", "row_label", "exposure", "temperature_c", "reference_c", "cumulative_lag", "or", "lower95", "upper95", "cases", "comparison", "comparison_p", "contextual_comparison_p", "contextual_comparison_p_bh"))
  write_csv(data, "results/estimates/figures/Figure3_estimates.csv", na = "")
  data
}

# Aggregate person-years and analytic-cohort deaths by age and sex.
prepare_figure_s1_data <- function() {
  population <- read_csv(file.path(OUTPUT_ROOT, "data/derived/age_sex_population.csv"))

  display_age_levels <- c("0-19", paste0(seq(20, 85, 5), "-", seq(24, 89, 5)), "90+")

  collapse_age <- function(age) {
    lower <- pmin(5L * floor(age / 5), 90L)
    fcase(lower < 20, "0-19", lower >= 90, "90+", default = paste0(lower, "-", lower + 4))
  }

  death_file <- file.path(OUTPUT_ROOT, "data/derived/deaths.rds")
  deaths <- analysis_deaths()

  deaths[, year := as.integer(format(death_date, "%Y"))]

  deaths <- deaths[year %between% c(2000L, 2019L)]

  deaths[, sex := fcase(SEXO == 1, "Male", SEXO == 2, "Female", default = NA_character_)]

  deaths[, age_group := fifelse(is.finite(age_years) & age_years >= 0, collapse_age(age_years), NA_character_)]

  known <- deaths[!is.na(sex) & !is.na(age_group)]
  total_analytic_deaths <- nrow(deaths)
  known_deaths <- nrow(known)

  annual <- population

  stopifnot(nrow(annual) == 20 * 2 * 16, !anyDuplicated(annual[, .(year, sex, age_group)]))

  death_summary <- known[, .(deaths = .N), by = .(sex, age_group)]

  population_summary <- annual[, .(person_years = sum(population)), by = .(sex, age_group)]

  plot_data <- merge(population_summary, death_summary, by = c("sex", "age_group"), all.x = TRUE)

  plot_data[is.na(deaths), deaths := 0L]

  stopifnot(sum(plot_data$deaths) == 52362, nrow(plot_data) == 32)

  plot_data[, `:=`(
    mortality_rate_per_100k = 1e+05 * deaths / person_years, deaths_thousands = deaths / 1000,
    sex = factor(sex, levels = c("Female", "Male")), age_group = factor(age_group, levels = display_age_levels),
    study_period = "2000-2019", analytic_k550_deaths = total_analytic_deaths, analytic_deaths_with_known_age_and_sex = known_deaths
  )]

  setorder(plot_data, age_group, sex)

  write_csv(plot_data, "results/estimates/figures/FigureS1.csv")
  plot_data
}

# Predict the temperature-by-lag surface from the retained national model.
prepare_figure_s2_data <- function() {
  suppressPackageStartupMessages(library(dlnm))

  matched <- as.data.table(read_rds(file.path(OUTPUT_ROOT, "data/derived/case_crossover.rds")))

  setorder(matched, case_id, ref_date)

  retained <- read_rds(file.path(OUTPUT_ROOT, "results/models/national_association.rds"))
  summary <- read_csv(file.path(OUTPUT_ROOT, "results/estimates/main/national_association.csv"))
  mmt <- unique(summary$ReferenceC)

  comparison_values <- c(P1 = summary[Exposure == "P1 cold", TemperatureC], P99 = summary[
    Exposure == "P99 heat",
    TemperatureC
  ])

  maximum_lag <- primary_max_lag

  history <- as.matrix(matched[, paste0("tmean_lag", 0:maximum_lag), with = FALSE])

  basis <- crossbasis(history, lag = c(0L, maximum_lag), argvar = list(
    fun = "ns", knots = retained$exposure_knots,
    Boundary.knots = retained$exposure_bounds
  ), arglag = list(
    fun = "ns", knots = retained$lag_knots,
    Boundary.knots = c(0L, maximum_lag)
  ))

  temperature_grid <- seq(comparison_values[["P1"]], comparison_values[["P99"]], by = 0.05)

  surface_prediction <- crosspred(basis,
    coef = retained$coefficients, vcov = retained$vcov, at = temperature_grid,
    cen = mmt, model.link = "log"
  )

  surface <- data.table(
    temperature_c = rep(surface_prediction$predvar, times = maximum_lag + 1L), lag = rep(0:maximum_lag,
      each = length(surface_prediction$predvar)
    ), log_or = as.numeric(surface_prediction$matfit), se = as.numeric(surface_prediction$matse),
    or = as.numeric(surface_prediction$matRRfit), lower95 = as.numeric(surface_prediction$matRRlow),
    upper95 = as.numeric(surface_prediction$matRRhigh), reference_c = mmt
  )

  write_csv(surface, "results/estimates/figures/FigureS2.csv")
  surface
}

# Merge municipal indicators and apply the locked subgroup cutoffs.
prepare_figure_s3_data <- function() {
  cutoffs <- read_csv(file.path(OUTPUT_ROOT, "results/estimates/supplement/TableS1.csv"))

  atlas <- read_csv(file.path(OUTPUT_ROOT, "data/derived/municipal_income.csv"))[year == 2000L, .(Code = sprintf(
    "%06d",
    as.integer(municipality_code6)
  ), income_per_capita = income_per_capita_brl)]

  beds <- read_csv(file.path(OUTPUT_ROOT, "data/derived/hospital_beds.csv"))[, .(Code = sprintf(
    "%06d",
    as.integer(municipality_code6)
  ), hospital_beds = hospital_beds_per_1000)]

  house <- read_csv(file.path(OUTPUT_ROOT, "data/derived/household_indicators.csv"))[, .(Code = sprintf(
    "%06d",
    as.integer(municipality_code_6)
  ), air_conditioning = households_air_conditioning_pct, automobile = households_automobile_pct)]

  d <- Reduce(function(x, y) merge(x, y, by = "Code", all = TRUE), list(atlas, beds, house))

  for (i in seq_len(nrow(cutoffs))) {
    field <- cutoffs$value_column[i]
    threshold <- cutoffs$cutoff[i]
    d[, `:=`(paste0(field, "_group"), fifelse(is.na(get(field)), NA_character_, fifelse(get(field) <=
      threshold, "Low", "High")))]
  }

  write_csv(d, "results/estimates/figures/FigureS3.csv")
  list(values = d, cutoffs = cutoffs)
}

# Apply the retained national curve to the groups used in Table2 and FigureS4.
prepare_figure_s4_data <- function() {

  suppressPackageStartupMessages({
    library(dlnm)
    library(Matrix)
  })

  setDTthreads(1L)

  m <- as.data.table(read_rds(file.path(OUTPUT_ROOT, "data/derived/case_crossover.rds")))

  setorder(m, case_id, ref_date)

  d <- m[case == 1L]

  stopifnot(nrow(d) == 52565L, uniqueN(m$case_id) == nrow(d))

  d[, age_two_group := fifelse(is.na(age_years), NA_character_, fifelse(age_years < 65, "<65", ">=65"))]

  d[, education_two_group := fcase(education == ">=8 years", ">=8 years", education %in% c(
    "1-7 years",
    "Never educated"
  ), "<8 years", default = NA_character_)]

  d[, race_two_group := fcase(race == "White", "White", race %in% c("Mixed", "Black", "Asian", "Indigenous"),
    "non-White",
    default = NA_character_
  )]

  definitions <- list(list(subgroup = "Overall", variable = NULL, levels = "Overall"), list(
    subgroup = "Sex",
    variable = "sex", levels = c("Male", "Female")
  ), list(
    subgroup = "Age (years)", variable = "age_two_group",
    levels = c("<65", ">=65")
  ), list(
    subgroup = "Education (years)", variable = "education_two_group",
    levels = c(">=8 years", "<8 years")
  ), list(
    subgroup = "Race/skin colour", variable = "race_two_group",
    levels = c("White", "non-White")
  ))

  group_rows <- list()

  for (definition in definitions) {
    for (level in definition$levels) {
      count <- if (is.null(definition$variable)) {
        nrow(d)
      } else {
        sum(d[[definition$variable]] == level, na.rm = TRUE)
      }
      group_rows[[length(group_rows) + 1L]] <- data.table(
        subgroup = definition$subgroup, level = level,
        display_order = length(group_rows) + 1L, deaths = count
      )
    }
  }

  groups <- rbindlist(group_rows)

  setorder(groups, display_order)

  ii <- jj <- integer()

  for (j in seq_len(nrow(groups))) {
    idx <- switch(groups$subgroup[j],
      Overall = seq_len(nrow(d)),
      Sex = which(d$sex == groups$level[j]),
      `Age (years)` = which(d$age_two_group == groups$level[j]),
      `Education (years)` = which(d$education_two_group ==
        groups$level[j]),
      `Race/skin colour` = which(d$race_two_group == groups$level[j])
    )
    stopifnot(length(idx) == groups$deaths[j])
    ii <- c(ii, idx)
    jj <- c(jj, rep(j, length(idx)))
  }

  member <- sparseMatrix(i = ii, j = jj, x = 1, dims = c(nrow(d), nrow(groups)))
  model_path <- file.path(OUTPUT_ROOT, "results/models/national_association.rds")
  model <- read_rds(model_path)

  set.seed(20260912L)

  draws_obj <- list(
    seed = 20260912L, coefficient_draws = MASS::mvrnorm(1000L, model$coefficients, model$vcov),
    source_model_md5 = unname(tools::md5sum(model_path))
  )

  save_rds(draws_obj, "results/models/burden_coefficient_draws.rds")

  draws <- draws_obj$coefficient_draws

  stopifnot(model$maximum_lag == 21L)

  h <- as.matrix(d[, paste0("tmean_lag", 0:21), with = FALSE])
  x <- h[, 1:6, drop = FALSE]
  reference <- 23.14
  q <- quantile(h[, 1], c(0.01, 0.99), type = 7, names = FALSE)
  basis <- crossbasis(h[1:20, ], lag = c(0L, 21L), argvar = list(
    fun = "ns", knots = model$exposure_knots,
    Boundary.knots = model$exposure_bounds
  ), arglag = list(fun = "ns", knots = model$lag_knots, Boundary.knots = c(
    0L,
    21L
  )))

  masks <- list(`Total non-optimal` = matrix(TRUE, nrow(x), 6), Cold = x <= reference, `Moderate cold` = x >=
    q[1] & x <= reference, `Extreme cold` = x < q[1], Heat = x > reference, `Moderate heat` = x > reference &
    x <= q[2], `Extreme heat` = x > q[2])

  stopifnot(all(Reduce(`+`, lapply(
    masks[c("Moderate cold", "Extreme cold", "Moderate heat", "Extreme heat")],
    function(z) z * 1L
  )) == 1L))

  results <- list()
  simulation_list <- list()

  for (s in names(masks)) {
    tx <- x
    tx[!masks[[s]]] <- reference
    n <- nrow(tx)
    raw <- dlnm:::mkXpred(type = "cb", basis = basis, at = tx, predvar = seq_len(n), predlag = 0:5, cen = reference)
    design <- matrix(0, n, ncol(raw))
    for (l in 0:5) design <- design + raw[seq_len(n) + n * l, , drop = FALSE]
    point <- as.numeric(crossprod(member, 1 - exp(-drop(design %*% model$coefficients))))
    sims <- matrix(NA_real_, nrow(groups), nrow(draws))
    for (first in seq.int(1L, nrow(draws), 50L)) {
      id <- first:min(first + 49L, nrow(draws))
      sims[, id] <- as.matrix(crossprod(member, 1 - exp(-design %*% t(draws[id, , drop = FALSE]))))
    }
    results[[s]] <- cbind(groups, data.table(scenario = s, attributable_deaths = point, attributable_deaths_lower95 = apply(sims,
      1, quantile,
      probs = 0.025, type = 7
    ), attributable_deaths_upper95 = apply(sims, 1, quantile,
      probs = 0.975, type = 7
    )))
    simulation_list[[s]] <- sims
    rm(raw, design)
    invisible(gc())
  }

  burden <- rbindlist(results)

  burden[, `:=`(
    attributable_fraction_percent = 100 * attributable_deaths / deaths, attributable_fraction_lower95 = 100 *
      attributable_deaths_lower95 / deaths, attributable_fraction_upper95 = 100 * attributable_deaths_upper95 / deaths,
    reference_c = reference, p1_c = q[1], p99_c = q[2], cumulative_lag = "0-5", common_national_curve = TRUE
  )]

  burden[, scenario_order := match(scenario, names(masks))]

  setorder(burden, display_order, scenario_order)

  {
    write_csv(burden[subgroup == "Overall"], "results/estimates/main/Table2.csv")
    write_csv(burden[subgroup != "Overall"], "results/estimates/figures/FigureS4.csv")
  }

  save_rds(list(groups = groups, simulations = simulation_list, draws_source = output_file("results/models/burden_coefficient_draws.rds")), "results/models/subgroup_burden_draws.rds")

  stopifnot(nrow(burden) == 63L, all(is.finite(burden$attributable_deaths)))

  burden[subgroup != "Overall"]
}

if (sys.nframe()==0L) {
  args <- commandArgs(trailingOnly=TRUE)
  if(length(args)!=2L) stop("Usage: Rscript code/02_models.R PRIVATE_INPUTS PRIVATE_WORK")
  script <- sub("^--file=", "", grep("^--file=",commandArgs(),value=TRUE)[1])
  here <- dirname(normalizePath(script))
  source(file.path(here,"01_data_processing.R"), encoding="UTF-8")
  RAW_ROOT <- normalizePath(args[1], mustWork = TRUE)
  OUTPUT_ROOT <- work_directory(args[2], RAW_ROOT, dirname(here))
  setDTthreads(1L)
  fit_primary_model()
  invisible(prepare_figure2_data())
  fit_individual_subgroups()
  fit_contextual_subgroups()
  make_table1()
  make_table_s2()
  make_table_s3()
  fit_occurrence_sensitivity()
  fit_exposure_sensitivity()
  fit_humidity_sensitivity()
  compute_window_sensitivity()
  fit_state_first_stage()
  fit_state_sensitivity()
  fit_percentile_sensitivity()
  make_table_s4()
  invisible(prepare_figure1_data())
  invisible(prepare_figure3_data())
  invisible(prepare_figure_s1_data())
  invisible(prepare_figure_s2_data())
  invisible(prepare_figure_s3_data())
  invisible(prepare_figure_s4_data())
  message("Aggregate estimates written to ", file.path(OUTPUT_ROOT, "results/estimates"))
}
