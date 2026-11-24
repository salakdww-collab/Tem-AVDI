# Rscript code/01_data_processing.R RAW_INPUT_DIR NEW_WORK_DIR
library(data.table)

study_start <- as.Date("2000-01-01")
study_end <- as.Date("2019-12-31")
primary_max_lag <- 21L
primary_cumulative_lag <- 5L
holiday_mmdd <- c("01-01", "04-21", "05-01", "09-07", "10-12", "11-02", "11-15", "11-20", "12-25")
state_to_region <- c(
  `11` = "North", `12` = "North", `13` = "North", `14` = "North", `15` = "North", `16` = "North", `17` = "North",
  `21` = "Northeast", `22` = "Northeast", `23` = "Northeast", `24` = "Northeast", `25` = "Northeast", `26` = "Northeast", `27` = "Northeast", `28` = "Northeast", `29` = "Northeast",
  `31` = "Southeast", `32` = "Southeast", `33` = "Southeast", `35` = "Southeast",
  `41` = "South", `42` = "South", `43` = "South",
  `50` = "CenterWest", `51` = "CenterWest", `52` = "CenterWest", `53` = "CenterWest"
)

resolve_input <- function(file) {
  actual <- normalizePath(file, mustWork=TRUE)
  allowed <- startsWith(actual, paste0(RAW_ROOT,"/")) || startsWith(actual,paste0(OUTPUT_ROOT,"/"))
  if (!allowed) stop("Input outside the study directories: ",actual)
  actual
}
read_csv <- function(file, ...) data.table::fread(resolve_input(file), ...)
read_rds <- function(file, ...) base::readRDS(resolve_input(file), ...)

output_file <- function(relative) {
  stopifnot(!grepl("(^|/)\\.\\.(/|$)", relative))
  target <- file.path(OUTPUT_ROOT, relative)
  dir.create(dirname(target), recursive = TRUE, showWarnings = FALSE)
  target
}

write_csv <- function(data, relative, ...) {
  data.table::fwrite(data, output_file(relative), ...)
}

save_rds <- function(model, relative, ...) {
  base::saveRDS(model, output_file(relative), ...)
}

# Return the fixed cohort retained by the primary case-crossover analysis.
analysis_deaths <- function() {
  deaths <- as.data.table(read_rds(file.path(OUTPUT_ROOT, "data/derived/deaths.rds")))
  matched <- as.data.table(read_rds(file.path(OUTPUT_ROOT, "data/derived/case_crossover.rds")))
  ids <- matched[case == 1L, unique(case_id)]
  stopifnot(!anyDuplicated(deaths$case_id), all(ids %in% deaths$case_id))
  deaths[case_id %in% ids]
}

# Prepare registered K55.0 deaths (2000--2019).
prepare_deaths <- function() {

  loaded <- load(resolve_input(file.path(OUTPUT_ROOT, "data/derived/SIM_raw_K550_export.RData")))

  if (!identical(loaded, "data_death_all")) {
    stop("Unexpected object(s) in death source: ", paste(loaded, collapse = ", "))
  }

  raw <- data_death_all

  rm(data_death_all)

  keep_rows <- which(raw$CAUSABAS == "K550")
  deaths <- as.data.table(raw[keep_rows, ])

  deaths[, source_row_id := keep_rows]

  setcolorder(deaths, "source_row_id")

  rm(raw, keep_rows)
  deaths[, death_date := as.IDate(DTOBITO, format = "%d%m%Y")]

  deaths[, birth_date := as.IDate(DTNASC, format = "%d%m%Y")]

  deaths[, state_code := substr(as.character(Code), 1L, 2L)]

  deaths[, region := unname(state_to_region[state_code])]

  deaths[, age_years := NA_real_]

  deaths[substr(IDADE, 1L, 1L) == "4", age_years := as.numeric(substr(IDADE, 2L, 3L))]

  deaths[substr(IDADE, 1L, 1L) == "5", age_years := 100 + as.numeric(substr(IDADE, 2L, 3L))]

  deaths[, age_group := fcase(age_years <= 49, "<=49", age_years >= 50 & age_years <= 64, "50-64", age_years >=
    65, ">=65", default = NA_character_)]

  deaths[, sex := fcase(SEXO == 1, "Male", SEXO == 2, "Female", default = NA_character_)]

  deaths[, race := fcase(RACACOR == 1, "White", RACACOR == 4, "Mixed", RACACOR == 2, "Black", RACACOR ==
    3, "Asian", RACACOR == 5, "Indigenous", default = NA_character_)]

  deaths[, education := fcase(ESC %in% c(4, 5), ">=8 years", ESC %in% c(2, 3), "1-7 years", ESC == 1,
    "Never educated",
    default = NA_character_
  )]

  setorder(deaths, death_date, Code)

  deaths[, case_id := .I]

  analysis_deaths <- deaths[!is.na(death_date) & death_date >= study_start & death_date <= study_end]

  save_rds(analysis_deaths, "data/derived/deaths.rds", compress = "gzip")
}

# Build time-stratified case-crossover records.
prepare_case_crossover <- function() {
  suppressPackageStartupMessages(library(data.table))

  cases <- as.data.table(read_rds(file.path(OUTPUT_ROOT, "data/derived/deaths.rds")))

  cases[, `:=`(
    year = as.integer(format(death_date, "%Y")), month = as.integer(format(death_date, "%m")),
    dow = as.integer(format(death_date, "%u")), case_date = death_date
  )]

  calendar <- data.table(ref_date = as.IDate(seq(study_start, study_end, by = "day")))

  calendar[, `:=`(
    year = as.integer(format(ref_date, "%Y")), month = as.integer(format(ref_date, "%m")),
    dow = as.integer(format(ref_date, "%u"))
  )]

  matched <- merge(cases, calendar, by = c("year", "month", "dow"), allow.cartesian = TRUE, sort = FALSE)

  matched[, case := as.integer(ref_date == case_date)]

  matched[, holiday := as.integer(format(ref_date, "%m-%d") %in% holiday_mmdd)]

  setorder(matched, case_id, ref_date)

  stratum_case_check <- matched[, .(ok = sum(case) == 1L), by = case_id]

  if (!all(stratum_case_check$ok)) {
    stop("At least one case-crossover stratum does not contain exactly one case day")
  }

  rm(stratum_case_check)

  loaded <- load(resolve_input(file.path(RAW_ROOT, "weather/municipal_daily.RData")))

  if (!identical(loaded, "data_meteorology_BRA_list")) {
    stop("Unexpected object(s) in weather source: ", paste(loaded, collapse = ", "))
  }

  weather <- data_meteorology_BRA_list

  rm(data_meteorology_BRA_list)

  n <- nrow(matched)

  temperature_history <- matrix(NA_real_, nrow = n, ncol = primary_max_lag + 1L, dimnames = list(
    NULL,
    paste0("tmean_lag", 0:primary_max_lag)
  ))

  rh_history <- matrix(NA_real_, nrow = n, ncol = 6L, dimnames = list(NULL, paste0("rh_lag", 0:5)))

  row_groups <- split(seq_len(n), matched$Code)

  for (code in names(row_groups)) {
    rows <- row_groups[[code]]
    city <- weather[[code]]
    if (is.null(city) || !nrow(city)) {
      next
    }
    city_dates <- as.Date(city$date)
    target_dates <- as.Date(matched$ref_date[rows])
    for (lag in 0:primary_max_lag) {
      temperature_history[rows, lag + 1L] <- city$Tmean[match(target_dates - lag, city_dates)]
    }
    for (lag in 0:5) {
      rh_history[rows, lag + 1L] <- city$RH[match(target_dates - lag, city_dates)]
    }
  }

  rm(weather, row_groups)

  invisible(gc())

  temperature_columns <- colnames(temperature_history)

  matched[, `:=`((temperature_columns), as.data.table(temperature_history))]

  matched[, rh_lag03 := rowMeans(rh_history[, 1:4, drop = FALSE], na.rm = FALSE)]

  matched[, rh_lag04 := rowMeans(rh_history[, 1:5, drop = FALSE], na.rm = FALSE)]

  matched[, rh_lag05 := rowMeans(rh_history[, 1:6, drop = FALSE], na.rm = FALSE)]

  rm(temperature_history, rh_history)

  exposure_cols <- c(paste0("tmean_lag", 0:primary_max_lag), "rh_lag03")

  matched[, complete_main_exposure := complete.cases(.SD), .SDcols = exposure_cols]

  complete_data <- matched[complete_main_exposure == TRUE]

  eligible <- complete_data[, .(has_case = sum(case) == 1L, has_control = sum(case == 0L) >= 1L), by = case_id][has_case &
    has_control, case_id]

  complete_data <- complete_data[case_id %in% eligible]

  setorder(complete_data, case_id, ref_date)

  save_rds(complete_data, "data/derived/case_crossover.rds", compress = "xz")
}

extract_sim <- function() {
  suppressPackageStartupMessages({library(data.table);library(read.dbc)})
  setDTthreads(1L)
  raw_dir <- file.path(RAW_ROOT, "SIM")
  required <- file.path(raw_dir,paste0("DOBR",2000:2019,".dbc"))
  stopifnot(all(file.exists(required)))
  parts <- vector("list",20)
  needed<-c("CODMUNRES","CODMUNOCOR","DTNASC","DTOBITO","CAUSABAS","LOCOCOR","IDADE","SEXO","RACACOR","ESTCIV","ESC","OCUP","LINHAA","LINHAB","LINHAC","LINHAD","LINHAII","CAUSABAS_O","NUMERODO")
  for(i in seq_along(required)){
   x<-as.data.table(read.dbc::read.dbc(required[i]))
   selected<-toupper(trimws(as.character(x$CAUSABAS)))=="K550"
   selected[is.na(selected)]<-FALSE
   keep<-intersect(needed,names(x))
   y<-x[selected,..keep]
   for(n in names(y))set(y,j=n,value=as.character(y[[n]]))
   y[,SOURCE_YEAR:=as.integer(1999+i)]
   parts[[i]]<-y
   rm(x,y);invisible(gc())
  }
  d<-rbindlist(parts,use.names=TRUE,fill=TRUE);rm(parts);invisible(gc())
  d[,death_date:=as.IDate(DTOBITO,format="%d%m%Y")]
  save_rds(d,"data/derived/death_certificates.rds",compress="gzip")
  # Recompute the provided analysis Code convention from raw residence/occurrence.
  d[,Code:=substr(CODMUNRES,1L,6L)]
  d[is.na(Code)|nchar(Code)!=6L|substr(Code,3L,6L)=="0000",Code:=substr(CODMUNOCOR,1L,6L)]
  export_cols<-c("CODMUNRES","CODMUNOCOR","DTNASC","DTOBITO","CAUSABAS","LOCOCOR","IDADE","SEXO","RACACOR","ESTCIV","ESC","OCUP","Code")
  data_death_all<-as.data.frame(d[,..export_cols])
  save(data_death_all,file=file.path(OUTPUT_ROOT,"data/derived/SIM_raw_K550_export.RData"),compress="gzip")
  rm(data_death_all,d);invisible(gc())
  prepare_deaths()
}

# Decode municipality population tables, census aggregates and contextual indicators.
prepare_population_context <- function() {
  derived <- file.path(OUTPUT_ROOT, "data/derived")
  save_summary <- function(x, filename) fwrite(x, file.path(derived, filename))
  tabnet_year <- function(year, sex) {
    path <- file.path(RAW_ROOT, "population", sprintf("population_%d_%s_age5.html", year, sex))
    html <- xml2::read_html(path, encoding = "ISO-8859-1", options = c("RECOVER", "NONET"))
    pre <- xml2::xml_text(xml2::xml_find_first(html, "//pre"))
    lines <- trimws(strsplit(pre, "\n", fixed = TRUE)[[1]])
    # Municipality rows start with the six-digit code; the printed total is checked separately.
    lines <- gsub('"', '', lines, fixed = TRUE)
    selected <- lines[grepl("^[0-9]{6}[[:space:]]", lines)]
    cells <- strsplit(selected, ";", fixed = TRUE)
    population <- vapply(cells, function(z) {
      a <- trimws(z[-1]); a[a == "-"] <- "0"
      a <- as.numeric(a)
      stopifnot(!anyNA(a), sum(head(a, -1)) == tail(a, 1))
      tail(a, 1)
    }, numeric(1))
    national <- strsplit(lines[grepl("^Total;", lines)], ";", fixed = TRUE)
    stopifnot(length(national) == 1L, length(cells) == 5570L,
              sum(population) == as.numeric(tail(national[[1]], 1)))
    data.table(code6 = substr(selected, 1, 6), year = year, population)
  }
  pop <- rbindlist(lapply(2000:2019, function(y) {
    rbind(tabnet_year(y, "male"), tabnet_year(y, "female"))[, .(population = sum(population)), by=.(year, code6)]
  }))
  setorder(pop, year, code6)
  stopifnot(nrow(pop) == 111400L)
  save_summary(pop[, .(code6, year, population)], "municipal_population.csv")

  age_file <- file.path(RAW_ROOT, "population", "age_sex.xlsx")
  cells <- as.data.frame(readxl::read_excel(age_file, sheet = "2) POP_GRUPO QUINQUENAL",
                     range = readxl::cell_limits(c(1,1), c(NA,25)), col_names = FALSE,
                     col_types = "text", .name_repair = "minimal"))
  names(cells) <- paste0("v", seq_len(ncol(cells)))
  years <- as.integer(unlist(cells[6,6:25]))
  stopifnot(identical(years, 2000:2019))
  cells <- cells[!is.na(cells$v5) & cells$v5 == "Brasil" & cells$v2 %in% c("Homens","Mulheres","Ambos") &
                 !is.na(cells$v1) & grepl("^[0-9]+", trimws(cells$v1)), , drop=FALSE]
  lower <- as.integer(sub("^([0-9]+).*", "\\1", trimws(cells$v1)))
  groups <- ifelse(lower < 20, "0-19", ifelse(lower >= 90, "90+", paste0(lower, "-", lower + 4)))
  age <- rbindlist(lapply(seq_along(years), function(j) {
    data.table(year=years[j], sex=cells$v2, age_group=groups, population=as.numeric(cells[[j+5]]))
  }))[, .(population=sum(population)), by=.(year,sex,age_group)]
  national <- age[sex=="Ambos", .(total=sum(population)), by=year]
  control <- merge(national, pop[, .(total_city=sum(population)), by=year], by="year")
  stopifnot(nrow(control)==20L, all(control$total==control$total_city))
  age <- age[sex!="Ambos"]; age[, sex:=fifelse(sex=="Homens","Male","Female")]
  stopifnot(nrow(age)==640L)
  save_summary(age[order(year,sex,age_group)], "age_sex_population.csv")
  weights <- rbindlist(lapply(c(2000,2002), function(y) {
    p <- fread(file.path(RAW_ROOT, "population", sprintf("POPTBR%d.csv",y)), colClasses=c(MUNIC_RES="character"))
    p[, .(municipality_code6=substr(MUNIC_RES,1,6), year=as.integer(ANO), population=as.numeric(POPULACAO))]
  }))
  stopifnot(!anyDuplicated(weights, by=c("municipality_code6","year")))
  save_summary(weights, "population_weights.csv")

  # SAS-layout byte positions in the IBGE 2000 fixed-width household file.
  # Stream one chunk at a time, retaining only municipal weighted sums.
  sum_fields <- c("households","phone","computer","automobile","air_conditioning")
  household_file <- function(path) {
    con <- file(path, open="r", encoding="latin1"); on.exit(close(con))
    chunks <- list()
    repeat {
      z <- readLines(con, n=50000L, warn=FALSE)
      if (!length(z)) break
      z <- z[nchar(z)>=167L & substr(z,72,72)=="1"]
      code <- substr(z,12,18)
      raw_weight <- trimws(substr(z,157,167))
      keep <- grepl("^[0-9]{7}$",code) & grepl("^[0-9]+$",raw_weight)
      if (!any(keep)) next
      z <- z[keep]; w <- as.numeric(raw_weight[keep])/1e8
      d <- data.table(code=code[keep], households=w, phone=w*(substr(z,109,109)=="1"),
          computer=w*(substr(z,111,111)=="1"), automobile=w*(substr(z,115,115) %in% as.character(1:9)),
          air_conditioning=w*(substr(z,117,117) %in% as.character(1:9)))
      chunks[[length(chunks)+1L]] <- d[, lapply(.SD,sum), by=code]
    }
    rbindlist(chunks)[, lapply(.SD,sum), by=code, .SDcols=sum_fields]
  }
  household_files <- list.files(file.path(RAW_ROOT, "households"), pattern="^[A-Z]{2}_2000[.]txt$", full.names=TRUE)
  stopifnot(length(household_files)==27L)
  hh <- rbindlist(lapply(household_files, household_file))[, lapply(.SD,sum), by=code, .SDcols=sum_fields]
  hh <- hh[order(code), .(municipality_code_7=code, municipality_code_6=substr(code,1,6), year=2000L,
        weighted_permanent_private_households=households, weighted_households_phone_line=phone,
        weighted_households_microcomputer=computer, weighted_households_automobile=automobile,
        weighted_households_air_conditioning=air_conditioning, households_air_conditioning_pct=100*air_conditioning/households,
        households_automobile_pct=100*automobile/households)]
  stopifnot(!anyDuplicated(hh$municipality_code_6))
  save_summary(hh, "household_indicators.csv")

  crosswalk <- Sys.getenv("AVDI_MUNICIPALITY_CROSSWALK", file.path(RAW_ROOT, "geography/municipalities.csv"))
  geo <- fread(crosswalk, colClasses = c(CD_MUN = "character"))
  normalize_name <- function(x) gsub("[^A-Z0-9]", "", toupper(stringi::stri_trans_general(x, "Latin-ASCII")))
  keys <- setNames(geo$CD_MUN,paste(normalize_name(geo$NM_MUN),geo$SIGLA,sep="|"))
  # Historical municipality spellings used by Atlas Brasil.
  aliases <- c(
    "AMPARODESAOFRANCISCO|SE"="2800100",
    "ARES|RN"="2401206",
    "AUGUSTOSEVERO|RN"="2401305",
    "ACU|RN"="2400208",
    "BARAODEMONTEALTO|MG"="3105509",
    "BRASOPOLIS|MG"="3108909",
    "DONAEUSEBIA|MG"="3122900",
    "ELDORADODOSCARAJAS|PA"="1502954",
    "EMBU|SP"="3515004",
    "FLORINIA|SP"="3516101",
    "FORTALEZADOTABOCAO|TO"="1708254",
    "GRACHOCARDOSO|SE"="2802601",
    "IGUARACI|PE"="2606903",
    "ITAPAGE|CE"="2306306",
    "MUQUEMDESAOFRANCISCO|BA"="2922250",
    "POXOREO|MT"="5107008",
    "PRESIDENTEJUSCELINO|RN"="2410306",
    "SANTAISABELDOPARA|PA"="1506500",
    "SANTATERESINHA|BA"="2928505",
    "SANTAREM|PB"="2513653",
    "SANTOANTONIODOLEVERGER|MT"="5107800",
    "SERIDO|PB"="2515401",
    "SAOLUIZ|RR"="1400605",
    "SAOLUISDOPARAITINGA|SP"="3550001",
    "SAOTHOMEDASLETRAS|MG"="3165206")
  keys[names(aliases)] <- aliases
  inc <- readxl::read_excel(file.path(RAW_ROOT,"context","income.xlsx"))
  labels <- trimws(as.character(inc[[1]]))
  keep <- grepl("\\([A-Z]{2}\\)$",labels)
  labels <- labels[keep]
  uf <- sub("^.*\\(([A-Z]{2})\\)$","\\1",labels)
  city <- trimws(sub("\\([A-Z]{2}\\)$","",labels))
  codes <- unname(keys[paste(normalize_name(city),uf,sep="|")])
  stopifnot(!anyNA(codes), !anyDuplicated(substr(codes,1,6)))
  income <- data.table(municipality_code6=substr(codes,1,6),year=2000L,
       income_per_capita_brl=as.numeric(inc[["Renda per capita 2000"]][keep]),source_label=labels)
  save_summary(income,"municipal_income.csv")
}

municipal_indicators <- function() {
  atlas <- read_csv(file.path(OUTPUT_ROOT, "data/derived/municipal_income.csv"))[year ==
    2000L, .(Code = sprintf("%06d", as.integer(municipality_code6)), income_per_capita = income_per_capita_brl)]

  households <- read_csv(file.path(OUTPUT_ROOT, "data/derived/household_indicators.csv"))[
    ,
    .(
      Code = sprintf("%06d", as.integer(municipality_code_6)), air_conditioning = households_air_conditioning_pct,
      automobile = households_automobile_pct
    )
  ]

  beds <- read_csv(file.path(OUTPUT_ROOT, "data/derived/hospital_beds.csv"))[, .(Code = sprintf(
    "%06d",
    as.integer(municipality_code6)
  ), hospital_beds = hospital_beds_per_1000)]

  population <- read_csv(file.path(OUTPUT_ROOT, "data/derived/population_weights.csv"))

  population[, Code := sprintf("%06d", as.integer(municipality_code6))]

  population_wide <- merge(population[year == 2000L, .(Code, population_2000 = population)], population[year ==
    2002L, .(Code, population_2002 = population)], by = "Code", all = TRUE)

  municipal <- Reduce(function(x, y) merge(x, y, by = "Code", all = TRUE), list(
    atlas, beds, households,
    population_wide
  ))

  if (municipal[, anyDuplicated(Code)] != 0L) {
    stop("Duplicate municipality code in contextual source data")
  }
  municipal
}

prepare_beds_and_cutoffs <- function() {
  suppressPackageStartupMessages(library(readxl))
  raw <- as.data.table(read_excel(file.path(RAW_ROOT, "context/beds.xls"),col_names=FALSE,skip=8))
  setnames(raw,paste0("V",seq_len(ncol(raw))))
  raw[, `:=`(uf_code=suppressWarnings(as.integer(V3)),municipality_component=suppressWarnings(as.integer(V4)),municipality_name_ams=as.character(V5),hospital_beds=suppressWarnings(as.numeric(gsub(" ","",as.character(V6)))))]
  published <- raw[uf_code==0L & municipality_component==0L,hospital_beds][1]
  a <- raw[!is.na(uf_code) & uf_code>0 & !is.na(municipality_component) & municipality_component>0]
  a[,municipality_code7:=sprintf("%02d%05d",uf_code,municipality_component)]
  a[,municipality_code6:=substr(municipality_code7,1,6)]
  a[is.na(hospital_beds),hospital_beds:=0]
  a<-unique(a[,.(municipality_code6,municipality_code7,municipality_name_ams,hospital_beds)],by="municipality_code6")
  pop<-read_csv(file.path(OUTPUT_ROOT,"data/derived/population_weights.csv"),colClasses=c(municipality_code6="character"))[year==2002]
  beds<-a[pop,on="municipality_code6"]
  stopifnot(!anyNA(beds$hospital_beds), sum(beds$hospital_beds)==published)
  beds[, `:=`(hospital_beds_per_1000=hospital_beds/population*1000,source_year=2002L)]
  write_csv(beds,"data/derived/hospital_beds.csv")
  weighted_median <- function(x, w) {
    ok <- is.finite(x) & is.finite(w) & w > 0
    if (!any(ok)) {
      stop("No valid observations for weighted median")
    }
    x <- x[ok]
    w <- w[ok]
    o <- order(x)
    x <- x[o]
    w <- w[o]
    x[which(cumsum(w) >= sum(w) / 2)[1L]]
  }

  municipal <- municipal_indicators()
  spec <- data.table(modifier = c("economic", "hospital_beds", "air_conditioning", "transport"), value_column = c(
    "income_per_capita",
    "hospital_beds", "air_conditioning", "automobile"
  ), label = c(
    "Per-capita income", "Hospital-bed density",
    "Household air-conditioning prevalence", "Household automobile prevalence"
  ), unit = c(
    "BRL (2010)",
    "beds per 1000 residents", "% of households", "% of households"
  ), source_year = c(
    2000L, 2002L, 2000L,
    2000L
  ), weight_column = c("population_2000", "population_2002", "population_2000", "population_2000"))

  spec[, cutoff := vapply(seq_len(.N), function(i) {
    weighted_median(municipal[[value_column[i]]], municipal[[weight_column[i]]])
  }, numeric(1))]

  spec[, cutoff_method := paste0(source_year, " population-weighted municipal median")]
  write_csv(spec, "results/estimates/supplement/TableS1.csv")
}

prepare_weather_summaries <- function() {
  e<-new.env();load(file.path(RAW_ROOT, "weather/municipal_daily.RData"),envir=e);daily_weather<-e$data_meteorology_BRA_list;rm(e)
  d<-read_rds(file.path(OUTPUT_ROOT,"data/derived/deaths.rds"))
  codes<-intersect(unique(as.character(d$Code)),names(daily_weather))
  municipal<-rbindlist(lapply(names(daily_weather),function(code){z<-daily_weather[[code]];data.table(Code=code,MeanTemperatureC=mean(z$Tmean,na.rm=TRUE),region=unname(state_to_region[substr(code,1,2)]))}))
  # Only official 5570 municipality population keys enter the map; lake geometry is separate.
  pop<-read_csv(file.path(OUTPUT_ROOT,"data/derived/municipal_population.csv"),colClasses=c(code6="character"))
  municipal<-municipal[Code %in% unique(pop$code6)]
  write_csv(municipal,"data/derived/municipal_exposure_municipal_means.csv")
  calendar<-data.table(date=seq(study_start,study_end,by="day"))
  calendar[,`:=`(year=as.integer(format(date,"%Y")),month=as.integer(format(date,"%m")))]
  citymonths<-rbindlist(lapply(codes,function(code){z<-daily_weather[[code]];ix<-match(as.Date(z$date),calendar$date);t<-data.table(year=calendar$year[ix],month=calendar$month[ix],Tmean=z$Tmean);t[,.(municipality_mean=mean(Tmean[is.finite(Tmean)]),valid_days=sum(is.finite(Tmean))),by=.(year,month)][,Code:=code]}))
  monthly<-citymonths[,.(mean_temperature_c=mean(municipality_mean[is.finite(municipality_mean)]),represented_municipalities=sum(is.finite(municipality_mean)),valid_municipal_days=sum(valid_days)),by=.(year,month)]
  monthly[,date:=as.Date(sprintf("%04d-%02d-01",year,month))]
  monthly[,population_coverage:=NA_real_] # Legacy field was stale after equal-municipality C revision; deliberately unclaimed.
  write_csv(monthly,"data/derived/municipal_exposure_monthly_temperature.csv")
  rm(daily_weather,citymonths);gc()
}

work_directory <- function(path, input, package, create = FALSE) {
  path <- if (file.exists(path)) normalizePath(path, mustWork = TRUE) else
    file.path(normalizePath(dirname(path), mustWork = TRUE), basename(path))
  roots <- normalizePath(c(input, package), mustWork = TRUE)
  if (any(startsWith(paste0(path, "/"), paste0(roots, "/"))))
    stop("Choose a work directory outside the raw inputs and the submission package.")
  if (create && !dir.create(path, mode = "0700", showWarnings = FALSE))
    stop("Data preparation requires a new work directory.")
  if (!dir.exists(path)) stop("Work directory not found: ", path)
  path
}

if (sys.nframe()==0L) {
  args <- commandArgs(trailingOnly=TRUE)
  if(length(args)!=2L) stop("Usage: Rscript code/01_data_processing.R PRIVATE_INPUTS NEW_PRIVATE_WORK")
  script <- sub("^--file=", "", grep("^--file=",commandArgs(),value=TRUE)[1])
  RAW_ROOT <- normalizePath(args[1], mustWork = TRUE)
  OUTPUT_ROOT <- work_directory(args[2], RAW_ROOT, dirname(dirname(normalizePath(script))), create = TRUE)
  extract_sim()
  prepare_population_context()
  prepare_beds_and_cutoffs()
  prepare_case_crossover()
  prepare_weather_summaries()
  message("Private analysis inputs prepared.")
}
