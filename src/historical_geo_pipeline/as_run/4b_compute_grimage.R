library(ClockBasis)
library(geograbi)
library(wateRmelon)
tidylog()
conflict_prefer("colSds", "MatrixGenerics") 


grimage1 = read.csv("/Lasagna/kying/p_2022/GEO_data/DNAmGrimAgeGitHub/input/ElasticNet_DNAmProtein_Vars_model4.csv")

do_grimage1 = \(mat, cpg = grimage1, keep = F, sex = NULL, age = NULL, same_age = FALSE, same_sex = FALSE) {
  mat = mat |>
    as.data.frame()
  rownames(mat) = as.character(mat[, 1])
  mat = mat[, -1]
  if (is.null(sex)) {
    sex = estimateSex(mat) |>
      as.data.frame() |>
      mutate(Female = ifelse(predicted_sex == "Female", 1, 0)) |>
      pull(Female)
  }
  if (is.null(age)) {
  age = as_tibble(mat, rownames = "name") |>
    ClockBasis::do_dnam_clock_human(clock = "HorvathAge", ncores = 10) |>
    pull(HorvathAge)
  }

  meta = data.frame(Age = age, Female = sex)
  if (same_age) {
    meta = meta |>
      mutate(Age = mean(Age))
  }
  if (same_sex) {
    meta = meta |>
      mutate(Female = round(mean(Female)))
  }
  meta = meta |>
    t() |>
    as.data.frame()

  clock_ls = cpg |>
    group_by(Y.pred) |>
    group_split(.keep = F) |>
    set_names(unique(sort((cpg$Y.pred))))
  alldata = do_dnam_clock_human(as_tibble(mat, rownames = "name"), clock_ls, check_scale = F, ncores = 10) |>
    cbind(as_tibble(t(meta)))

# step2: generate a raw variable of DNAmGrimAge(variable name=COX)
  alldata$COX = alldata$DNAmGDF_15 * 0.000348777412272004 + alldata$DNAmB2M * 4.59105969389204e-07 + alldata$DNAmCystatin_C * 3.49816671441537e-06 + alldata$DNAmTIMP_1 * 0.000143661105491888 + alldata$DNAmadm * 0.00790270975255529 + alldata$DNAmpai_1 * 2.55560382039825e-05 + alldata$DNAmleptin * -7.32066983502079e-06 + alldata$DNAmPACKYRS * 0.0303981613409142 + alldata$Age * 0.0300823182194075 + alldata$Female * -0.228468475622039

# step3: scaling COX to Age with units of year
  m.age = 59.63951
  sd.age = 9.049608
  m.COX = 13.20127
  sd.COX = 1.086805

  Y = (alldata$COX - m.COX) / sd.COX
  alldata$DNAmGrimAge = (Y * sd.age) + m.age
  alldata$DNAmGrimAge_missing = apply(alldata[, str_detect(base::colnames(alldata), "missing")], 1, mean, na.rm = T)
  if (! keep) {
    alldata = alldata |>
      select(sample_id, DNAmGrimAge, DNAmGrimAge_missing)
  }

  return(as_tibble(alldata))
}

F_scale = function(INPUT0, Y.pred0.name, Y.pred.name, gold) {
  out.para = subset(gold, var == "COX")
  out.para.age = subset(gold, var == "Age")
  m.age = out.para.age$mean
  sd.age = out.para.age$sd
  Y0 = INPUT0[, Y.pred0.name]
  Y = (Y0 - out.para$mean) / out.para$sd
  INPUT0[, Y.pred.name] = as.numeric((Y * sd.age) + m.age)
  return(INPUT0)
}

grimage2.rds = "/Lasagna/kying/p_2022/GEO_data/DNAmGrimAgeGitHub/input/DNAmGrimAge2_final.Rds"
grimage2 = readRDS(grimage2.rds)

do_grimage2 = \(mat, cpg = grimage2, keep = F, age = NULL, sex = NULL, same_age = FALSE, same_sex = FALSE) {
  mat = mat |>
    as.data.frame()
  rownames(mat) = as.character(mat[, 1])
  mat = mat[, -1]
  if (is.null(sex)) {
    sex = estimateSex(mat) |>
      as.data.frame() |>
      mutate(Female = ifelse(predicted_sex == "Female", 1, 0)) |>
      pull(Female)
  }
  if (is.null(age)) {
  age = as_tibble(mat, rownames = "name") |>
    ClockBasis::do_dnam_clock_human(clock = "HorvathAge", ncores = 10) |>
    pull(HorvathAge)
  }


  meta = data.frame(Age = age, Female = sex)
  if (same_age) {
    meta = meta |>
      mutate(Age = mean(Age))
  }
  if (same_sex) {
    meta = meta |>
      mutate(Female = round(mean(Female)))
  }
  meta = meta |>
    t() |>
    as.data.frame()


  grimage2 = cpg
  cpgs = grimage2[[1]]
  glmnet.final1 = grimage2[[2]]
  gold = grimage2[[3]]

  clock_ls = cpgs |>
    group_by(Y.pred) |>
    group_split(.keep = F) |>
    set_names(unique(sort((cpgs$Y.pred))))
  alldata = do_dnam_clock_human(as_tibble(mat, rownames = "name"), clock_ls, check_scale = F, ncores = 10) |>
    cbind(as_tibble(t(meta)))


# step2: generate a raw variable of DNAmGrimAge(variable name=COX)
  alldata$COX = as.numeric(as.matrix(subset(alldata, select = glmnet.final1$var)) %*% glmnet.final1$beta)
  alldata = F_scale(alldata, "COX", 'DNAmGrimAge2', gold)
  alldata$COX = NULL
  alldata$DNAmGrimAge2_missing = apply(alldata[, str_detect(base::colnames(alldata), "missing")], 1, mean, na.rm = T)
  if (! keep) {
    alldata = alldata |>
      select(sample_id, DNAmGrimAge2, DNAmGrimAge2_missing)
  }
  return(as_tibble(alldata))
}

run_human_clock = function(file, out_folder) {
  gse_id = str_extract(file, "GSE[0-9]+")
  if (file.exists(glue("{out_folder}/{gse_id}.csv"))) {
    print(glue("{gse_id} already exists"))
  } else {
    message(glue("Starts reading {gse_id}"))
    mat = geograbi::geograbi.read.gse.matrix(file)
    mat = mat$data |>
      as.data.frame() |>
      rownames_to_column("probe_id")
    if (nrow(mat) == 0) {
      print(glue("{gse_id} has no data"))
      return(NULL)
    } else {
      print(glue("{gse_id} has {nrow(mat)} probes"))
      g1 = do_grimage1(mat, keep = F)
      g2 = do_grimage2(mat, keep = T)

      dat = cbind(g2, g1[, -1])
      write_csv(dat, glue("{out_folder}/{gse_id}.csv"))
    }
  }
}

if (sys.nframe() == 0) {
  possible_run_human_clock = possibly(run_human_clock, otherwise = NULL)

  files = dir("./GEO_hs_dnam/series_matrix/", full.names = T)
  out_folder = "./GEO_hs_dnam/sample_age_grim"


  if (!dir.exists(out_folder)) {
    dir.create(out_folder)
  }

  # run_human_clock(files[1000], out_folder)
  pbmcapply::pbmclapply(
    files,
    possible_run_human_clock,
    out_folder = out_folder,
    mc.cores = 3
  )

  # rewrite to a for loop
  # geoid = "GSE140038"
  # i = which(str_detect(files, geoid))
  # for (i in 1:length(files)) {
  #   print(i)
  #   file = files[i]
  #   run_human_clock(file, out_folder)
  # }
}
