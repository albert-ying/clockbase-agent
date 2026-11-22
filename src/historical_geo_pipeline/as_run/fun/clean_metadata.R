require(geograbi)

## Variant Function --------------------------------------------------------------------------------------
create_variants = function(fields_list) {
  variants <- c()
  for (field in fields_list) {
    space_variants <- field
    if (grepl(" ", field)) {
      space_variants <- c(space_variants, gsub(" ", "_", field), gsub(" ", "-", field), gsub(" ", "", field))
    }
    for (space_variant in space_variants) {
      variants <- c(
        variants, space_variant, paste0("characteristic [", space_variant, "]"),
        paste0(
          "characteristic[", space_variant, "]"
        ), paste0("characteristics [", space_variant, "]"), paste0("characteristics[", space_variant, "]"),

        paste0(
          "comment [", space_variant, "]"
        ), paste0("comment[", space_variant, "]"), paste0("comments [", space_variant, "]"),

        paste0(
          "comments[", space_variant, "]"
        ), paste0("factorvalue[", space_variant, "]"), paste0("factor value[", space_variant, "]"),

        paste0(
          "factorvalue [", space_variant, "]"), paste0(
            "factor value [", space_variant, "]"
          ), paste0("sample_", space_variant), paste0("sample_host", space_variant),

        paste0("sample_sample_", space_variant)
      )
    }
  }
  variants
}

## Define the variant dictionary -------------------------------------------------------------------
sex_fields = c(
    "sex",
    "gender",
    "subject gender",
    "subject sex",
    # This looks reduntant, but there are some samples which use
    # Characteristic[Characteristic[sex]]
    "characteristic [sex]",
    "characteristics [sex]"
)
sex_repl = c("^sex ", " sex$", "gender ", " gender")

age_fields = c(
    "age",
    "patient age",
    "age of patient",
    "age (years)",
    "age (yrs)",
    "age (months)",
    "age (days)",
    "age (hours)",
    "age at diagnosis",
    "age at diagnosis years",
    "age at diagnosis months",
    "age at diagnosis days",
    "age at diagnosis hours",
    "characteristic [age]",
    "characteristics [age]"
)
age_repl = c("^age ", " age$")

specimen_part_fields = c(
    # AE
    "organism part",
    "cell type",
    "tissue",
    "tissue type",
    "tissue source",
    "tissue origin",
    "source tissue",
    "tissue subtype",
    "tissue/cell type",
    "tissue region",
    "tissue compartment",
    "tissues",
    "tissue of origin",
    "tissue-type",
    "tissue harvested",
    "cell/tissue type",
    "tissue subregion",
    "organ",
    "characteristic [organism part]",
    "characteristics [organism part]",
    # SRA
    "cell_type",
    "organismpart",
    # GEO
    "isolation source",
    "tissue sampled",
    "cell description"
)


genetic_information_fields = c(
    "strain/background",
    "strain",
    "strain or line",
    "background strain",
    "genotype",
    "genetic background",
    "genetic information",
    "genotype/variation",
    "ecotype",
    "cultivar",
    "strain/genotype"
)

disease_fields = c(
    "disease",
    "disease state",
    "disease status",
    "diagnosis",
    "infection with",
    "sample type"
)

disease_stage_fields = c(
    "disease state",
    "disease staging",
    "disease stage",
    "grade",
    "tumor grade",
    "who grade",
    "histological grade",
    "tumor grading",
    "disease outcome",
    "subject status"
)

cell_line_fields = c(
    "cell line",
    "sample strain"
)

treatment_fields = c(
    "treatment",
    "treatment group",
    "treatment protocol",
    "drug treatment",
    "clinical treatment",
    "drug",
    "drug name",
    "transfection"
)

race_fields = c(
    "race",
    "ethnicity",
    "race/ethnicity"
)

subject_fields = c(
    # AE
    "subject",
    "subject id",
    "subject/sample source id",
    "subject identifier",
    "human subject anonymized id",
    "individual",
    "individual identifier",
    "individual id",
    "patient",
    "patient id",
    "patient identifier",
    "patient number",
    "patient no",
    "donor id",
    "donor",
    # SRA
    "sample_source_name"
)

developmental_stage_fields = c(
    "developmental stage",
    "development stage",
    "development stages"
)

compound_fields = c(
    "compound",
    "compound1",
    "compound2",
    "compound name",
    "drug",
    "drugs",
    "immunosuppressive drugs"
)

time_fields = c(
    "time",
    "initial time point",
    "start time",
    "stop time",
    "time point",
    "sampling time point",
    "sampling time",
    "time post infection"
)

field.list = list(
  harmonized.sex = create_variants(sex_fields),
  harmonized.age = create_variants(age_fields),
  harmonized.specimen_part = create_variants(specimen_part_fields),
  harmonized.genetic_information = create_variants(genetic_information_fields),
  harmonized.disease = create_variants(disease_fields),
  harmonized.disease_stage = create_variants(disease_stage_fields),
  harmonized.cell_line = create_variants(cell_line_fields),
  harmonized.treatment = create_variants(treatment_fields),
  harmonized.race = create_variants(race_fields),
  harmonized.subject = create_variants(subject_fields),
  harmonized.developmental_stage = create_variants(developmental_stage_fields),
  harmonized.compound = create_variants(compound_fields),
  harmonized.time = create_variants(time_fields)
)

## Define Functions ------------------------------------------------------------------------------------
harmonize_value = function(field_name, value) {
  value = as.character(value)
  if (field_name == "age") {
    readr::parse_number(value)
  } else if (field_name == "sex") {
    purrr::map_chr(value, ~{
      if (tolower(.x) %in% c("f", "female", "woman")) {
        "female"
      } else if (tolower(.x) %in% c("m", "male", "man")) {
        "male"
      } else {
        tolower(.x)
      }
    })
  } else if (field_name == "subject") {
    as.numeric(as.factor(value))
  } else {
    value
  }
}
# test
if (FALSE) {
  field_name <- "age"
  value <- c("1.5 years", "2 years", "3 months", "4 days", "5 hours", "6")
  harmonize_value(field_name, value)
  value = c(1.2, 2, 3.1313)
  harmonize_value(field_name, value)
  field_name = "sex"
  value = c("M", "F", "1")
  harmonize_value(field_name, value)
}


pos_str_detect = purrr::possibly(str_detect, otherwise = FALSE)

harmonize_field = function(sample_metadata, field = field.list, debug = FALSE) {
  field_name = substr(names(field), 12, nchar(names(field)))
  unmatched_metadata = sample_metadata
  harmonized_meta = matrix(NA, nrow = nrow(sample_metadata), ncol = length(field)) |>
    as_tibble()
  colnames(harmonized_meta) = names(field)
  if (debug) message("Harmonizing ")
  match_id_vec = c()
  for (j in seq_along(field)) {
    matched_id = NA
    for (i in seq_along(sample_metadata)) {
      lower_key = tolower(colnames(sample_metadata)[i])
      if (lower_key %in% field[[j]]) {
        if (debug) message("Found ", lower_key, " in ", names(field)[j])
        harmonized_value = harmonize_value(field_name[[j]], sample_metadata[[i]])
        harmonized_meta[[j]] = harmonized_value
        matched_id = i # record the matched id
      } else if (
        any(map_lgl(rex::rex(field[[j]]), ~pos_str_detect(lower_key, .x))) &&
        all(is.na(harmonized_meta[[j]]))
      ) {
        if (debug) message("Found pattern of", lower_key, " in ", names(field)[j])
        harmonized_value = harmonize_value(field_name[[j]], sample_metadata[[i]])
        harmonized_meta[[j]] = harmonized_value
        matched_id = i
      } else {
        harmonized_value = NA
      }
    }
    if (!is.na(matched_id)) {
      match_id_vec = c(match_id_vec, matched_id)
    }
  }
  if (length(match_id_vec) > 0) {
    unmatched_metadata = unmatched_metadata[, -match_id_vec]
  }
  if (ncol(unmatched_metadata) > 0) {
    message("Unmatched metadata:")
    message(colnames(unmatched_metadata))
    harmonized_meta = bind_cols(unmatched_metadata, harmonized_meta) |>
      select(GSM, title, source_name, starts_with("harmonized"), everything())
  } else {
    message("All metadata are matched.")
  }
  return(harmonized_meta)
}

extract_metadata = function(sample_metadata) { ## Extract the key info from raw metadata
  char_meta = geograbi.extract.characteristics(sample_metadata)
  if (nrow(char_meta) > 0) {
    sample_metadata |>
      select(GSM, title, source_name = source_name_ch1) |>
      bind_cols(char_meta)
  } else {
    message("No metadata found.")
    sample_metadata |>
      select(GSM, title, source_name = source_name_ch1)
  }
}

simplify_metadata = function(sample_metadata) {
  #* Remove the column with all NAs
  sample_metadata[!map_lgl(sample_metadata, ~ all(is.na(.)))]
}

# test
if (FALSE) {
  sample_metadata = read_csv("./GEO_hs_dnam/metadata/GSE101961.csv")
  sample_metadata = sample_metadata |>
    extract_metadata() |>
    harmonize_field(debug = T) |>
    simplify_metadata()
}

## ZOOMA (not implemented yet) ---------------------------------------------
library(httr)
library(jsonlite)

term_to_top_match = function(term) {  
  term <- gsub("[^[:alnum:]]", "+", term)
# Define the ZOOMA API endpoint
  endpoint <- "http://www.ebi.ac.uk/spot/zooma/v2/api/services/annotate"

# Send a GET request to the ZOOMA API
  response <- GET(endpoint, query = list(propertyValue = term))

# Print the response status code
  print(status_code(response))

# Parse the response content as JSON
  json_content <- fromJSON(content(response, as = "text"))

  res_df = json_content$annotatedProperty |>
    mutate(confidence = json_content$confidence) |>
    as_tibble()
  res_df$propertyValue[[1]]
}

possible_term_to_top_match = possibly(term_to_top_match, otherwise = NULL)
term_to_top_match_mc = function(term, ncore) {  
  pbmcapply::pbmclapply(
    term,
    possible_term_to_top_match,
    mc.cores = ncore
  ) |>
    unlist()
}

get_desc = function(geo, retry = 50) {
  require(rvest)
  webpage = str_c("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", geo)
  content = NULL
  res = NULL
  while (is.null(content) & retry > 0) {
    tryCatch({
      content = rvest::read_html(webpage)
      content_text = as.character(content)
      title = str_extract(content_text, 'Title</td>\n<td style=\"text-align: justify\">(.*)</td>') |>
        str_replace_all('Title</td>\n<td style=\"text-align: justify\">', '') |>
        str_replace_all('</td>', '')
      summary = str_extract(content_text, 'Summary</td>\n<td style=\"text-align: justify\">(.*)<br>\n</td>') |>
        str_replace_all('Summary</td>\n<td style=\"text-align: justify\">', '') |>
        str_replace_all('<br>\n</td>', '')
        res = tibble(accession = geo, title = title, summary = summary)
    }, error = function(e) {
      message("Error: ", e)
      retry = retry - 1
      message("Retry: ", retry)
    })
  }
  return(res)
}

## modify the function to add the retry option
geograbi.retrieve.datasets <- function(gpl = NULL, retmax = 9000, max_retries = 10) {
  url <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
  lines <- NULL
  retstart <- 0
  done <- F
  retries <- 0
  if (is.null(gpl)) {
    term <- "gse[ETYP]"
  } else {
    term <- paste0(gpl, "[ACCN]+AND+gse[ETYP]")
    term <- paste0(term, collapse = "+OR+")
  }
  while (!done) {
    cat(date(), "retrieving dataset UID", retstart + 1, "...\n")
    result <- tryCatch({
      gds <- readLines(paste0(
        url, "esearch.fcgi", "?db=gds",
        "&term=", term, "&retstart=", format(retstart, scientific = F),
        "&retmax=", format(retmax, scientific = F), "&usehistory=y"
      ))
      gds <- xmlToList(xmlParse(gds, asText = T))
      id.n <- length(gds$IdList)
      new <- readLines(paste0(
        url, "efetch.fcgi", "?db=gds",
        "&query_key=", gds$QueryKey, "&WebEnv=", gds$WebEnv,
        "&retstart=", gds$RetStart, "&retmax=", gds$RetMax
      ))
      new.n <- length(grep("^[0-9]+[.]+ ", new))
      if (new.n != id.n) {
        warning("Retrieved only", new.n, "of", id.n, "summaries.")
      }
      lines <- c(lines, new)
      cat(date(), "retrieved", new.n, "\n")
      retstart <- retstart + new.n
      done <- new.n < retmax
      retries <- 0
      TRUE
    }, error = function(e) {
      cat("Error occurred:", e$message, "\n")
      retries <- retries + 1
      if (retries >= max_retries) {
        cat("Max retries reached. Exiting.\n")
        return(FALSE)
      } else {
        cat("Retrying...\n")
        return(TRUE)
      }
    })
    if (!result) {
      break
    }
  }
  lines <- gsub("Platforms:", "Platform:", lines)
  lines <- gsub("[[:space:]]+Platform:", "XXXXPlatform:", lines)
  lines <- gsub("[[:space:]]+ID:", "XXXXID:", lines)
  lines <- unlist(strsplit(lines, "XXXX"))
  lines <- gsub(
    "Series[[:space:]]+Accession:", "Accession:",
    lines
  )
  lines <- gsub("Datasets:", "Dataset:", lines)
  columns <- unique(sub("^([[:alpha:]]+:).*", "\\1", lines[grep(
    "^[[:alpha:]]+:",
    lines
  )]))
  idx <- grep("^[0-9]+[.]+ ", lines)
  rows <- c(tail(idx, -1), length(lines) + 1) - idx
  idx <- rep(1:length(idx), rows)
  filter.and.insert <- function(x, pattern, idx) {
    ret <- rep("", max(idx))
    p.idx <- grep(pattern, x)
    idx <- idx[p.idx]
    x <- x[p.idx]
    ret[idx] <- sub(pattern, "", x)
    ret
  }
  ret <- data.frame(title = filter.and.insert(
    lines, "^[0-9]+[.]+ ",
    idx
  ), description = filter.and.insert(
    lines, "^\\(Submitter supplied\\) ",
    idx
  ), stringsAsFactors = F)
  for (column in columns) {
    ret[[sub(":", "", column)]] <- filter.and.insert(
      lines,
      paste0("^", column, "[^A-Za-z0-9]*"), idx
    )
  }
  colnames(ret) <- tolower(colnames(ret))
  ret$samples <- as.integer(sub(
    ".* ([0-9]+) Sample[s]*$",
    "\\1", ret$platform
  ))
  ret$platform <- sub(" [0-9]+ Sample[s]*$", "", ret$platform)
  ret
}
download.series.files = function(path = ".", gses, gpls = NULL) {
    require(RCurl)
    filenames <- sapply(seq_along(gses), function(gse) {
        urls <- get.gse.matrix.url(gses[gse])
        if (!is.null(gpls)) {
            urls <- c(urls, gsub("_series_", paste0("-", gpls[gse], 
                "_series_"), urls[1]))
        }
        base::cat(base::date(), gses[gse], urls[1], "\n")
        for (url in urls) {
            destination <- file.path(path, basename(url))
            if (file.exists(destination)) 
                return(destination)
        }
        for (url in urls) {
            destination <- file.path(path, basename(url))
            tryCatch(download.file(url, destfile = destination), 
                error = function(e) {
                  unlink(destination)
                  print(e)
                })
            Sys.sleep(3)
            if (file.exists(destination)) 
                return(destination)
        }
        NA
    })
    filenames
}
