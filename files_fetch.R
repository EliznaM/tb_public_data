

# script to run on hemera

library(furrr)

## set up parallel

wrkrs <- 10
plan(multisession, workers = wrkrs) 


## Parameters

tmp <- readr::read_delim("geo_acc_list_manual_filter.txt")
access_nrs <- tmp$GEO
# access_nrs <- c("GSE89403", "GSE84076", "GSE107991")

sm_files <- paste0(access_nrs, "_series_matrix.txt")
# path_sm_files <- "sm_files"
path_sm_files <- "/home/emaasdorp/public_tb_data/sm_files"

# create dir

if(sum(grepl("sm_files", list.dirs(recursive = FALSE))) == 0){
  dir.create("sm_files")
}

## get files
# furrr functions

poss_getgeo <- purrr:::possibly(.f = GEOquery::getGEO, otherwise = "Error")

gse_list <- future_map(access_nrs, 
                       ~ poss_getgeo(.x, destdir = path_sm_files, AnnotGPL = FALSE))

# gse_list <- future_map(access_nrs, 
#                        ~ GEOquery::getGEO(.x, destdir = path_sm_files, AnnotGPL = FALSE))

# errors

sm_fldr_list <- list.files("sm_files")
sm_fldr_list <- gsub("_series_matrix.txt.gz", "", sm_fldr_list[stringr::str_detect(sm_fldr_list, "GSE")])
sm_fldr_list <- gsub("-GPL[0-9]+", "", sm_fldr_list)


err_acc_nrs <- access_nrs[!access_nrs %in% sm_fldr_list]

saveRDS(gse_list, "gse_list.rds")
saveRDS(err_acc_nrs, "err_acc_nrs.rds")

## Use phenoData.  Remove redundant columns and make long set. 
  # sort out lists in list
x <- which(purrr::map_dbl(gse_list, ~length(.x)) != 1)
if(length(x != 0)){
gse_list <- unlist(gse_list, recursive = FALSE)
}


# pData(gse_list[[1]][[1]])

pheno <- future_map2(gse_list, access_nrs,
                     ~pData(.x[[1]]) %>%
                       as_tibble() %>% 
                       mutate(geo_set = .y) %>% 
                       pivot_longer(-c(geo_accession, geo_set), names_to = "key", values_to = "value") %>% 
                       filter(!str_detect(key, "data_processing") &
                                !str_detect(value, "data_processing"))) %>% 
  bind_rows()

saveRDS(pheno, "pheno.rds")

