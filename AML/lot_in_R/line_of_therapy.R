lot_func <- function(tbl, allowed_days = 45, L1_days = 30) {
  require(data.table)
  require(lubridate)

  # input table consists of four columns
  # --------------------------------------------
  # | ptid | drug_name | drug_start | drug_end |
  # --------------------------------------------
  #
  # Note: drug end may be a calculated field. ie. drug_start + days supply
  
  # 1. Look for continuous episodes of each therapy per patient using allowable days
  #    and collapse using minimum drug_start and maximum drug end
  epi_tbl <- as.data.table(tbl)
  setorder(epi_tbl, ptid, drug_name, drug_start)
  epi_tbl[, day_diff := as.integer(drug_start) - as.integer(shift(drug_end)), by = .(ptid, drug_name)]
  epi_tbl[, flag := ifelse(day_diff > allowed_days, 1, 0), by = .(ptid, drug_name)][is.na(flag), flag := 0]
  epi_tbl[, csum := cumsum(flag), by = .(ptid, drug_name)][, `:=`(day_diff = NULL, flag = NULL)]
  epi_tbl[, `:=`(drug_start = min(drug_start), drug_end = max(drug_end)), by = .(ptid, drug_name, csum)]
  episode_tbl <- unique(epi_tbl[, .(ptid, drug_name, drug_start, drug_end)])
  setorder(episode_tbl, ptid, drug_start)
  
  # 2. Look for overlapping drugs using drug start and drug end dates and 
  #    combine on seperate rows using '+' seperator
  comb_tbl <- episode_tbl
  setorder(comb_tbl, ptid, drug_name, drug_start)
  comb_tbl[, `:=`(drug_start = ymd(drug_start), drug_end = ymd(drug_end))]
  
  comb_tbl[, nrow := 1:.N, by = ptid]
  comb_tbl <- comb_tbl[rep(seq_len(.N), (drug_end - drug_start)+1), .SD, by = .(ptid, nrow)]
  comb_tbl[, day := seq(min(drug_start), max(drug_end), 'day'), by = .(ptid, nrow)]
  comb_tbl <- comb_tbl[, .(ptid, drug_name, day)]
  comb_tbl <- unique(comb_tbl[, drug_name := paste0(sort(unique(drug_name)), collapse = " + "), keyby = .(ptid, day)])
  comb_tbl[, id := rleid(drug_name), by = .(ptid)]
  comb_tbl[, `:=`(drug_start = min(day), drug_end = max(day)), keyby = .(ptid, drug_name, id)]
  comb_tbl <- unique(comb_tbl[order(ptid, drug_start), .(ptid, drug_name, drug_start, drug_end)])
  setnames(comb_tbl, 'drug_name', 'treatment')

  # Count number of drugs per row 
  comb_tbl[, drug_count := stringr::str_count(treatment, "\\+") + 1]
  
  # 3. Asign line of therapy numbers to combination table using L1_days to
  #    identify line one drugs
  lot <- comb_tbl
  lot_min <- copy(lot)
  lot_min <- lot_min[, .(drug_start = min(drug_start)), by = ptid]
  lot_tbl <- data.table('ptid' = character(), 
                       'treatment' = character(),
                       'drug_count' = numeric(),
                       'treatment_start' = ymd(), 
                       'treatment_end' = ymd(), 
                       'lot' = numeric(),
                       'regimen' = character())
  
  for(i in 1:nrow(lot_min)){
    # 1. subset on ptid
    # 2. minimum start date + L1 days
    # 3. find all unique drugs in first L1 days
    # 4. Assign new line when new drug introduced
    lot_i <- lot[ptid == lot_min[i, ptid]]
    l1_day <- as.Date(lot_min[i, drug_start] + L1_days) 
    line_drugs <- unique(unlist(strsplit(lot_i[drug_start <= l1_day, treatment], ' \\+ ')))
    
    lot_num <- 1 # set lot with initial value
    for(j in 1:nrow(lot_i)){
      extra_drugs <- setdiff(unique(unlist(strsplit(lot_i[j, treatment], ' \\+ '))), line_drugs)
      if(length(extra_drugs) == 0){
        lot_temp <- data.table('ptid' = lot_i[j, ptid], 
                              'treatment' = lot_i[j, treatment],
                              'drug_count' = lot_i[j, drug_count], 
                              'treatment_start' = lot_i[j, drug_start], 
                              'treatment_end' = lot_i[j, drug_end], 
                              'lot' = lot_num,
                              'regimen' = paste(sort(line_drugs), collapse = " + "))
        lot_tbl <- rbind(lot_tbl, lot_temp)
      } else {
        lot_num <- lot_num + 1
        line_drugs <- unique(unlist(strsplit(lot_i[j, treatment], ' \\+ ')))
        lot_temp <- data.table('ptid' = lot_i[j, ptid], 
                              'treatment' = lot_i[j, treatment],
                              'drug_count' = lot_i[j, drug_count], 
                              'treatment_start' = lot_i[j, drug_start], 
                              'treatment_end' = lot_i[j, drug_end], 
                              'lot' = lot_num,
                              'regimen' = paste(sort(line_drugs), collapse = " + "))
        lot_tbl <- rbind(lot_tbl, lot_temp)    
      }
    }
    
  }
  
  lot_tbl[, max_lot := max(lot), by = ptid]
  
  setorder(lot_tbl, ptid, lot)
  return(lot_tbl) 
}
