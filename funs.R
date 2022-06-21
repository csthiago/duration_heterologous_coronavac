##############################################################################
# Name of file: auxiliary functions
# Original author(s): Thiago Cerqueira
# Latest update author (if not using version control) - thiago.c.silva@fiocruz.br
# Written/run onL: R Studio SERVER
# Version of R that the script was most recently run on: R 4.1.1
# Description of content: functions
# Approximate run time: Unknown
##############################################################################

# Define functions
fun_coef_gam <- function(mod){
  round(exp(cbind("Odds Ratio" = coef(mod), confint.default(mod, level = 0.95))), digits = 3) %>% 
    as.data.frame() %>%
    rownames_to_column(var = "term")}

fct_case_when <- function(...) {
  # uses dplyr::case_when but converts the output to a factor,
  # with factors ordered as they appear in the case_when's  ... argument
  args <- as.list(match.call())
  levels <- sapply(args[-1], function(f) f[[3]])  # extract RHS of formula
  levels <- levels[!is.na(levels)]
  factor(dplyr::case_when(...), levels=levels)
}

models_fit <- function(dataset) {
  mod <- bam(confirmado == "Confirmado" ~ vs_type2 +
    s(temporal_trend, bs = "cr") +
    s(idade, bs = "cr") +
    prev_infec +
    sexo +
    uf +
    diabetes +
    obesidade +
    imunossupressao +
    drespiratoria +
    dcardiaca +
    drc +
    q_measure_1f_12 +
    capital_not,
  family = binomial, nthreads = 24, discrete = T, data = dataset
  )
tab <- fun_coef_gam(mod) %>% filter(str_detect(term,"vs_type2"))
  table_results <- tab %>%
    rename(conf.low = `2.5 %`,
           conf.high= `97.5 %`,
           or = `Odds Ratio`) %>%
    mutate(
      term = str_remove(term, "vs_type2"),
      age_group = paste(min(dataset$idade), max(dataset$idade), sep = "-")
    )
  table_results
}

fun_table_ve <- function(dataset=overall,dataset2=age,title="Title"){
  dataset %>% bind_rows(dataset2 %>%  select(fit) %>% unnest(names_repair = "check_unique",cols = c(fit))) %>% 
    mutate(across(or:conf.high, ~ (1 - .x) * 100)) %>%
    mutate(across(or:conf.high, ~ if_else(.x < (-1000) | .x > 1000, Inf, .x))) %>%
    mutate(type = case_when(
      str_detect(term, "CV$") ~ "CV-Booster",
      str_detect(term, "BNT162b2$") ~ "BNT162b2-Booster",
      str_detect(term, "BNT162b2$|CV$", negate = T) ~ "No booster"
    )) %>%
    mutate(
      term = str_remove(term, "CV_"),
      term = str_remove(term, "_CV$"),
      term = str_remove(term, "_BNT162b2$")
    ) %>%
    group_by(type, age_group) %>%
    mutate(term = fct_relevel(
      factor(term),
      "v1_0:1",
      "v1_2",
      "v2_0:13",
      "v2_14:90",
      "v2_91:180",
      "v2_181",
      "v3_0:13",
      "v3_14:30",
      "v3_31:60",
      "v3_61:90",
      "v3_91:120",
      "v3_121"
    )) %>%
    arrange(term, type, age_group) %>%
    gt() %>%
    fmt_number(columns = 2:4, decimals = 1) %>%
    tab_style(
      style = list(
        cell_text(align = "center")
      ),
      locations = cells_column_labels(columns = everything())
    ) %>%
    tab_style(
      style = cell_text(color = "black", weight = "bold", align = "left"),
      locations = cells_row_groups()
    ) %>%
    cols_merge(
      columns = c(conf.high, conf.low),
      pattern = "{1} &mdash; {2}"
    ) %>%
    cols_label(
      or = md("**VE**"),
      conf.high = md("**95% CI**"),
      term = md("**Period**")
    ) %>%
    opt_row_striping() %>% 
    tab_header(title = title)
}

fun_table_or <- function(dataset=overall,dataset2=age,title="Title"){
  dataset %>% bind_rows(dataset2 %>%  select(fit) %>% unnest(names_repair = "check_unique",cols = c(fit))) %>% 
    mutate(across(or:conf.high, ~ if_else(.x < (-1000) | .x > 1000, Inf, .x))) %>%
    mutate(type = case_when(
      str_detect(term, "CV$") ~ "CV-Booster",
      str_detect(term, "BNT162b2$") ~ "BNT162b2-Booster",
      str_detect(term, "BNT162b2$|CV$", negate = T) ~ "No booster"
    )) %>%
    mutate(
      term = str_remove(term, "CV_"),
      term = str_remove(term, "_CV$"),
      term = str_remove(term, "_BNT162b2$")
    ) %>%
    group_by(type, age_group) %>%
    mutate(term = fct_relevel(
      factor(term),
      "v2_181",
      "v1_0:1",
      "v1_2",
      "v2_0:13",
      "v2_14:90",
      "v2_91:180",
      "v3_0:13",
      "v3_14:30",
      "v3_31:60",
      "v3_61:90",
      "v3_91:120",
      "v3_121"
    )) %>%
    arrange(term, type, age_group) %>%
    gt() %>%
    fmt_number(columns = 2:4, decimals = 3) %>%
    tab_style(
      style = list(
        cell_text(align = "center")
      ),
      locations = cells_column_labels(columns = everything())
    ) %>%
    tab_style(
      style = cell_text(color = "black", weight = "bold", align = "left"),
      locations = cells_row_groups()
    ) %>%
    cols_merge(
      columns = c(conf.high, conf.low),
      pattern = "{1} &mdash; {2}"
    ) %>%
    cols_label(
      or = md("**OR**"),
      conf.high = md("**95% CI**"),
      term = md("**Period**")
    ) %>%
    opt_row_striping() %>% 
    tab_header(title = title)
}


fun_plot_ve_type <- function(dataset,dataset2){
  dataset %>% bind_rows(dataset2 %>%  select(fit) %>% unnest(names_repair = "check_unique",cols = c(fit))) %>% 
    mutate(across(or:conf.high, ~ (1 - .x) * 100)) %>%
    mutate(across(or:conf.high, ~ if_else(.x < (-1000) | .x > 1000, Inf, .x))) %>%
    mutate(type = case_when(
      str_detect(term, "CV$") ~ "CV-Booster",
      str_detect(term, "BNT162b2$") ~ "BNT162b2-Booster",
      str_detect(term, "BNT162b2$|CV$", negate = T) ~ "No booster"
    )) %>%
    mutate(
      term = str_remove(term, "CV_"),
      term = str_remove(term, "_CV$"),
      term = str_remove(term, "_BNT162b2$")
    ) %>% filter(type != "No booster") %>% 
    group_by(type, age_group) %>%
    mutate(term = fct_relevel(
      factor(term),
      "v3_0:13",
      "v3_14:30",
      "v3_31:60",
      "v3_61:90",
      "v3_91:120",
      "v3_121"
    )) %>%
    arrange(term, type, age_group) %>% 
    ggplot(aes(x = term, y = or, ymin = conf.low, ymax = conf.high)) +
    geom_pointrange() +
    geom_hline(yintercept = 0, lty = 2) +
    labs(y = "VE (%)", x = "", title="Compared to Unvaccinated") +
    facet_wrap(~ type + age_group, ncol = 4, scales = "free_x") +
    theme(axis.text.x = element_text(angle = 90)) +
    theme(strip.background = element_rect(fill = "gray90")) +
    theme(strip.text = element_text(colour = "black")) +
    coord_cartesian(ylim = c(-20, 100))
}

fun_plot_or_type <- function(dataset,dataset2){
  dataset %>% bind_rows(dataset2 %>%  select(fit) %>% unnest(names_repair = "check_unique",cols = c(fit))) %>% 
    mutate(across(or:conf.high, ~ if_else(.x < (-1000) | .x > 1000, Inf, .x))) %>%
    mutate(type = case_when(
      str_detect(term, "CV$") ~ "CV-Booster",
      str_detect(term, "BNT162b2$") ~ "BNT162b2-Booster",
      str_detect(term, "BNT162b2$|CV$", negate = T) ~ "No booster"
    )) %>%
    mutate(
      term = str_remove(term, "CV_"),
      term = str_remove(term, "_CV$"),
      term = str_remove(term, "_BNT162b2$")
    ) %>% filter(type != "No booster") %>% 
    group_by(type, age_group) %>%
    mutate(term = fct_relevel(
      factor(term),
      "v3_0:13",
      "v3_14:30",
      "v3_31:60",
      "v3_61:90",
      "v3_91:120",
      "v3_121"
    )) %>%
    arrange(term, type, age_group) %>% 
    ggplot(aes(x = term, y = or, ymin = conf.low, ymax = conf.high)) +
    geom_pointrange() +
    geom_hline(yintercept = 1, lty = 2) +
    labs(y = "OR", x = "", title="Compared to 2nd dose 181+") +
    facet_wrap(~ type + age_group, ncol = 4, scales = "free_x") +
    theme(axis.text.x = element_text(angle = 90)) +
    theme(strip.background = element_rect(fill = "gray90")) +
    theme(strip.text = element_text(colour = "black")) +
    coord_cartesian(ylim = c(0, 2))
}


extract_variant_period <- function(path,country){
  data_dt <- path %>%
    readxl::excel_sheets() %>%
    set_names() %>%
    map_df(~ readxl::read_excel(path = path, sheet = .x), .id = "sheet") %>%
    janitor::clean_names() %>%
    fill(x1, .direction = "down") %>%
    filter(!is.na(x1)) %>%
    mutate(count_total = coalesce(!!!select(., starts_with("voc")))) %>%
    relocate(count_total, .after = x1) %>%
    select(!matches("^voc|^voi"), -readme_variants_download_package) %>%
    pivot_longer(cols = -c(sheet:count_total), names_to = "date", values_to = "value") %>%
    filter(x1 == country, !is.na(count_total)) %>%
    pivot_wider(names_from = count_total, values_from = value) %>%
    filter(!is.na(count)) %>%
    mutate(
      date = str_remove(date, "^x"),
      date = lubridate::ymd(date)
    ) %>%
    pivot_wider(names_from = sheet, values_from = count) %>%
    mutate(across(3:last_col(), ~ .x * 100 / total)) %>% 
    select(-x1) %>%
    janitor::clean_names() %>%
    group_by(date) %>%
    left_join(
      select(., -total) %>%
        pivot_longer(names_to = "voc", values_to = "pct", cols = voc_omicron:voc_gamma) %>%
        group_by(date) %>%
        slice(which.max(pct)),
      by = "date"
    ) %>%
    select(date, voc, pct) %>%
    mutate(
      pct = if_else(pct < 51, NA_real_, pct),
      voc = if_else(is.na(pct), "wt", voc)
    ) %>%
    group_by(voc) %>%
    arrange(date) %>%
    slice(c(1, n())) %>%
    mutate(startend = c("start", "end")) %>%
    mutate(date = if_else(voc == "wt" & startend == "start", as.Date("2020-02-20"), date)) %>% 
    pivot_wider(id_cols = voc,names_from = startend, values_from = c(date,pct)) %>% 
    mutate(date_start=date_start-6) %>% 
    group_by(voc) %>% mutate(intervalo=interval(date_start,date_end)) %>% 
    column_to_rownames(var = "voc")
  return(data_dt)}


# 
# hp_table <- function(x){
#   gt(x) %>% 
#     as_raw_html()
# }
# 
# mod_cv %>% bind_rows(infection %>%  select(fit) %>% unnest(names_repair = "check_unique",cols = c(fit))) %>% 
#   mutate(across(or:conf.high, ~ (1 - .x) * 100)) %>%
#   mutate(across(or:conf.high, ~ if_else(.x < (-1000) | .x > 1000, Inf, .x))) %>%
#   mutate(type = case_when(
#     str_detect(term, "CV$") ~ "CV-Booster",
#     str_detect(term, "BNT162b2$") ~ "BNT162b2-Booster",
#     str_detect(term, "BNT162b2$|CV$", negate = T) ~ "No booster"
#   )) %>%
#   mutate(
#     term = str_remove(term, "CV_"),
#     term = str_remove(term, "_BNT162b2$")
#   ) %>%
#   group_by(type, age_group) %>%
#   mutate(term = fct_relevel(
#     factor(term),
#     "v1_0:1",
#     "v1_2",
#     "v2_0:13",
#     "v2_14:90",
#     "v2_91:180",
#     "v2_181",
#     "v3_0:13",
#     "v3_14:30",
#     "v3_31:60",
#     "v3_61:90",
#     "v3_91:120",
#     "v3_121"
#   )) %>%
#   arrange(term, type, age_group) %>% 
#   group_map(~hp_table(.x)) %>% 
#   setNames(., c("High mileage", "Low mileage")) %>% 
#   data.frame(.) %>% 
#   gt() %>% 
#   fmt_markdown(columns = TRUE)


