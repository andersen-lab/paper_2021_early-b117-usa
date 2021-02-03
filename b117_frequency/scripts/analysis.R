library(tidyverse)
library(lubridate)
library(zoo)
library(gridExtra)
library(sf)

## Combine sequencing and testing data
helix_b117 <- read_tsv("../data/covid_baseline_for_b117_paper.20210127_update.txt") %>%
    select(state, collection_date, n_b117, n_sgtf_seq)

helix_sgtf <- read_tsv("../data/covid_baseline_for_b117_paper.20210201_klados20211029_phyloseq.txt") %>%
    select(state, collection_date, n, n_sgtf)

helix_metadata <- left_join(helix_sgtf, helix_b117, by=c("state", "collection_date"))


tmp <- helix_metadata %>%
    group_by(collection_date) %>%
    summarise(n_sgtf = sum(n_sgtf), n = sum(n)) %>%
    mutate(state = "USA")

tmp <- bind_rows(tmp, helix_metadata)
states_gt_500 <- tmp %>% group_by(state) %>% summarise(n = sum(n), n_sgtf = sum(n_sgtf)) %>% filter(n > 500 & n_sgtf > 0) %>% select(state) %>% as_vector()

p <- tmp %>%
    mutate(
        sgtf_pct = (n_sgtf/n) * 100
    ) %>%
    filter(state %in% states_gt_500) %>%
    group_by(state) %>%
    group_map(~{
        .x %>%
            ggplot(aes(collection_date, sgtf_pct)) + geom_col() + theme_bw() + ggtitle(.y) + ylab("% SGTF") + xlab("Date") +
            scale_x_date(limits=c(as.Date("2020-12-15"), as.Date("2021-02-01")))
    })
pdf("../plots/sgtf_drop_all_states.pdf", w = 15, h =12)
do.call("grid.arrange", c(p, ncol=5))
dev.off()

tmp <- helix_metadata %>%
    group_by(collection_date) %>%
    summarise(n_b117 = sum(n_b117, na.rm=TRUE), n_sgtf_seq = sum(n_sgtf_seq, na.rm=TRUE), n =sum(n, na.rm=TRUE)) %>%
    mutate(state = "USA")

tmp <- bind_rows(tmp, helix_metadata)
states_gt_500 <- tmp %>% group_by(state) %>% summarise(n = sum(n, na.rm=TRUE), n_sgtf_seq = sum(n_sgtf_seq, na.rm=TRUE), n_b117 = sum(n_b117, na.rm=TRUE)) %>% filter(n > 500 & n_sgtf_seq > 0 & n_b117 > 0) %>% select(state) %>% as_vector()

p <- tmp %>%
    mutate(
        b117_sgtf_pct = (n_b117/n_sgtf_seq) * 100
    ) %>%
    filter(state %in% states_gt_500) %>%
    group_by(state) %>%
    group_map(~{
        .x %>%
            ggplot(aes(collection_date, b117_sgtf_pct)) + geom_col(fill="indianred") + theme_bw() + ylab("% B.1.1.7 in sequenced SGTF") + xlab("Date") + ggtitle(.y) +             scale_x_date(limits=c(as.Date("2020-12-15"), as.Date("2021-02-01"))) + scale_y_continuous(limits=c(0,100))
    })
pdf("../plots/b117_sgtf_pct.pdf", w = 15, h =12)
do.call("grid.arrange", c(p, ncol=5))
dev.off()

## Calculate rolling means
helix_metadata_rolling_avg <- helix_metadata %>%
    filter(state %in% c("CA", "FL", "GA")) %>%
    group_by(state) %>%
    group_modify(~{
        tmp <- .x
        last_date <- .x %>% filter(n_b117 > 0) %>% tail(n = 1) %>% select(collection_date) %>% first()
        first_date <- .x %>% filter(n_b117 > 0) %>% head(n = 1) %>% select(collection_date) %>% first()
        avg_last_5day_b117_seq_pct = tmp %>% filter(n_b117 > 0 & n_sgtf_seq > 0) %>% filter(collection_date >= last_date - 5) %>% summarise(b117_sgtf_pct = mean(n_b117/n_sgtf_seq)) %>% select(b117_sgtf_pct) %>% last() #Averge of last 5 days of B117
        tmp <- tmp %>%
            mutate(
                n_b117 = ifelse(collection_date > last_date, NA, n_b117)
            )
        tmp <- tmp %>%
            mutate(
                type = ifelse(is.na(n_b117), "extrapolated", "data"),
                b117_est_pct_raw = ifelse(!is.na(n_b117), (n_b117/n_sgtf_seq) * (n_sgtf/n), avg_last_5day_b117_seq_pct * (n_sgtf/n))
            )
        tmp %>%
            filter(!is.na(b117_est_pct_raw)) %>%
            group_by(collection_date) %>%
            group_modify(~{
                tmp_date <- .x
                rolling_values <- tmp %>%
                    filter(collection_date >= .y - 2 & collection_date <= .y + 2) %>%
                    summarise(
                        b117_est_pct = mean(b117_est_pct_raw, na.rm=TRUE)
                    )
                tmp_date %>%
                    mutate(
                        b117_est_pct = rolling_values$b117_est_pct
                    )
            })
    })

## Plot rolling means
rm_ca_p <- helix_metadata_rolling_avg %>%
    filter(state == "CA") %>%
    ggplot() +
    geom_col(aes(collection_date, b117_est_pct_raw, fill = type)) +
    geom_line(aes(collection_date, b117_est_pct)) +
    geom_point(aes(collection_date, b117_est_pct)) +
    theme_bw() +
    scale_y_continuous(limits = c(0, 0.06), name="Proportion of B117 cases") + xlab("Days") +
    scale_x_date(limits=c(as.Date("2020-12-15"), as.Date("2021-02-01")))

rm_fl_p <- helix_metadata_rolling_avg %>%
    filter(state == "FL") %>%
    ggplot() +
    geom_col(aes(collection_date, b117_est_pct_raw, fill = type)) +
    geom_line(aes(collection_date, b117_est_pct)) +
    geom_point(aes(collection_date, b117_est_pct)) +
    theme_bw() +
    scale_y_continuous(limits = c(0, 0.06), name="Proportion of B117 cases") + xlab("Days") +
    scale_x_date(limits=c(as.Date("2020-12-15"), as.Date("2021-02-01")))

## USA
## Calculate rolling means
helix_metadata_us <- helix_metadata %>%
    group_by(collection_date) %>%
    summarise(n_b117 = sum(n_b117), n_sgtf_seq = sum(n_sgtf_seq), n_sgtf = sum(n_sgtf), n = sum(n))

last_date <- helix_metadata_us %>% filter(n_b117 > 0) %>% tail(n = 1) %>% select(collection_date)
first_date <- helix_metadata_us %>% filter(n_b117 > 0) %>% head(n = 1) %>% select(collection_date)
avg_last_5day_b117_seq_pct = helix_metadata_us %>% filter(n_b117 > 0 & n_sgtf_seq > 0) %>% filter(collection_date >= last_date - 5) %>% summarise(b117_sgtf_pct = mean(n_b117/n_sgtf_seq)) %>% select(b117_sgtf_pct) %>% last() #Averge of last 5 days of B117

helix_metadata_us <- helix_metadata_us %>%
    mutate(
        n_b117 = ifelse(collection_date > last_date, NA, n_b117)
    )

helix_metadata_us <- helix_metadata_us %>%
    mutate(
        type = ifelse(is.na(n_b117), "extrapolated", "data"),
        b117_est_pct_raw = ifelse(!is.na(n_b117), (n_b117/n_sgtf_seq) * (n_sgtf/n), avg_last_5day_b117_seq_pct * (n_sgtf/n))
    )
helix_metadata_us_rolling_avg <- helix_metadata_us %>%
    filter(!is.na(b117_est_pct_raw)) %>%
    group_by(collection_date) %>%
    group_modify(~{
        tmp_date <- .x
        rolling_values <- helix_metadata_us %>%
            filter(collection_date >= .y - 2 & collection_date <= .y + 2) %>%
            summarise(
                b117_est_pct = mean(b117_est_pct_raw, na.rm=TRUE)
            )

        tmp_date %>%
                    mutate(
                        b117_est_pct = rolling_values$b117_est_pct
                    )
    })



## Plot rolling means
rm_usa_p <- helix_metadata_us_rolling_avg %>%
    ggplot() +
    geom_col(aes(collection_date, b117_est_pct_raw, fill = type)) +
    geom_line(aes(collection_date, b117_est_pct)) +
    geom_point(aes(collection_date, b117_est_pct)) +
    theme_bw() +
    scale_y_continuous(limits = c(0, 0.06), name="Proportion of B117 cases") + xlab("Days") +
    scale_x_date(limits=c(as.Date("2020-12-15"), as.Date("2021-02-01")))

svg("../plots/b117_rolling_mean.svg", w= 15, h = 5)
grid.arrange(rm_usa_p, rm_fl_p, rm_ca_p, ncol = 3, nrow =1)
dev.off()

helix_metadata_us_rolling_avg_pct %>%
    select(b117_sgtf_pct, sgtf_pct, b117_est_pct, collection_date) %>%
    gather(key, val, -collection_date) %>%
    ggplot(aes(collection_date, val, fill = key)) + geom_line() + geom_col() + facet_grid(key ~ ., scales="free") + theme_bw()
ggsave("../plots/pct_from_rolling_mean_usa.pdf", w = 10, h =10)

## Fit NLS
fit_plot_logistic_growth <- function(df, name){
    tmp <- df %>%
        arrange(collection_date)
    first_date <- tmp %>%
        filter(b117_est_pct > 0) %>%
        arrange(collection_date) %>%
        ungroup()  %>%
        select(collection_date) %>%
        head(n = 1) %>%
        first()
    ## x0 <- tmp %>%
    ##     filter(b117_est_pct > 0) %>%
    ##     arrange(collection_date) %>%
    ##     ungroup()  %>%
    ##     select(b117_est_pct) %>%
    ##     head(n = 1) %>%
    ##     first()
    tmp <- tmp %>%
        filter(collection_date >= first_date)
    tmp$ndays = (tmp$collection_date - tmp$collection_date[1]) %>% as.integer()
    fit <- nls(b117_est_pct ~ 1/(1+(((1/x0) - 1) * exp( -1 * r * ndays))), data = tmp, start=list(r = 0.08, x0 = 0.001))
    tmp_newdata <- data.frame(ndays = seq(0,300), collection_date = first_date + seq(0,300))
    tmp_newdata$predicted = predict(fit, newdata = tmp_newdata)
    above_50_date <- tmp_newdata %>% filter(predicted >= 0.5)  %>% select(collection_date) %>% first() %>% first()
    serial_interval <- c(5,6.5)
    print(coef(fit))
    ggplot() +
        geom_point(data = tmp, aes(collection_date, b117_est_pct, color = type)) +
        geom_line(data = tmp_newdata, aes(collection_date, predicted)) +
        geom_text(data = tmp, x = first_date, y = 1, label = paste0("Logistic growth rate = ", round(coef(fit)[1], 2), " per day"), hjust=0) +
        geom_text(data=tmp, x = first_date, y = 0.9, label = paste0("Increase in transmission = ", paste(paste0(round(serial_interval * coef(fit)[1] *100, 1), "%"), collapse = " - ")), hjust=0) +
        ## geom_text(x = first_date, y = max(tmp$b117_est_pct, na.rm=T), label = paste0("Doubling time ~ 1 week "), hjust=0) +
        geom_text(data = tmp, x = first_date, y = 0.8, label = paste0("Rough doubling time = ", round(log(2)/coef(fit)[1], 2)), hjust=0) +
        geom_text(data = tmp, x = above_50_date, y = 0.1, label = paste0("50% date = ", above_50_date), hjust=0) +
        geom_vline(xintercept = above_50_date, color="red") +
    theme_bw() + scale_x_date(limits=c(as.Date("2020-12-15"), as.Date("2021-02-15")), date_breaks="15 days", date_labels="%b %d") + scale_y_continuous(limits = c(0, 0.06), name="Proportion of B117 cases") + xlab("Days") + ggtitle(name)
}

p1 <- fit_plot_logistic_growth(helix_metadata_us_rolling_avg, "USA")
p2 <- helix_metadata_rolling_avg %>%
    filter(state == "FL") %>%
    fit_plot_logistic_growth("FL")
p4 <- helix_metadata_rolling_avg %>%
    filter(state == "CA") %>%
    fit_plot_logistic_growth("CA")

svg("../plots/fit_b117_rolling_mean_inset.svg", w= 20, h = 5)
## grid.arrange(p1, p2, p3, p4, ncol =2, nrow =2)
grid.arrange(p1, p2, p4, ncol = 3, nrow =1)
dev.off()

## Testing coverage
testing_counts <- helix_sgtf %>% group_by(state) %>% summarise(n = sum(n))
us_map <- st_read("..//geo/US_states.json")
testing_us_map <- left_join(us_map, testing_counts, by=c("state_id" = "state"))
testing_us_centroid <- testing_us_map %>%
    filter(!(state_id %in% c("GU", "PR", "HI", "AK"))) %>%
    st_centroid

testing_us_map %>%
    filter(!(state_id %in% c("GU", "PR", "HI", "AK"))) %>%
    ggplot() +
    geom_sf(fill="#FFFFFF") +
    geom_sf(data = testing_us_centroid, mapping = aes(size = n, fill = n), shape = 21, alpha = 0.6, color = "#333333") +
    scale_size_continuous(range = c(1, 20), breaks = c(10, 100, 1000, 10000, 50000), trans="log") +
    scale_fill_distiller(trans="log", palette="YlGnBu", breaks = c(10, 100, 1000, 10000, 50000), direction = 1) +
    guides(fill=guide_legend(), size = guide_legend()) +
    coord_sf(crs = "+proj=aea +lat_1=25 +lat_2=50 +lon_0=-100") +
    theme_void()

ggsave("../plots/testing_coverage_bubble.svg", w= 10, h = 7.5)

testing_us_map %>%
    filter(!(state_id %in% c("GU", "PR", "HI", "AK"))) %>%
    ggplot() +
    geom_sf(aes(fill=n)) +
    scale_fill_distiller(trans="log", palette="YlGnBu", breaks = c(10, 100, 1000, 10000, 100000), direction = 1) +
    coord_sf(crs = "+proj=aea +lat_1=25 +lat_2=50 +lon_0=-100") +
    theme_void()

ggsave("../plots/testing_coverage.svg", w= 10, h = 7.5)
