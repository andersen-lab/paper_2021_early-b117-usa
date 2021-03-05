library(tidyverse)
library(lubridate)
library(zoo)
library(gridExtra)
library(sf)
library(nlstools)

## Combine sequencing and testing data
helix_metadata <- read_tsv("../data/covid_baseline_for_b117_paper.20210222.txt") %>%
    group_by(state) %>%
    group_modify(~{
        tmp <- .x %>%
            filter(collection_date <= max(collection_date) - 1)
        if(.y != "CA"){
            tmp %>%
                filter(collection_date <= as.Date("2021-02-11"))
        } else {
            tmp
        }
    })

state_codes <- read_tsv("./state_mapping.tsv")
seq_metadata <- read_csv("./us_b117_seq_metadata.csv") %>%
    filter(str_detect(accession, "STM")) %>% #Filter only Helix
    left_join(state_codes, by="state") %>%
    group_by(date, code) %>%
    count() %>%
    rename(
        seq_count = n
    )

seq_counts <- seq_metadata %>% inner_join(helix_metadata, by=c("code"= "state", "date" = "collection_date"))

usa_seq_counts <- seq_counts %>%
    group_by(date) %>%
    summarise(seq_count = sum(seq_count), n_sgtf_seq = sum(n_sgtf_seq)) %>%
    mutate(code = "USA") %>%
    ungroup()

seq_counts <- seq_counts %>%
    bind_rows(usa_seq_counts)

## Sequencing counts

tmp <- helix_metadata %>%
    group_by(collection_date) %>%
    summarise(n_sgtf = sum(n_sgtf), n = sum(n), n_b117 = sum(n_b117), n_sgtf_seq = sum(n_sgtf_seq)) %>%
    mutate(state = "USA")

tmp <- bind_rows(tmp, helix_metadata)

p <- tmp %>%
    filter(state %in% c("USA", "CA", "FL", "GA", "TX")) %>%
    mutate(
        prop_sequenced = ifelse(n_sgtf_seq >= n_b117, n_b117/n_sgtf_seq, 1) #Accounts for 1 day with more B117 sequences due to CQ < 27,
    ) %>%
    group_by(state) %>%
    group_map(~{
        max_seq_count <- .x %>% summarise(n = max(n_sgtf_seq)) %>% first()
        tmp <- .x
        .x %>%
            mutate(
                n_sgtf_seq = ifelse(n_sgtf_seq > n_b117, n_sgtf_seq - n_b117, 0)
            ) %>%
            select(n_b117, n_sgtf_seq, collection_date) %>%
            gather(key, val, -collection_date) %>%
            ggplot() +
            geom_col(aes(collection_date, val, fill = key)) +
            ## geom_line(data = tmp, aes(date, prop_sequenced * 50), fill="#000000") +
            ## geom_point(data = tmp, aes(date, prop_sequenced * 50)) +
            scale_y_reverse(limits=c(50, 0)) + theme_bw() + scale_x_date(limits = c(as.Date("2020-12-15"), as.Date("2021-03-01"))) + ggtitle(.y)
    })

svg("../plots/b117_num_seqs.svg", w = 20, h =10)
do.call("grid.arrange", c(p, ncol=3))
dev.off()

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
            scale_x_date(limits=c(as.Date("2020-12-15"), as.Date("2021-03-01")))  + scale_y_continuous(trans="pseudo_log", limits=c(0, 100), breaks = c(0, 25, 50, 75, 100))
    })
pdf("../plots/sgtf_drop_all_states_log.pdf", w = 15, h =12)
do.call("grid.arrange", c(p, ncol=5))
dev.off()

## Sampling
tmp <- helix_metadata %>%
    group_by(collection_date) %>%
    summarise(n_sgtf = sum(n_sgtf), n = sum(n), n_b117 = sum(), n_sgtf_seq = sum(n_sgtf_seq)) %>%
    mutate(state = "USA")

tmp <- bind_rows(tmp, helix_metadata)

states_with_sampling <- tmp %>% group_by(state) %>% summarise(n = sum(n_sgtf)) %>% filter( n > 0) %>% select(state) %>% as_vector()
## Sequence sampling
p <- tmp %>%
    filter(state %in% states_with_sampling) %>%
    mutate(
        sampling_50_pct = (n_sgtf * 0.5),
        sampling_40_pct = (n_sgtf * 0.4),
        sampling_30_pct = (n_sgtf * 0.3),
        sampling_20_pct = (n_sgtf * 0.2),
        sampling_10_pct = (n_sgtf * 0.1)
    ) %>%
    group_by(state) %>%
    group_map(~{
        .x %>%
            select(n_sgtf, n_sgtf_seq, sampling_10_pct, sampling_20_pct, sampling_30_pct, sampling_40_pct, sampling_50_pct) %>%
            gather(key, val, -n_sgtf, -n_sgtf_seq) %>%
            ggplot(aes(n_sgtf, n_sgtf_seq)) + geom_point() + theme_bw() + ggtitle(.y) + geom_line(aes(n_sgtf, val, color = key), linetype="dashed") + geom_smooth(method="lm", se=F, color="#1c86ee", size = 0.5) +
            scale_color_manual(labels = c("10% sampling", "20% sampling", "30% sampling", "40% sampling", "50% sampling"), values=c("indianred", "#ee9a49", "#d15fee", "#228b22", "#ff1493")) +
            xlab("SGTF dropout samples") + ylab("Sequenced SGTF samples") +
            theme(text = element_text(size = 20))
    })

pdf("../plots/sgtf_seq_sampling.pdf", w = 25, h =27.5)
do.call("grid.arrange", c(p, ncol=4))
dev.off()

## B117 out of SGTF
p <- tmp %>%
    filter(state %in% states_with_sampling) %>%
    mutate(
        sampling_80_pct = (n_sgtf_seq * 0.8),
        sampling_70_pct = (n_sgtf_seq * 0.7),
        sampling_60_pct = (n_sgtf_seq * 0.6),
        sampling_50_pct = (n_sgtf_seq * 0.5),
        sampling_40_pct = (n_sgtf_seq * 0.4)
    ) %>%
    group_by(state) %>%
    group_map(~{
        .x %>%
            select(n_b117, n_sgtf_seq, sampling_40_pct, sampling_50_pct, sampling_60_pct, sampling_70_pct, sampling_80_pct) %>%
            gather(key, val, -n_b117, -n_sgtf_seq) %>%
            ggplot(aes(n_sgtf_seq, n_b117)) + geom_point() + theme_bw() + ggtitle(.y) + geom_line(aes(n_sgtf_seq, val, color = key), linetype="dashed") + geom_smooth(method="lm", se=F, color="#1c86ee", size = 0.5) +
            scale_color_manual(labels = c("40% proportion", "50% proportion", "60% proportion", "70% proportion", "80% proportion"), values=c("indianred", "#ee9a49", "#d15fee", "#228b22", "#ff1493")) +
            xlab("SGTF dropout sequences") + ylab("B117 sequences") +
            theme(text = element_text(size = 20))
    })

pdf("../plots/sgtf_b117_proportion.pdf", w = 25, h =27.5)
do.call("grid.arrange", c(p, ncol=4))
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
            ggplot(aes(collection_date, b117_sgtf_pct)) + geom_col(fill="indianred") + theme_bw() + ylab("% B.1.1.7 in sequenced SGTF") + xlab("Date") + ggtitle(.y) +             scale_x_date(limits=c(as.Date("2020-12-15"), as.Date("2021-03-01"))) + scale_y_continuous(limits=c(0,100))
    })
pdf("../plots/b117_sgtf_pct.pdf", w = 15, h =12)
do.call("grid.arrange", c(p, ncol=5))
dev.off()

## Calculate rolling means for states
helix_metadata_rolling_avg <- helix_metadata %>%
    filter(state %in% c("CA", "FL", "GA", "PA", "TX")) %>%
    group_by(state) %>%
    group_modify(~{
        tmp <- .x
        last_date <- .x %>% filter(n_b117 > 0) %>% tail(n = 1) %>% select(collection_date) %>% first()
        first_date <- .x %>% filter(n_b117 > 0) %>% head(n = 1) %>% select(collection_date) %>% first()
        avg_last_5day_b117_seq_pct = tmp %>% filter(n_b117 > 0 & n_sgtf_seq > 0) %>% filter(collection_date >= last_date - 5) %>% summarise(b117_sgtf_pct = mean(n_b117/n_sgtf_seq)) %>% select(b117_sgtf_pct) %>% last() #Averge of last 5 days of B117
        print(paste0(.y, " ",avg_last_5day_b117_seq_pct))
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

## Calculate rolling means for USA
states_sampling <- helix_metadata_rolling_avg %>%
    group_by(state) %>%
    summarise(n = sum(n_b117, na.rm=T)) %>%
    filter(n >= 20) %>%
    select(state) %>%
    as_vector()

helix_metadata_us <- helix_metadata %>%
    ## filter(state %in% states_sampling) %>%
    filter(collection_date <= as.Date("2021-02-11")) %>%
    group_by(collection_date) %>%
    summarise(n_b117 = sum(n_b117), n_sgtf_seq = sum(n_sgtf_seq), n_sgtf = sum(n_sgtf), n = sum(n))

last_date <- helix_metadata_us %>% filter(n_b117 > 0) %>% tail(n = 1) %>% select(collection_date) %>% first()
first_date <- helix_metadata_us %>% filter(n_b117 > 0) %>% head(n = 1) %>% select(collection_date) %>% first()
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
    }) %>%
    mutate(
        state= "USA"
    )

## Merge USA into another "state"
helix_metadata_rolling_avg <- bind_rows(helix_metadata_us_rolling_avg, helix_metadata_rolling_avg)

## Plot rolling means
p <- helix_metadata_rolling_avg %>%
    filter(state %in% c("USA", "CA", "FL", "GA", "TX", "PA")) %>%
    group_by(state) %>%
    group_map(~{
        .x %>%
            ggplot() +
            geom_col(aes(collection_date, b117_est_pct_raw, fill = type)) +
            geom_line(aes(collection_date, b117_est_pct)) +
            geom_point(aes(collection_date, b117_est_pct)) +
            theme_bw() +
            scale_y_continuous(limits = c(0, 0.3), name="Proportion of B117 cases") + xlab("Days") +
            scale_x_date(limits=c(as.Date("2020-12-15"), as.Date("2021-03-01"))) + ggtitle(.y)
    })

svg("../plots/b117_rolling_mean.svg", w = 20, h =10)
do.call("grid.arrange", c(p, ncol=3))
dev.off()

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
    coef_fit <- coef(fit)
    ci_fit <- confint2(fit, method="asymptotic")
    tmp_newdata$predicted_lower = 1/(1+(((1/ci_fit[2,1]) - 1) * exp( -1 * ci_fit[1,1] * tmp_newdata$ndays)))
    tmp_newdata$predicted_upper = 1/(1+(((1/ci_fit[2,2]) - 1) * exp( -1 * ci_fit[1,2] * tmp_newdata$ndays)))

    above_50_date <- tmp_newdata %>% filter(predicted >= 0.5)  %>% select(collection_date) %>% first() %>% first()
    above_50_date_lower <- tmp_newdata %>% filter(predicted_lower >= 0.5)  %>% select(collection_date) %>% first() %>% first()
    above_50_date_upper <- tmp_newdata %>% filter(predicted_upper >= 0.5)  %>% select(collection_date) %>% first() %>% first()

    serial_interval <- 5.5

    print(paste(c(round(coef_fit[1]* 100, 4), paste(paste0(round(ci_fit[1,] * 100, 2), "%"), collapse="-"), round(coef_fit[2], 5), paste(round(ci_fit[2,], 5), collapse="-"), paste0(round(serial_interval * coef_fit[1] * 100, 4), " ", paste(paste0(round(serial_interval * ci_fit[1,] * 100, 1), "%"), collapse = " - ")), round(log(2)/coef(fit)[1], 2), paste(round(log(2)/ci_fit[1,], 2), collapse = " - ")), collapse=" | "))

    ggplot() +
        geom_rect(data = data.frame(xmin = above_50_date_upper,
                                    xmax = above_50_date_lower,
                                    ymin = -Inf,
                                    ymax = Inf),
                  aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                  fill = "#ff6a6a", alpha = 0.5) +
        geom_ribbon(data = tmp_newdata, aes(x = collection_date, ymin = predicted_lower, ymax = predicted_upper), alpha = 0.25, fill="#1e90ff") +
        geom_point(data = tmp, aes(collection_date, b117_est_pct, color = type)) +
        geom_line(data = tmp_newdata, aes(collection_date, predicted)) +
        geom_text(data = tmp, x = first_date, y = 1, label = paste0("Logistic growth rate = ", round(coef(fit)[1], 2), " per day"), hjust=0) +
        geom_text(data=tmp, x = first_date, y = 0.9, label = paste0("Increase in transmission = ", paste(paste0(round(serial_interval * ci_fit[1,] * 100, 1), "%"), collapse = " - ")), hjust=0) +
        ## geom_text(x = first_date, y = max(tmp$b117_est_pct, na.rm=T), label = paste0("Doubling time ~ 1 week "), hjust=0) +
        geom_text(data = tmp, x = first_date, y = 0.8, label = paste0("Rough doubling time = ", round(log(2)/coef(fit)[1], 2)), hjust=0) +
        geom_text(data = tmp, x = above_50_date, y = 0.1, label = paste0("50% date = ", above_50_date, " [",above_50_date_upper, " - ", above_50_date_lower,"]"), hjust=0) +
        geom_vline(xintercept = above_50_date, color="#ff6a6a") +
        theme_bw() + scale_x_date(limits=c(as.Date("2020-12-15"), as.Date("2021-02-20")), date_breaks="1 month", date_labels="%b %d") + scale_y_continuous(limits = c(0, 0.36), name="Proportion of B117 cases") + xlab("Days") + ggtitle(name)
}

p <- helix_metadata_rolling_avg %>%
    filter(state %in% c("USA", "CA", "FL", "GA", "TX")) %>%
    group_by(state) %>%
    group_map(~{
        print(.y)
        fit_plot_logistic_growth(.x, .y)
    })
## p1 <- fit_plot_logistic_growth(helix_metadata_us_rolling_avg, "USA")
## p2 <- helix_metadata_rolling_avg %>%
##     filter(state == "FL") %>%
##     fit_plot_logistic_growth("FL")
## p4 <- helix_metadata_rolling_avg %>%
##     filter(state == "CA") %>%
##     fit_plot_logistic_growth("CA")
svg("../plots/fit_b117_rolling_mean_inset.svg", w = 25, h =10)
do.call("grid.arrange", c(p, ncol=3))
dev.off()

## Get proportions for text
## USA
helix_metadata_rolling_avg %>%
    filter(!is.na(n_b117) & state == "USA" & collection_date >= as.Date("2021-01-29")) %>%
    mutate(
        prop_b117 = n_b117/n_sgtf_seq
    ) %>%
    ungroup() %>%
    summarise(prop_b117 = mean(prop_b117))

## SGTF dropout across USA
helix_metadata %>%
    ungroup() %>%
    filter(collection_date >= as.Date("2021-01-01") & collection_date <= as.Date("2021-01-07")) %>%
    summarise(n = sum(n_sgtf)/sum(n) * 100)

helix_metadata %>%
    ungroup() %>%
    filter(collection_date >= as.Date("2021-02-14")) %>%
    summarise(n = sum(n_sgtf)/sum(n) * 100)

## SGTF % in October
helix_metadata %>%
    ungroup() %>%
    filter(n_sgtf > 0 & collection_date >= as.Date("2020-10-18") & collection_date <= as.Date("2020-10-24")) %>%
    summarise(prop_sgtf = sum(n_sgtf)/sum(n))

## SGTF % in January
helix_metadata %>%
    ungroup() %>%
    filter(n_sgtf > 0 & collection_date >= as.Date("2021-01-01") & collection_date <= as.Date("2021-01-31")) %>%
    mutate(
        week = isoweek(collection_date)
    ) %>%
    group_by(week) %>%
    summarise(prop_sgtf = sum(n_sgtf)/sum(n))

## CA
helix_metadata_rolling_avg %>%
    filter(!is.na(n_b117) & state == "CA" & collection_date >= as.Date("2021-01-28")) %>%
    mutate(
        prop_b117 = n_b117/n_sgtf_seq
    ) %>%
    ungroup() %>%
    summarise(prop_b117 = mean(prop_b117))

## Florida
helix_metadata_rolling_avg %>%
    filter(!is.na(n_b117) & state == "FL" & collection_date >= as.Date("2021-01-29")) %>%
    mutate(
        prop_b117 = n_b117/n_sgtf_seq
    ) %>%
    ungroup() %>%
    summarise(prop_b117 = mean(prop_b117))

## MA
helix_metadata %>%
    filter(!is.na(n_b117) & !is.na(n_sgtf_seq) & state == "MA" & collection_date >= as.Date("2021-01-27") & collection_date <= as.Date("2021-02-02")) %>%
    mutate(
        prop_b117 = ifelse(n_sgtf_seq >0, n_b117/n_sgtf_seq, 0)
    ) %>%
    ungroup() %>%
    summarise(prop_b117 = mean(prop_b117))

## B117 estimated proportion
helix_metadata_rolling_avg %>%
    filter(collection_date >= as.Date("2021-02-04")) %>%
    group_by(state) %>%
    summarise(b117_est_pct = mean(b117_est_pct))

## Testing coverage
testing_counts <- helix_metadata %>% group_by(state) %>% summarise(n = sum(n))
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

## Sensitivity analysis
sensitivity_df <- seq(0.5, 1, 0.1) %>%
    map(~{
        b117_seq_pct_threshold <- .x

        ## Calculate rolling means for states
        helix_metadata_rolling_avg <- helix_metadata %>%
            filter(state %in% c("CA", "FL", "GA", "PA", "TX")) %>%
            group_by(state) %>%
            group_modify(~{
                tmp <- .x
                last_date <- .x %>% filter(n_b117 > 0) %>% tail(n = 1) %>% select(collection_date) %>% first()
                first_date <- .x %>% filter(n_b117 > 0) %>% head(n = 1) %>% select(collection_date) %>% first()
                avg_last_5day_b117_seq_pct = tmp %>% filter(n_b117 > 0 & n_sgtf_seq > 0) %>% filter(collection_date >= last_date - 5) %>% summarise(b117_sgtf_pct = mean(n_b117/n_sgtf_seq)) %>% select(b117_sgtf_pct) %>% last() #Averge of last 5 days of B117
                avg_last_5day_b117_seq_pct <- b117_seq_pct_threshold
                print(paste0(avg_last_5day_b117_seq_pct, " ", .y))
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

        ## Calculate rolling means for USA
        helix_metadata_us <- helix_metadata %>%
            group_by(collection_date) %>%
            summarise(n_b117 = sum(n_b117), n_sgtf_seq = sum(n_sgtf_seq), n_sgtf = sum(n_sgtf), n = sum(n))

        last_date <- helix_metadata_us %>% filter(n_b117 > 0) %>% tail(n = 1) %>% select(collection_date) %>% first()
        first_date <- helix_metadata_us %>% filter(n_b117 > 0) %>% head(n = 1) %>% select(collection_date) %>% first()
        avg_last_5day_b117_seq_pct = helix_metadata_us %>% filter(n_b117 > 0 & n_sgtf_seq > 0) %>% filter(collection_date >= last_date - 5) %>% summarise(b117_sgtf_pct = mean(n_b117/n_sgtf_seq)) %>% select(b117_sgtf_pct) %>% last() #Averge of last 5 days of B117
        avg_last_5day_b117_seq_pct <- b117_seq_pct_threshold

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
            }) %>%
            mutate(
                state= "USA"
            )

        ## Merge USA into another "state"
        helix_metadata_rolling_avg <- bind_rows(helix_metadata_us_rolling_avg, helix_metadata_rolling_avg)

        ## Plot rolling means
        p <- helix_metadata_rolling_avg %>%
            filter(state %in% c("USA", "CA", "FL", "GA", "TX", "PA")) %>%
            group_by(state) %>%
            group_map(~{
                .x %>%
                    ggplot() +
                    geom_col(aes(collection_date, b117_est_pct_raw, fill = type)) +
                    geom_line(aes(collection_date, b117_est_pct)) +
                    geom_point(aes(collection_date, b117_est_pct)) +
                    theme_bw() +
                    scale_y_continuous(limits = c(0, 0.16), name="Proportion of B117 cases") + xlab("Days") +
                    scale_x_date(limits=c(as.Date("2020-12-15"), as.Date("2021-03-01"))) + ggtitle(.y)
            })

        pdf(paste0("../plots/b117_rolling_mean_",b117_seq_pct_threshold,".pdf"), w = 20, h =10)
        do.call("grid.arrange", c(p, ncol=3))
        dev.off()

        ## Fit
        fit_logistic_growth <- function(df, name){
            tmp <- df %>%
                arrange(collection_date)
            first_date <- tmp %>%
                filter(b117_est_pct > 0) %>%
                arrange(collection_date) %>%
                ungroup()  %>%
                select(collection_date) %>%
                head(n = 1) %>%
                first()
            tmp <- tmp %>%
                filter(collection_date >= first_date)
            tmp$ndays = (tmp$collection_date - tmp$collection_date[1]) %>% as.integer()
            fit <- nls(b117_est_pct ~ 1/(1+(((1/x0) - 1) * exp( -1 * r * ndays))), data = tmp, start=list(r = 0.08, x0 = 0.001))
            tmp_newdata <- data.frame(ndays = seq(0,300), collection_date = first_date + seq(0,300))
            tmp_newdata$predicted = predict(fit, newdata = tmp_newdata)
            coef_fit <- coef(fit)
            ci_fit <- confint2(fit, method="asymptotic")
            tmp_newdata$predicted_lower = 1/(1+(((1/ci_fit[2,1]) - 1) * exp( -1 * ci_fit[1,1] * tmp_newdata$ndays)))
            tmp_newdata$predicted_upper = 1/(1+(((1/ci_fit[2,2]) - 1) * exp( -1 * ci_fit[1,2] * tmp_newdata$ndays)))

            above_50_date <- tmp_newdata %>% filter(predicted >= 0.5)  %>% select(collection_date) %>% first() %>% first()
            above_50_date_lower <- tmp_newdata %>% filter(predicted_lower >= 0.5)  %>% select(collection_date) %>% first() %>% first()
            above_50_date_upper <- tmp_newdata %>% filter(predicted_upper >= 0.5)  %>% select(collection_date) %>% first() %>% first()

            serial_interval <- 5.5

            data.frame(
                r = coef_fit[1],
                r_lower = ci_fit[1,1],
                r_upper= ci_fit[1,2],
                x0 = coef_fit[2],
                x0_lower = ci_fit[2,1],
                x0_upper = ci_fit[2,2],
                trans_increase = serial_interval * coef_fit[1] * 100,
                trans_increase_lower = serial_interval * ci_fit[1,1] * 100,
                trans_increase_upper = serial_interval * ci_fit[1,2] * 100,
                b117_seq_pct = b117_seq_pct_threshold,
                location = name
            )
        }


        helix_metadata_rolling_avg %>%
            filter(state %in% c("USA", "CA", "FL", "GA", "TX")) %>%
            group_by(state) %>%
            group_map(~{
                print(.y)
                fit_logistic_growth(.x, .y)
            })
    }) %>%
    bind_rows()

sensitivity_df %>%
    mutate(
        growth_rate = paste0(round(r*100, 2), "%", " [", paste0(round(r_lower * 100, 2)), "%-", paste0(round(r_upper * 100, 2)) ,"%]")
    ) %>%
    select(state, growth_rate, b117_seq_pct) %>%
    spread(key = state, val = growth_rate) %>%
    write_csv("../sentivitiy_analysis.csv")

## Read growth rate from csv with values from above
final_growth_rates <- read_csv("../growth_rate.csv") %>%
    mutate(
        b117_seq_pct = paste0("mean last 5 days (", round(b117_seq_pct,2), ")")
    )

p <- sensitivity_df %>%
    mutate(
        r = r * 100,
        r_lower = r_lower * 100,
        r_upper = r_upper * 100
    ) %>%
    mutate(
        b117_seq_pct = as.character(b117_seq_pct)
    ) %>%
    bind_rows(
        final_growth_rates
    ) %>%
    filter(state %in% c("CA", "FL", "GA", "USA")) %>%
    group_by(state) %>%
    group_map(~{
        .x %>%
            ggplot(aes(b117_seq_pct, r)) +
            geom_errorbar(aes(ymin = r_lower, ymax = r_upper), width = 0.05) +
            geom_point() +
            ggtitle(.y) +
            theme_bw() +
            xlab("(# of B117 sequences)/(# of sequenced SGTF samples) used to infer (# of B117)/(# of total positive tests)") +
            ylab("% Growth rate")
    })

pdf("../plots/sensitivity.pdf", w = 7.5, h =10)
do.call("grid.arrange", c(p, ncol=1))
dev.off()
