### Comprehensive evaluation of streamflow performance (17/06/2022)
## between OL and EnKF, covering:
## 1) long-term metrics
## 2) MAE disaggregated by different flow levels and different months 
## 3) flow duration curve: long-term distribuiton/hydrologcial signature

catch_key_feat_df = read.csv("June6th_Catch_KeyFeatures.csv",header = T)

catch_key_feat_df$ID  = as.character(catch_key_feat_df$ID)

###############

catnames_current = c("110003", "121002", "130319", "137201", 
                     "142001", "143110", "206014","206025",
                     "226222", "237200", "238235", "401210", 
                     "401212", "401217", "401230", "402206",
                     "403217", "403222", "403232", "403244",
                     "405218", "405227", "405264", "408202", 
                     "416305", "419053", "422319", "422338", "497")

OL_perf = data.frame(ID = catnames_current, rbind(Perf_110003_1518, Perf_121002_1518, Perf_130319_1518,
                                                  Perf_137201_1518, Perf_142001_1518, Perf_143110_1518, 
                                                  Perf_206014_1518, Perf_206025_1518, Perf_226222_1518,
                                                  Perf_237200_1518, Perf_238235_1518, Perf_401210_1518, 
                                                  Perf_401212_1518, Perf_401217_1518, Perf_401230_1518, 
                                                  Perf_402206_1518, Perf_403217_1518, Perf_403222_1518, 
                                                  Perf_403232_1518, Perf_403244_1518, Perf_405218_1518, 
                                                  Perf_405227_1518, Perf_405264_1518, Perf_408202_1518, 
                                                  Perf_416305_1518, Perf_419053_1518, Perf_422319_1518, 
                                                  Perf_422338_1518, Perf_497_1518))%>%dplyr::select(-F1)

EnKF_perf = data.frame(ID = catnames_current, rbind(Perf_110003_EnKF, Perf_121002_EnKF, Perf_130319_EnKF,
                                                    Perf_137201_EnKF, Perf_142001_EnKF, Perf_143110_EnKF, 
                                                    Perf_206014_EnKF, Perf_206025_EnKF, Perf_226222_EnKF,
                                                    Perf_237200_EnKF, Perf_238235_EnKF, Perf_401210_EnKF, 
                                                    Perf_401212_EnKF, Perf_401217_EnKF, Perf_401230_EnKF, 
                                                    Perf_402206_EnKF, Perf_403217_EnKF, Perf_403222_EnKF, 
                                                    Perf_403232_EnKF, Perf_403244_EnKF, Perf_405218_EnKF, 
                                                    Perf_405227_EnKF, Perf_405264_EnKF, Perf_408202_EnKF, 
                                                    Perf_416305_EnKF, Perf_419053_EnKF, Perf_422319_EnKF, 
                                                    Perf_422338_EnKF, Perf_497_EnKF))%>%dplyr::select(-F1)


OL_MAE= data.frame(ID = catnames_current, rbind(MAE_OL_110003, MAE_OL_121002, MAE_OL_130319,
                                                MAE_OL_137201, MAE_OL_142001, MAE_OL_143110, 
                                                MAE_OL_206014, MAE_OL_206025, MAE_OL_226222,
                                                MAE_OL_237200, MAE_OL_238235, MAE_OL_401210, 
                                                MAE_OL_401212, MAE_OL_401217, MAE_OL_401230, 
                                                MAE_OL_402206, MAE_OL_403217, MAE_OL_403222, 
                                                MAE_OL_403232, MAE_OL_403244, MAE_OL_405218, 
                                                MAE_OL_405227, MAE_OL_405264, MAE_OL_408202, 
                                                MAE_OL_416305, MAE_OL_419053, MAE_OL_422319, 
                                                MAE_OL_422338, MAE_OL_497))


EnKF_MAE= data.frame(ID = catnames_current, rbind(MAE_EnKF_110003, MAE_EnKF_121002, MAE_EnKF_130319,
                                                  MAE_EnKF_137201, MAE_EnKF_142001, MAE_EnKF_143110, 
                                                  MAE_EnKF_206014, MAE_EnKF_206025, MAE_EnKF_226222,
                                                  MAE_EnKF_237200, MAE_EnKF_238235, MAE_EnKF_401210, 
                                                  MAE_EnKF_401212, MAE_EnKF_401217, MAE_EnKF_401230, 
                                                  MAE_EnKF_402206, MAE_EnKF_403217, MAE_EnKF_403222, 
                                                  MAE_EnKF_403232, MAE_EnKF_403244, MAE_EnKF_405218, 
                                                  MAE_EnKF_405227, MAE_EnKF_405264, MAE_EnKF_408202, 
                                                  MAE_EnKF_416305, MAE_EnKF_419053, MAE_EnKF_422319, 
                                                  MAE_EnKF_422338, MAE_EnKF_497))

OL_FDC= data.frame(ID = catnames_current, rbind(FDC_OL_110003, FDC_OL_121002, FDC_OL_130319,
                                                FDC_OL_137201, FDC_OL_142001, FDC_OL_143110, 
                                                FDC_OL_206014, FDC_OL_206025, FDC_OL_226222,
                                                FDC_OL_237200, FDC_OL_238235, FDC_OL_401210, 
                                                FDC_OL_401212, FDC_OL_401217, FDC_OL_401230, 
                                                FDC_OL_402206, FDC_OL_403217, FDC_OL_403222, 
                                                FDC_OL_403232, FDC_OL_403244, FDC_OL_405218, 
                                                FDC_OL_405227, FDC_OL_405264, FDC_OL_408202, 
                                                FDC_OL_416305, FDC_OL_419053, FDC_OL_422319, 
                                                FDC_OL_422338, FDC_OL_497))

EnKF_FDC= data.frame(ID = catnames_current, rbind(FDC_EnKF_110003, FDC_EnKF_121002, FDC_EnKF_130319,
                                                  FDC_EnKF_137201, FDC_EnKF_142001, FDC_EnKF_143110, 
                                                  FDC_EnKF_206014, FDC_EnKF_206025, FDC_EnKF_226222,
                                                  FDC_EnKF_237200, FDC_EnKF_238235, FDC_EnKF_401210, 
                                                  FDC_EnKF_401212, FDC_EnKF_401217, FDC_EnKF_401230, 
                                                  FDC_EnKF_402206, FDC_EnKF_403217, FDC_EnKF_403222, 
                                                  FDC_EnKF_403232, FDC_EnKF_403244, FDC_EnKF_405218, 
                                                  FDC_EnKF_405227, FDC_EnKF_405264, FDC_EnKF_408202, 
                                                  FDC_EnKF_416305, FDC_EnKF_419053, FDC_EnKF_422319, 
                                                  FDC_EnKF_422338, FDC_EnKF_497))



OL_MAE_month = data.frame(ID = catnames_current, rbind(monthMAE_OL_110003, monthMAE_OL_121002, monthMAE_OL_130319,
                                                       monthMAE_OL_137201, monthMAE_OL_142001, monthMAE_OL_143110, 
                                                       monthMAE_OL_206014, monthMAE_OL_206025, monthMAE_OL_226222,
                                                       monthMAE_OL_237200, monthMAE_OL_238235, monthMAE_OL_401210, 
                                                       monthMAE_OL_401212, monthMAE_OL_401217, monthMAE_OL_401230, 
                                                       monthMAE_OL_402206, monthMAE_OL_403217, monthMAE_OL_403222, 
                                                       monthMAE_OL_403232, monthMAE_OL_403244, monthMAE_OL_405218, 
                                                       monthMAE_OL_405227, monthMAE_OL_405264, monthMAE_OL_408202, 
                                                       monthMAE_OL_416305, monthMAE_OL_419053, monthMAE_OL_422319, 
                                                       monthMAE_OL_422338, monthMAE_OL_497))
names(OL_MAE_month)[-1] = names(EnKF_MAE_month)[-1] = month.abb
rownames(OL_MAE_month) = rownames(EnKF_MAE_month)  =  catnames_current
EnKF_MAE_month = data.frame(ID = catnames_current, rbind(monthMAE_EnKF_110003, monthMAE_EnKF_121002, monthMAE_EnKF_130319,
                                                         monthMAE_EnKF_137201, monthMAE_EnKF_142001, monthMAE_EnKF_143110, 
                                                         monthMAE_EnKF_206014, monthMAE_EnKF_206025, monthMAE_EnKF_226222,
                                                         monthMAE_EnKF_237200, monthMAE_EnKF_238235, monthMAE_EnKF_401210, 
                                                         monthMAE_EnKF_401212, monthMAE_EnKF_401217, monthMAE_EnKF_401230, 
                                                         monthMAE_EnKF_402206, monthMAE_EnKF_403217, monthMAE_EnKF_403222, 
                                                         monthMAE_EnKF_403232, monthMAE_EnKF_403244, monthMAE_EnKF_405218, 
                                                         monthMAE_EnKF_405227, monthMAE_EnKF_405264, monthMAE_EnKF_408202, 
                                                         monthMAE_EnKF_416305, monthMAE_EnKF_419053, monthMAE_EnKF_422319, 
                                                         monthMAE_EnKF_422338, monthMAE_EnKF_497))

###### OL vs EnKF (16th June) ###### 

# library(reshape2) for melt function
KGE_diff = EnKF_perf$KGE  - OL_perf$KGE
NSE_diff = EnKF_perf$NSE  - OL_perf$NSE
NSElog_diff = EnKF_perf$NSE_log  - OL_perf$NSE_log
NSEsqrt_diff = EnKF_perf$NSE_sqrt  - OL_perf$NSE_sqrt
NSEbox_diff = EnKF_perf$NSE_box  - OL_perf$NSE_box
Bias_diff = abs(OL_perf$Bias) - abs(EnKF_perf$Bias) 

boxplot(KGE_diff,NSE_diff,NSElog_diff,NSEsqrt_diff,NSEbox_diff,ylim = c(-1,1))
abline(h = 0 , col = 'grey', lwd = 1.8)

DA_Diff_df = data.frame(ID = EnKF_perf$ID, KGE = KGE_diff, Bias = Bias_diff, 
                        NSE = NSE_diff, NSElog = NSElog_diff, 
                        NSEsqrt = NSEsqrt_diff, NSEbox = NSEbox_diff)%>%
  left_join(catch_key_feat_df%>%select(ID,status,humidity,lag1_SMAP_Q,lag1_SsQsim), by = "ID")%>%
  mutate(humid_fac = factor(ifelse(humidity>=0.7,"Wet","Dry")))%>%
  mutate(status = factor(status ,levels = c("Low-impact-routing","High-impact-routing") ))


### 1) difference in key quantitiatve metrics #### 

## wide to long format: same catchment, but different variable names ####

DA_Diff_melt = melt(DA_Diff_df, id.vars = c("ID","status", "humidity","lag1_SMAP_Q","lag1_SsQsim","humid_fac"))  
quantile(DA_Diff_df$KGE,seq(0,1,.1))

ggplot(DA_Diff_melt,aes(variable,value))  + 
  geom_hline(yintercept = 0, linetype = 'dashed', size = 1, col = 'grey70') + 
  stat_boxplot(geom = 'errorbar', width = 0.2, coef = 1.5,position = position_dodge(0.75)) + 
  geom_boxplot() + theme_classic() + graph.theme.beta + 
  labs(x="Performance metrics", y = "Change", title = "Performance at all catchments") + 
  scale_y_continuous(limits = c(-2,1)) 


ggplot(DA_Diff_melt,aes(variable,value,fill = humid_fac))  + 
  geom_hline(yintercept = 0, linetype = 'dashed', size = 1, col = 'grey70') + 
  stat_boxplot(geom = 'errorbar', width = 0.2, coef = 1.5,position = position_dodge(0.75)) + 
  geom_boxplot() +
  theme_classic() + graph.theme.beta + labs(x="Performance metrics", y = "Change") + 
  guides(fill = guide_legend(title = "Humidity",title.hjust = 0.5)) + 
  scale_y_continuous(limits = c(-2,1)) + scale_fill_manual(values = c("royalblue2","firebrick2"))


ggplot(DA_Diff_melt,aes(variable,value,fill = status))  + 
  geom_hline(yintercept = 0, linetype = 'dashed', size = 1, col = 'grey70') + 
  stat_boxplot(geom = 'errorbar', width = 0.2, coef = 1.5,position = position_dodge(0.75)) + 
  geom_boxplot() +
  theme_classic() + graph.theme.beta + labs(x="Performance metrics", y = "Change") + 
  guides(fill = guide_legend(title = "Routing impact",title.hjust = 0.5)) + 
  scale_y_continuous(limits = c(-2,1)) + 
  scale_fill_manual(values = c("royalblue2","firebrick2"),labels = c("Low","High"))


##### 2. MAE at different flow levels ##### 

Rel_MAE_all = (OL_MAE$all - EnKF_MAE$all)/(OL_MAE$all)
Rel_MAE_low  = (OL_MAE$low - EnKF_MAE$low)/(OL_MAE$low)
Rel_MAE_mid  = (OL_MAE$mid - EnKF_MAE$mid)/(OL_MAE$mid)
Rel_MAE_high  = (OL_MAE$high - EnKF_MAE$high)/(OL_MAE$high)
Rel_MAE_top = (OL_MAE$top - EnKF_MAE$top)/(OL_MAE$top)
Rel_MAE_peak = (OL_MAE$peak - EnKF_MAE$peak)/(OL_MAE$peak)

Rel_MAE_df  = data.frame(ID = EnKF_perf$ID, All = Rel_MAE_all, 
                         Low = Rel_MAE_low, Mid = Rel_MAE_mid, 
                         High = Rel_MAE_high, Peak = Rel_MAE_peak )%>%
  left_join(catch_key_feat_df%>%select(ID,status,humidity,lag1_SMAP_Q,lag1_SsQsim), by = "ID")%>%
  mutate(humid_fac = factor(ifelse(humidity>=0.7,"Wet","Dry")))%>%
  mutate(status = factor(status ,levels = c("Low-impact-routing","High-impact-routing") ))


Rel_MAE_melt = melt(Rel_MAE_df, id.vars = c("ID","status", "humidity","lag1_SMAP_Q","lag1_SsQsim","humid_fac"))  
Rel_MAE_melt = Rel_MAE_melt%>%mutate(value = 100*value)

Rel_MAE_df%>%filter(Peak >0)%>%group_by(humid_fac)%>%summarise(num = n())

ggplot(Rel_MAE_melt,aes(variable,value))  + 
  geom_hline(yintercept = 0, linetype = 'dashed', size = 1, col = 'grey70') + 
  stat_boxplot(geom = 'errorbar', width = 0.2, coef = 1.5,position = position_dodge(0.75)) + 
  geom_boxplot() + theme_classic() + graph.theme.beta + 
  labs(x="Flow level", y = "Relative MAE (%)", title = "Performance at all catchments") 


ggplot(Rel_MAE_melt,aes(variable,value,fill = humid_fac))  + 
  geom_hline(yintercept = 0, linetype = 'dashed', size = 1, col = 'grey70') + 
  stat_boxplot(geom = 'errorbar', width = 0.2, coef = 1.5,position = position_dodge(0.75)) + 
  geom_boxplot() +
  theme_classic() + graph.theme.beta + labs(x="Flow level", y = "Relative MAE (%)") + 
  guides(fill = guide_legend(title = "Humidity",title.hjust = 0.5)) +  
  scale_fill_manual(values = c("royalblue2","firebrick2")) +
  scale_y_continuous(limits = c(-70,20))


### Performance in flow duration curve - hydrological signature #### 

### compare DA vs OL: positive values represent better model skill ####

Rel_BiasFMS = (abs(OL_FDC$BiasFMS) - abs(EnKF_FDC$BiasFMS))/(abs(OL_FDC$BiasFMS))
Rel_BiasFHV = (abs(OL_FDC$BiasFHV) - abs(EnKF_FDC$BiasFHV))/(abs(OL_FDC$BiasFHV))
Rel_BiasFLV = (abs(OL_FDC$BiasFLV) - abs(EnKF_FDC$BiasFLV))/(abs(OL_FDC$BiasFLV))
Rel_BiasFMM = (abs(OL_FDC$BiasFMM) - abs(EnKF_FDC$BiasFMM))/(abs(OL_FDC$BiasFMM))
Rel_BiasFMV = (abs(OL_FDC$BiasFMV) - abs(EnKF_FDC$BiasFMV))/(abs(OL_FDC$BiasFMV))

Rel_FDC_df%>%arrange(desc(MidSegment))%>%filter(ID == "401230")

Rel_FDC_df  = data.frame(ID = EnKF_perf$ID, MidSlope = Rel_BiasFMS, MidSegment = Rel_BiasFMV,
                         HighSegment= Rel_BiasFHV, LowSegment = Rel_BiasFLV, Median = Rel_BiasFMM)%>%
  left_join(catch_key_feat_df%>%select(ID,status,humidity,lag1_SMAP_Q,lag1_SsQsim), by = "ID")%>%
  mutate(humid_fac = factor(ifelse(humidity>=0.7,"Wet","Dry")))%>%
  mutate(status = factor(status ,levels = c("Low-impact-routing","High-impact-routing") ))


Rel_FDC_melt = melt(Rel_FDC_df, id.vars = c("ID","status", "humidity","lag1_SMAP_Q","lag1_SsQsim","humid_fac"))  
Rel_FDC_melt = Rel_FDC_melt%>%mutate(value = 100*value)

ggplot(Rel_FDC_melt,aes(variable,value))  + 
  geom_hline(yintercept = 0, linetype = 'dashed', size = 1, col = 'grey70') + 
  stat_boxplot(geom = 'errorbar', width = 0.2, coef = 1.5,position = position_dodge(0.75)) + 
  geom_boxplot() + theme_classic() + graph.theme.beta + 
  labs(x="FDC Segment", y = "Skill in FDC signature measures (%)", title = "Performance at all catchments")  +
  scale_y_continuous(limits = c(-50,100))


ggplot(Rel_FDC_melt,aes(variable,value,fill = humid_fac))  + 
  geom_hline(yintercept = 0, linetype = 'dashed', size = 1, col = 'grey70') + 
  stat_boxplot(geom = 'errorbar', width = 0.2, coef = 1.5,position = position_dodge(0.75)) + 
  geom_boxplot() +
  theme_classic() + graph.theme.beta + labs(x="FDC Segment", y = "Skill in FDC signature measures (%)")  +
  guides(fill = guide_legend(title = "Humidity",title.hjust = 0.5)) +  
  scale_fill_manual(values = c("royalblue2","firebrick2"))+
  scale_y_continuous(limits = c(-50,100))


### Performance disaggregated by month ###### 


Rel_MAE_Jan = (OL_MAE_month$Jan - EnKF_MAE_month$Jan)/(OL_MAE_month$Jan)
Rel_MAE_Feb = (OL_MAE_month$Feb - EnKF_MAE_month$Feb)/(OL_MAE_month$Feb)
Rel_MAE_Mar = (OL_MAE_month$Mar - EnKF_MAE_month$Mar)/(OL_MAE_month$Mar)
Rel_MAE_Apr = (OL_MAE_month$Apr - EnKF_MAE_month$Apr)/(OL_MAE_month$Apr)
Rel_MAE_May = (OL_MAE_month$May - EnKF_MAE_month$May)/(OL_MAE_month$May)
Rel_MAE_Jun = (OL_MAE_month$Jun - EnKF_MAE_month$Jun)/(OL_MAE_month$Jun)
Rel_MAE_Jul = (OL_MAE_month$Jul - EnKF_MAE_month$Jul)/(OL_MAE_month$Jul)
Rel_MAE_Aug = (OL_MAE_month$Aug - EnKF_MAE_month$Aug)/(OL_MAE_month$Aug)
Rel_MAE_Sep = (OL_MAE_month$Sep - EnKF_MAE_month$Sep)/(OL_MAE_month$Sep)
Rel_MAE_Oct = (OL_MAE_month$Oct - EnKF_MAE_month$Oct)/(OL_MAE_month$Oct)
Rel_MAE_Nov = (OL_MAE_month$Nov - EnKF_MAE_month$Nov)/(OL_MAE_month$Nov)
Rel_MAE_Dec = (OL_MAE_month$Dec - EnKF_MAE_month$Dec)/(OL_MAE_month$Dec)


Rel_MAE_month_df  = data.frame(ID = EnKF_perf$ID, Jan = Rel_MAE_Jan, Feb = Rel_MAE_Feb, 
                               Mar = Rel_MAE_Mar, Apr = Rel_MAE_Apr, May = Rel_MAE_May,
                               Jun = Rel_MAE_Jun, Jul = Rel_MAE_Jul, Aug = Rel_MAE_Aug,
                               Sep = Rel_MAE_Sep, Oct = Rel_MAE_Oct, Nov = Rel_MAE_Nov, 
                               Dec = Rel_MAE_Dec)%>%
  left_join(catch_key_feat_df%>%select(ID,status,humidity,lag1_SMAP_Q,lag1_SsQsim), by = "ID")%>%
  mutate(humid_fac = factor(ifelse(humidity>=0.7,"Wet","Dry")))%>%
  mutate(status = factor(status ,levels = c("Low-impact-routing","High-impact-routing") ))

Rel_MAE_month_melt = melt(Rel_MAE_month_df, id.vars = c("ID","status", "humidity","lag1_SMAP_Q","lag1_SsQsim","humid_fac"))
Rel_MAE_month_melt = Rel_MAE_month_melt%>%mutate(value = 100*value)

ggplot(Rel_MAE_month_melt,aes(variable,value))  + 
  geom_hline(yintercept = 0, linetype = 'dashed', size = 1, col = 'grey70') + 
  stat_boxplot(geom = 'errorbar', width = 0.2, coef = 1.5,position = position_dodge(0.75)) + 
  geom_boxplot() + theme_classic() + graph.theme.beta + 
  labs(x="Month", y = "Relative MAE skill (%)", title = "Performance at all catchments")  +
  scale_y_continuous(limits = c(-100,30))

ggplot(Rel_MAE_month_melt,aes(variable,value,fill = humid_fac))  + 
  geom_hline(yintercept = 0, linetype = 'dashed', size = 1, col = 'grey70') + 
  stat_boxplot(geom = 'errorbar', width = 0.2, coef = 1.5,position = position_dodge(0.75)) + 
  geom_boxplot() +
  theme_classic() + graph.theme.beta + labs(x="Month", y = "Relative MAE skill (%)") + 
  guides(fill = guide_legend(title = "Humidity",title.hjust = 0.5)) +  
  scale_fill_manual(values = c("royalblue2","firebrick2")) +
  scale_y_continuous(limits = c(-70,20))
