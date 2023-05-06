#####GAMS COURSE#####
#https://github.com/gavinsimpson/physalia-gam-course


# update any installed R packages
#update.packages(ask = FALSE, checkBuilt = TRUE)

# packages to install
pkgs <- c("mgcv",  "brms", "qgam", "gamm4", "tidyverse", "readxl",
          "rstan", "mgcViz", "DHARMa", "gratia")
vapply(pkgs, library, logical(1), character.only = TRUE, logical.return = TRUE,
       quietly = TRUE)


data = read_csv("~/Dropbox (UFL)/Laughinghouse_Lab/MANUSCRIPTS/Lake_Okeechobee_Cyano/Apr2023/For_GAM_genus_RARE.csv") %>% dplyr::select(-c("SampleID"))

#####Temp#####
model_Microcystis_Temp <- gam(Microcystis ~ s(Temp,k=6, bs = "cr") + s(Location,bs="re") + s(Outing,bs="re"),
                              data = data, method = "REML", family= nb(link = "log"))
AIC(model_Microcystis_Temp)
k.check(model_Microcystis_Temp)
draw(model_Microcystis_Temp, residuals = TRUE, select = "s(Temp)")
summary(model_Microcystis_Temp)
appraise(model_Microcystis_Temp, method = 'simulate')
variance_comp(model_Microcystis_Temp)
resids = DHARMa::simulateResiduals(fittedModel = model_Microcystis_Temp, plot = FALSE)
plot(resids) #This is the better way to deteremine if a good model
#Abundance increases with temperature, with maximal abundance at ~23 and ~30C
#Lower abundances around 27.5C, decreases after 30C

sm <- smooth_estimates(model_Microcystis_Temp) %>%
  add_confint()
sm
eg1 <- data %>%
  add_partial_residuals(model_Microcystis_Temp)


plt_Microcystis_temp <- sm %>%
  filter(smooth == "s(Temp)") %>%
  ggplot() +
  geom_rug(aes(x = Temp),
           data = eg1,
           sides = "b", length = grid::unit(0.02, "npc")) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, x = Temp),
              alpha = 0.2) +
  geom_point(aes(x = Temp, y = `s(Temp)`),
             data = eg1, cex = 1.5, colour = "steelblue3") +
  geom_line(aes(x = Temp, y = est), lwd = 1.2) +
  labs(y = "Partial effect", title = "Microcystis ~ s(Temp)", x = "Temperature") +
  geom_point(aes(x = Temp, y = `s(Temp)`,
                 colour = Season, shape = Limitation), # <-- map fac to colour aesthetic
             data = eg1, cex = 1.5)
plt_Microcystis_temp



model_Raphidiopsis_Temp <- gam(Raphidiopsis ~ s(Temp,k=40, bs = "cr") + s(Location,bs="re") + s(Outing,bs="re"),
                              data = data, method = "REML", family= nb(link = "log"))
AIC(model_Raphidiopsis_Temp)
k.check(model_Raphidiopsis_Temp)
draw(model_Raphidiopsis_Temp, residuals = TRUE, select = "s(Temp)")
summary(model_Raphidiopsis_Temp)
appraise(model_Raphidiopsis_Temp, method = 'simulate')
variance_comp(model_Raphidiopsis_Temp)
resids = DHARMa::simulateResiduals(fittedModel = model_Raphidiopsis_Temp, plot = FALSE)
plot(resids) #This is the better way to deteremine if a good model
#Abundance increases with temperature, with maximal abundance at ~24 
#dips around 31 and then increases again

sm <- smooth_estimates(model_Raphidiopsis_Temp) %>%
  add_confint()
sm
eg1 <- data %>%
  add_partial_residuals(model_Raphidiopsis_Temp)


plt_Raphidiopsis_temp <- sm %>%
  filter(smooth == "s(Temp)") %>%
  ggplot() +
  geom_rug(aes(x = Temp),
           data = eg1,
           sides = "b", length = grid::unit(0.02, "npc")) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, x = Temp),
              alpha = 0.2) +
  geom_point(aes(x = Temp, y = `s(Temp)`),
             data = eg1, cex = 1.5, colour = "steelblue3") +
  geom_line(aes(x = Temp, y = est), lwd = 1.2) +
  labs(y = "Partial effect", title = "Raphidiopsis ~ s(Temp)") +
  geom_point(aes(x = Temp, y = `s(Temp)`,
                 colour = Season, shape = Limitation), # <-- map fac to colour aesthetic
             data = eg1, cex = 1.5)

plt_Raphidiopsis_temp



model_Cuspidothrix_Temp <- gam(Cuspidothrix ~ s(Temp,k=40, bs = "cr") + s(Location,bs="re") + s(Outing,bs="re"),
                               data = data, method = "REML", family= nb(link = "log"))
AIC(model_Cuspidothrix_Temp)
k.check(model_Cuspidothrix_Temp)
draw(model_Cuspidothrix_Temp, residuals = TRUE, select = "s(Temp)")
summary(model_Cuspidothrix_Temp)
appraise(model_Cuspidothrix_Temp, method = 'simulate')
variance_comp(model_Cuspidothrix_Temp)
resids = DHARMa::simulateResiduals(fittedModel = model_Cuspidothrix_Temp, plot = FALSE)
plot(resids) #This is the better way to deteremine if a good model
#Abundance increases with temperature, with maximal abundance at ~24 
#dips around 31 and then increases again

sm <- smooth_estimates(model_Cuspidothrix_Temp) %>%
  add_confint()
sm
eg1 <- data %>%
  add_partial_residuals(model_Cuspidothrix_Temp)


plt_Cuspidothrix_temp <- sm %>%
  filter(smooth == "s(Temp)") %>%
  ggplot() +
  geom_rug(aes(x = Temp),
           data = eg1,
           sides = "b", length = grid::unit(0.02, "npc")) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, x = Temp),
              alpha = 0.2) +
  geom_point(aes(x = Temp, y = `s(Temp)`),
             data = eg1, cex = 1.5, colour = "steelblue3") +
  geom_line(aes(x = Temp, y = est), lwd = 1.2) +
  labs(y = "Partial effect", title = "Cuspidothrix ~ s(Temp)") +
  geom_point(aes(x = Temp, y = `s(Temp)`,
                 colour = Season, shape = Limitation), # <-- map fac to colour aesthetic
             data = eg1, cex = 1.5)

plt_Cuspidothrix_temp





#WHY WONT YOU WORK#
model_Dolichospermum_Temp <- gam(Dolichospermum ~ s(Temp, k = 40, bs = "cr") + s(Location,bs="re") + s(Outing,bs="re"),
                               data = data, method = "REML", family= nb(link = "log"))
AIC(model_Dolichospermum_Temp)
k.check(model_Dolichospermum_Temp)
draw(model_Dolichospermum_Temp, residuals = TRUE, select = "s(Temp)")
summary(model_Dolichospermum_Temp)
appraise(model_Dolichospermum_Temp, method = 'simulate')
variance_comp(model_Dolichospermum_Temp)
resids = DHARMa::simulateResiduals(fittedModel = model_Dolichospermum_Temp, plot = FALSE)
plot(resids) #This is the better way to deteremine if a good model
#Abundance increases with temperature, with maximal abundance at ~23 and ~30C
#Lower abundances around 27.5C, decreases after 30C

sm <- smooth_estimates(model_Dolichospermum_Temp) %>%
  add_confint()
sm
eg1 <- data %>%
  add_partial_residuals(model_Dolichospermum_Temp)


plt_Dolichospermum_temp <- sm %>%
  filter(smooth == "s(Temp)") %>%
  ggplot() +
  geom_rug(aes(x = Temp),
           data = eg1,
           sides = "b", length = grid::unit(0.02, "npc")) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, x = Temp),
              alpha = 0.2) +
  geom_point(aes(x = Temp, y = `s(Temp)`),
             data = eg1, cex = 1.5, colour = "steelblue3") +
  geom_line(aes(x = Temp, y = est), lwd = 1.2) +
  labs(y = "Partial effect", title = "Dolichospermum ~ s(Temp)") +
  geom_point(aes(x = Temp, y = `s(Temp)`,
                 colour = Season, shape = Limitation), # <-- map fac to colour aesthetic
             data = eg1, cex = 1.5)
plt_Dolichospermum_temp




model_Cyanobium_Temp <- gam(Cyanobium ~ s(Temp) + s(Location,bs="re") + s(Outing,bs="re"),
                                 data = data, method = "REML", family= nb(link = "log"))
AIC(model_Cyanobium_Temp)
k.check(model_Cyanobium_Temp)
draw(model_Cyanobium_Temp, residuals = TRUE, select = "s(Temp)")
summary(model_Cyanobium_Temp)
appraise(model_Cyanobium_Temp, method = 'simulate')
variance_comp(model_Cyanobium_Temp)
resids = DHARMa::simulateResiduals(fittedModel = model_Cyanobium_Temp, plot = FALSE)
plot(resids) #This is the better way to deteremine if a good model
#temperatuer has no effect on cyanoboum 

sm <- smooth_estimates(model_Cyanobium_Temp) %>%
  add_confint()
sm
eg1 <- data %>%
  add_partial_residuals(model_Cyanobium_Temp)


plt_Cyanobium_temp <- sm %>%
  filter(smooth == "s(Temp)") %>%
  ggplot() +
  geom_rug(aes(x = Temp),
           data = eg1,
           sides = "b", length = grid::unit(0.02, "npc")) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, x = Temp),
              alpha = 0.2) +
  geom_point(aes(x = Temp, y = `s(Temp)`),
             data = eg1, cex = 1.5, colour = "steelblue3") +
  geom_line(aes(x = Temp, y = est), lwd = 1.2) +
  labs(y = "Partial effect", title = "Cyanobium ~ s(Temp)") +
  geom_point(aes(x = Temp, y = `s(Temp)`,
                 colour = Season, shape = Limitation), # <-- map fac to colour aesthetic
             data = eg1, cex = 1.5)
plt_Cyanobium_temp

#this one sucks too :(
model_Vulcanococcus_Temp <- gam(Vulcanococcus ~ s(Temp, k = 50, bs = "cr") + s(Location,bs="re") + s(Outing,bs="re"),
                            data = data, method = "REML", family= nb(link = "log"))
AIC(model_Vulcanococcus_Temp)
k.check(model_Vulcanococcus_Temp)
draw(model_Vulcanococcus_Temp, residuals = TRUE, select = "s(Temp)")
summary(model_Vulcanococcus_Temp)
appraise(model_Vulcanococcus_Temp, method = 'simulate')
variance_comp(model_Vulcanococcus_Temp)
resids = DHARMa::simulateResiduals(fittedModel = model_Vulcanococcus_Temp, plot = FALSE)
plot(resids) #This is the better way to deteremine if a good model
#Abundance increases with temperature, with maximal abundance at ~23 and ~30C
#Lower abundances around 27.5C, decreases after 30C

sm <- smooth_estimates(model_Vulcanococcus_Temp) %>%
  add_confint()
sm
eg1 <- data %>%
  add_partial_residuals(model_Vulcanococcus_Temp)


plt_Vulcanococcus_temp <- sm %>%
  filter(smooth == "s(Temp)") %>%
  ggplot() +
  geom_rug(aes(x = Temp),
           data = eg1,
           sides = "b", length = grid::unit(0.02, "npc")) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, x = Temp),
              alpha = 0.2) +
  geom_point(aes(x = Temp, y = `s(Temp)`),
             data = eg1, cex = 1.5, colour = "steelblue3") +
  geom_line(aes(x = Temp, y = est), lwd = 1.2) +
  labs(y = "Partial effect", title = "Vulcanococcus ~ s(Temp)") +
  geom_point(aes(x = Temp, y = `s(Temp)`,
                 colour = Season, shape = Limitation), # <-- map fac to colour aesthetic
             data = eg1, cex = 1.5)
plt_Vulcanococcus_temp

#####DIN####
model_Microcystis_DIN <- gam(Microcystis ~ s(DIN,k=20, bs = "cr") + s(Location,bs="re") + s(Outing,bs="re"),
                              data = data, method = "REML", family= nb(link = "log"))
AIC(model_Microcystis_DIN)
k.check(model_Microcystis_DIN)
draw(model_Microcystis_DIN, residuals = TRUE, select = "s(DIN)")
summary(model_Microcystis_DIN)
appraise(model_Microcystis_DIN, method = 'simulate')
variance_comp(model_Microcystis_DIN)
resids = DHARMa::simulateResiduals(fittedModel = model_Microcystis_DIN, plot = FALSE)
plot(resids) #This is the better way to deteremine if a good model
#

sm <- smooth_estimates(model_Microcystis_DIN) %>%
  add_confint()
sm
eg1 <- data %>% filter(DIN >= 0) %>%
  add_partial_residuals(model_Microcystis_DIN)


plt_Microcystis_DIN <- sm %>%
  filter(smooth == "s(DIN)") %>%
  ggplot() +
  geom_rug(aes(x = DIN),
           data = eg1,
           sides = "b", length = grid::unit(0.02, "npc")) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, x = DIN),
              alpha = 0.2) +
  geom_point(aes(x = DIN, y = `s(DIN)`),
             data = eg1, cex = 1.5, colour = "steelblue3") +
  geom_line(aes(x = DIN, y = est), lwd = 1.2) +
  labs(y = "Partial effect", title = "Microcystis ~ s(DIN)") +
  geom_point(aes(x = DIN, y = `s(DIN)`,
                 colour = Season, shape = Limitation), # <-- map fac to colour aesthetic
             data = eg1, cex = 1.5)
plt_Microcystis_DIN



model_Raphidiopsis_DIN <- gam(Raphidiopsis ~ s(DIN,k=25, bs = "cr") + s(Location,bs="re") + s(Outing,bs="re"),
                               data = data, method = "REML", family= nb(link = "log"))
AIC(model_Raphidiopsis_DIN)
k.check(model_Raphidiopsis_DIN)
draw(model_Raphidiopsis_DIN, residuals = TRUE, select = "s(DIN)")
summary(model_Raphidiopsis_DIN)
appraise(model_Raphidiopsis_DIN, method = 'simulate')
variance_comp(model_Raphidiopsis_DIN)
resids = DHARMa::simulateResiduals(fittedModel = model_Raphidiopsis_DIN, plot = FALSE)
plot(resids) #This is the better way to deteremine if a good model
#Abundance decreases with increasing concentarions of DIN

sm <- smooth_estimates(model_Raphidiopsis_DIN) %>%
  add_confint()
sm
eg1 <- data %>% filter(DIN >= 0) %>%
  add_partial_residuals(model_Raphidiopsis_DIN)


plt_Raphidiopsis_DIN <- sm %>%
  filter(smooth == "s(DIN)") %>%
  ggplot() +
  geom_rug(aes(x = DIN),
           data = eg1,
           sides = "b", length = grid::unit(0.02, "npc")) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, x = DIN),
              alpha = 0.2) +
  geom_point(aes(x = DIN, y = `s(DIN)`),
             data = eg1, cex = 1.5, colour = "steelblue3") +
  geom_line(aes(x = DIN, y = est), lwd = 1.2) +
  labs(y = "Partial effect", title = "Raphidiopsis ~ s(DIN)") +
  geom_point(aes(x = DIN, y = `s(DIN)`,
                 colour = Season, shape = Limitation), # <-- map fac to colour aesthetic
             data = eg1, cex = 1.5)
plt_Raphidiopsis_DIN


model_Cuspidothrix_DIN <- gam(Cuspidothrix ~ s(DIN,k=25, bs = "cr") + s(Location,bs="re") + s(Outing,bs="re"),
                              data = data, method = "REML", family= nb(link = "log"))
AIC(model_Cuspidothrix_DIN)
k.check(model_Cuspidothrix_DIN)
draw(model_Cuspidothrix_DIN, residuals = TRUE, select = "s(DIN)")
summary(model_Cuspidothrix_DIN)
appraise(model_Cuspidothrix_DIN, method = 'simulate')
variance_comp(model_Cuspidothrix_DIN)
resids = DHARMa::simulateResiduals(fittedModel = model_Cuspidothrix_DIN, plot = FALSE)
plot(resids) #This is the better way to deteremine if a good model
#Abundance decreases with increasing concentarions of DIN

sm <- smooth_estimates(model_Cuspidothrix_DIN) %>%
  add_confint()
sm
eg1 <- data %>% filter(DIN >= 0) %>%
  add_partial_residuals(model_Cuspidothrix_DIN)


plt_Cuspidothrix_DIN <- sm %>%
  filter(smooth == "s(DIN)") %>%
  ggplot() +
  geom_rug(aes(x = DIN),
           data = eg1,
           sides = "b", length = grid::unit(0.02, "npc")) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, x = DIN),
              alpha = 0.2) +
  geom_point(aes(x = DIN, y = `s(DIN)`),
             data = eg1, cex = 1.5, colour = "steelblue3") +
  geom_line(aes(x = DIN, y = est), lwd = 1.2) +
  labs(y = "Partial effect", title = "Cuspidothrix ~ s(DIN)") +
  geom_point(aes(x = DIN, y = `s(DIN)`,
                 colour = Season, shape = Limitation), # <-- map fac to colour aesthetic
             data = eg1, cex = 1.5)
plt_Cuspidothrix_DIN


#WHY WONT YOU WORK#
model_Dolichospermum_DIN <- gam(Dolichospermum ~ s(DIN, k=20, bs = "cr") + s(Location,bs="re") + s(Outing,bs="re"),
                                 data = data, method = "REML", family= nb(link = "log"))
AIC(model_Dolichospermum_DIN)
k.check(model_Dolichospermum_DIN)
draw(model_Dolichospermum_DIN, residuals = TRUE, select = "s(DIN)")
summary(model_Dolichospermum_DIN)
appraise(model_Dolichospermum_DIN, method = 'simulate')
variance_comp(model_Dolichospermum_DIN)
resids = DHARMa::simulateResiduals(fittedModel = model_Dolichospermum_DIN, plot = FALSE)
plot(resids) #This is the better way to deteremine if a good model
#Abundance decreases with DIN, but seems to increased between 0.3-0.5

sm <- smooth_estimates(model_Dolichospermum_DIN) %>%
  add_confint()
sm
eg1 <- data %>% filter(DIN >= 0) %>%
  add_partial_residuals(model_Dolichospermum_DIN)


plt_Dolichospermum_DIN <- sm %>%
  filter(smooth == "s(DIN)") %>%
  ggplot() +
  geom_rug(aes(x = DIN),
           data = eg1,
           sides = "b", length = grid::unit(0.02, "npc")) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, x = DIN),
              alpha = 0.2) +
  geom_point(aes(x = DIN, y = `s(DIN)`),
             data = eg1, cex = 1.5, colour = "steelblue3") +
  geom_line(aes(x = DIN, y = est), lwd = 1.2) +
  labs(y = "Partial effect", title = "Dolichospermum ~ s(DIN)") +
  geom_point(aes(x = DIN, y = `s(DIN)`,
                 colour = Season, shape = Limitation), # <-- map fac to colour aesthetic
             data = eg1, cex = 1.5)
plt_Dolichospermum_DIN




model_Cyanobium_DIN <- gam(Cyanobium ~ s(DIN) + s(Location,bs="re") + s(Outing,bs="re"),
                            data = data, method = "REML", family= nb(link = "log"))
AIC(model_Cyanobium_DIN)
k.check(model_Cyanobium_DIN)
draw(model_Cyanobium_DIN, residuals = TRUE, select = "s(DIN)")
summary(model_Cyanobium_DIN)
appraise(model_Cyanobium_DIN, method = 'simulate')
variance_comp(model_Cyanobium_DIN)
resids = DHARMa::simulateResiduals(fittedModel = model_Cyanobium_DIN, plot = FALSE)
plot(resids) #This is the better way to deteremine if a good model
#DINeratuer has no effect on cyanoboum 

sm <- smooth_estimates(model_Cyanobium_DIN) %>%
  add_confint()
sm
eg1 <- data %>% filter(DIN >= 0) %>%
  add_partial_residuals(model_Cyanobium_DIN)


plt_Cyanobium_DIN <- sm %>%
  filter(smooth == "s(DIN)") %>%
  ggplot() +
  geom_rug(aes(x = DIN),
           data = eg1,
           sides = "b", length = grid::unit(0.02, "npc")) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, x = DIN),
              alpha = 0.2) +
  geom_point(aes(x = DIN, y = `s(DIN)`),
             data = eg1, cex = 1.5, colour = "steelblue3") +
  geom_line(aes(x = DIN, y = est), lwd = 1.2) +
  labs(y = "Partial effect", title = "Cyanobium ~ s(DIN)") +
  geom_point(aes(x = DIN, y = `s(DIN)`,
                 colour = Season, shape = Limitation), # <-- map fac to colour aesthetic
             data = eg1, cex = 1.5)
plt_Cyanobium_DIN

#this one sucks too :(
model_Vulcanococcus_DIN <- gam(Vulcanococcus ~ s(DIN, bs = "tp") + s(Location,bs="re") + s(Outing,bs="re"),
                                data = data, method = "REML", family= nb(link = "log"))
AIC(model_Vulcanococcus_DIN)
k.check(model_Vulcanococcus_DIN)
draw(model_Vulcanococcus_DIN, residuals = TRUE, select = "s(DIN)")
summary(model_Vulcanococcus_DIN)
appraise(model_Vulcanococcus_DIN, method = 'simulate')
variance_comp(model_Vulcanococcus_DIN)
resids = DHARMa::simulateResiduals(fittedModel = model_Vulcanococcus_DIN, plot = FALSE)
plot(resids) #This is the better way to deteremine if a good model
#

sm <- smooth_estimates(model_Vulcanococcus_DIN) %>%
  add_confint()
sm
eg1 <- data %>% filter(DIN >= 0) %>%
  add_partial_residuals(model_Vulcanococcus_DIN)


plt_Vulcanococcus_DIN <- sm %>%
  filter(smooth == "s(DIN)") %>%
  ggplot() +
  geom_rug(aes(x = DIN),
           data = eg1,
           sides = "b", length = grid::unit(0.02, "npc")) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, x = DIN),
              alpha = 0.2) +
  geom_point(aes(x = DIN, y = `s(DIN)`),
             data = eg1, cex = 1.5, colour = "steelblue3") +
  geom_line(aes(x = DIN, y = est), lwd = 1.2) +
  labs(y = "Partial effect", title = "Vulcanococcus ~ s(DIN)") +
  geom_point(aes(x = DIN, y = `s(DIN)`,
                 colour = Season, shape = Limitation), # <-- map fac to colour aesthetic
             data = eg1, cex = 1.5)
plt_Vulcanococcus_DIN


#####Photic_Depth#####

model_Microcystis_Photic_Depth <- gam(Microcystis ~ s(Photic_Depth) + s(Location,bs="re") + s(Outing,bs="re"),
                             data = data, method = "REML", family= nb(link = "log"))
AIC(model_Microcystis_Photic_Depth)
k.check(model_Microcystis_Photic_Depth)
draw(model_Microcystis_Photic_Depth, residuals = TRUE, select = "s(Photic_Depth)")
summary(model_Microcystis_Photic_Depth)
appraise(model_Microcystis_Photic_Depth, method = 'simulate')
variance_comp(model_Microcystis_Photic_Depth)
resids = DHARMa::simulateResiduals(fittedModel = model_Microcystis_Photic_Depth, plot = FALSE)
plot(resids) #This is the better way to deteremine if a good model
#Abundance increses with incread pd

sm <- smooth_estimates(model_Microcystis_Photic_Depth) %>%
  add_confint()
sm
eg1 <- data %>%
  add_partial_residuals(model_Microcystis_Photic_Depth)


plt_Microcystis_Photic_Depth <- sm %>%
  filter(smooth == "s(Photic_Depth)") %>%
  ggplot() +
  geom_rug(aes(x = Photic_Depth),
           data = eg1,
           sides = "b", length = grid::unit(0.02, "npc")) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, x = Photic_Depth),
              alpha = 0.2) +
  geom_point(aes(x = Photic_Depth, y = `s(Photic_Depth)`),
             data = eg1, cex = 1.5, colour = "steelblue3") +
  geom_line(aes(x = Photic_Depth, y = est), lwd = 1.2) +
  labs(y = "Partial effect", title = "Microcystis ~ s(Photic_Depth)") +
  geom_point(aes(x = Photic_Depth, y = `s(Photic_Depth)`,
                 colour = Season, shape = Limitation), # <-- map fac to colour aesthetic
             data = eg1, cex = 1.5)
plt_Microcystis_Photic_Depth



model_Raphidiopsis_Photic_Depth <- gam(Raphidiopsis ~ s(Photic_Depth,k=10, bs = "cr") + s(Location,bs="re") + s(Outing,bs="re"),
                              data = data, method = "REML", family= nb(link = "log"))
AIC(model_Raphidiopsis_Photic_Depth)
k.check(model_Raphidiopsis_Photic_Depth)
draw(model_Raphidiopsis_Photic_Depth, residuals = TRUE, select = "s(Photic_Depth)")
summary(model_Raphidiopsis_Photic_Depth)
appraise(model_Raphidiopsis_Photic_Depth, method = 'simulate')
variance_comp(model_Raphidiopsis_Photic_Depth)
resids = DHARMa::simulateResiduals(fittedModel = model_Raphidiopsis_Photic_Depth, plot = FALSE)
plot(resids) #This is the better way to deteremine if a good model
#Abundance decreases with increasing Photic_Depth, it likes low photic depth??

sm <- smooth_estimates(model_Raphidiopsis_Photic_Depth) %>%
  add_confint()
sm
eg1 <- data %>%
  add_partial_residuals(model_Raphidiopsis_Photic_Depth)


plt_Raphidiopsis_Photic_Depth <- sm %>%
  filter(smooth == "s(Photic_Depth)") %>%
  ggplot() +
  geom_rug(aes(x = Photic_Depth),
           data = eg1,
           sides = "b", length = grid::unit(0.02, "npc")) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, x = Photic_Depth),
              alpha = 0.2) +
  geom_point(aes(x = Photic_Depth, y = `s(Photic_Depth)`),
             data = eg1, cex = 1.5, colour = "steelblue3") +
  geom_line(aes(x = Photic_Depth, y = est), lwd = 1.2) +
  labs(y = "Partial effect", title = "Raphidiopsis ~ s(Photic_Depth)") +
  geom_point(aes(x = Photic_Depth, y = `s(Photic_Depth)`,
                 colour = Season, shape = Limitation), # <-- map fac to colour aesthetic
             data = eg1, cex = 1.5)
plt_Raphidiopsis_Photic_Depth


model_Cuspidothrix_Photic_Depth <- gam(Cuspidothrix ~ s(Photic_Depth,k=10, bs = "cr") + s(Location,bs="re") + s(Outing,bs="re"),
                                       data = data, method = "REML", family= nb(link = "log"))
AIC(model_Cuspidothrix_Photic_Depth)
k.check(model_Cuspidothrix_Photic_Depth)
draw(model_Cuspidothrix_Photic_Depth, residuals = TRUE, select = "s(Photic_Depth)")
summary(model_Cuspidothrix_Photic_Depth)
appraise(model_Cuspidothrix_Photic_Depth, method = 'simulate')
variance_comp(model_Cuspidothrix_Photic_Depth)
resids = DHARMa::simulateResiduals(fittedModel = model_Cuspidothrix_Photic_Depth, plot = FALSE)
plot(resids) #This is the better way to deteremine if a good model
#Abundance decreases with increasing Photic_Depth, it likes low photic depth??

sm <- smooth_estimates(model_Cuspidothrix_Photic_Depth) %>%
  add_confint()
sm
eg1 <- data %>%
  add_partial_residuals(model_Cuspidothrix_Photic_Depth)


plt_Cuspidothrix_Photic_Depth <- sm %>%
  filter(smooth == "s(Photic_Depth)") %>%
  ggplot() +
  geom_rug(aes(x = Photic_Depth),
           data = eg1,
           sides = "b", length = grid::unit(0.02, "npc")) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, x = Photic_Depth),
              alpha = 0.2) +
  geom_point(aes(x = Photic_Depth, y = `s(Photic_Depth)`),
             data = eg1, cex = 1.5, colour = "steelblue3") +
  geom_line(aes(x = Photic_Depth, y = est), lwd = 1.2) +
  labs(y = "Partial effect", title = "Cuspidothrix ~ s(Photic_Depth)") +
  geom_point(aes(x = Photic_Depth, y = `s(Photic_Depth)`,
                 colour = Season, shape = Limitation), # <-- map fac to colour aesthetic
             data = eg1, cex = 1.5)
plt_Cuspidothrix_Photic_Depth





#
model_Dolichospermum_Photic_Depth <- gam(Dolichospermum ~ s(Photic_Depth, k = 5, bs="cr") + s(Location,bs="re") + s(Outing,bs="re"),
                                data = data, method = "REML", family= nb(link = "log"))
AIC(model_Dolichospermum_Photic_Depth)
k.check(model_Dolichospermum_Photic_Depth)
draw(model_Dolichospermum_Photic_Depth, residuals = TRUE, select = "s(Photic_Depth)")
summary(model_Dolichospermum_Photic_Depth)
appraise(model_Dolichospermum_Photic_Depth, method = 'simulate')
variance_comp(model_Dolichospermum_Photic_Depth)
resids = DHARMa::simulateResiduals(fittedModel = model_Dolichospermum_Photic_Depth, plot = FALSE)
plot(resids) #This is the better way to deteremine if a good model
#how confusing lol

sm <- smooth_estimates(model_Dolichospermum_Photic_Depth) %>%
  add_confint()
sm
eg1 <- data %>%
  add_partial_residuals(model_Dolichospermum_Photic_Depth)


plt_Dolichospermum_Photic_Depth <- sm %>%
  filter(smooth == "s(Photic_Depth)") %>%
  ggplot() +
  geom_rug(aes(x = Photic_Depth),
           data = eg1,
           sides = "b", length = grid::unit(0.02, "npc")) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, x = Photic_Depth),
              alpha = 0.2) +
  geom_point(aes(x = Photic_Depth, y = `s(Photic_Depth)`),
             data = eg1, cex = 1.5, colour = "steelblue3") +
  geom_line(aes(x = Photic_Depth, y = est), lwd = 1.2) +
  labs(y = "Partial effect", title = "Dolichospermum ~ s(Photic_Depth)") +
  geom_point(aes(x = Photic_Depth, y = `s(Photic_Depth)`,
                 colour = Season, shape = Limitation), # <-- map fac to colour aesthetic
             data = eg1, cex = 1.5)
plt_Dolichospermum_Photic_Depth




model_Cyanobium_Photic_Depth <- gam(Cyanobium ~ s(Photic_Depth) + s(Location,bs="re") + s(Outing,bs="re"),
                           data = data, method = "REML", family= nb(link = "log"))
AIC(model_Cyanobium_Photic_Depth)
k.check(model_Cyanobium_Photic_Depth)
draw(model_Cyanobium_Photic_Depth, residuals = TRUE, select = "s(Photic_Depth)")
summary(model_Cyanobium_Photic_Depth)
appraise(model_Cyanobium_Photic_Depth, method = 'simulate')
variance_comp(model_Cyanobium_Photic_Depth)
resids = DHARMa::simulateResiduals(fittedModel = model_Cyanobium_Photic_Depth, plot = FALSE)
plot(resids) #This is the better way to deteremine if a good model
#Photic_Deptheratuer has no effect on cyanoboum 

sm <- smooth_estimates(model_Cyanobium_Photic_Depth) %>%
  add_confint()
sm
eg1 <- data %>%
  add_partial_residuals(model_Cyanobium_Photic_Depth)


plt_Cyanobium_Photic_Depth <- sm %>%
  filter(smooth == "s(Photic_Depth)") %>%
  ggplot() +
  geom_rug(aes(x = Photic_Depth),
           data = eg1,
           sides = "b", length = grid::unit(0.02, "npc")) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, x = Photic_Depth),
              alpha = 0.2) +
  geom_point(aes(x = Photic_Depth, y = `s(Photic_Depth)`),
             data = eg1, cex = 1.5, colour = "steelblue3") +
  geom_line(aes(x = Photic_Depth, y = est), lwd = 1.2) +
  labs(y = "Partial effect", title = "Cyanobium ~ s(Photic_Depth)") +
  geom_point(aes(x = Photic_Depth, y = `s(Photic_Depth)`,
                 colour = Season, shape = Limitation), # <-- map fac to colour aesthetic
             data = eg1, cex = 1.5)
plt_Cyanobium_Photic_Depth




#
model_Vulcanococcus_Photic_Depth <- gam(Vulcanococcus ~ s(Photic_Depth, bs = "cr") + s(Location,bs="re") + s(Outing,bs="re"),
                               data = data, method = "REML", family= nb(link = "log"))
AIC(model_Vulcanococcus_Photic_Depth)
k.check(model_Vulcanococcus_Photic_Depth)
draw(model_Vulcanococcus_Photic_Depth, residuals = TRUE, select = "s(Photic_Depth)")
summary(model_Vulcanococcus_Photic_Depth)
appraise(model_Vulcanococcus_Photic_Depth, method = 'simulate')
variance_comp(model_Vulcanococcus_Photic_Depth)
resids = DHARMa::simulateResiduals(fittedModel = model_Vulcanococcus_Photic_Depth, plot = FALSE)
plot(resids) #This is the better way to deteremine if a good model
#Vulcanococcus likes it clear

sm <- smooth_estimates(model_Vulcanococcus_Photic_Depth) %>%
  add_confint()
sm
eg1 <- data %>%
  add_partial_residuals(model_Vulcanococcus_Photic_Depth)


plt_Vulcanococcus_Photic_Depth <- sm %>%
  filter(smooth == "s(Photic_Depth)") %>%
  ggplot() +
  geom_rug(aes(x = Photic_Depth),
           data = eg1,
           sides = "b", length = grid::unit(0.02, "npc")) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, x = Photic_Depth),
              alpha = 0.2) +
  geom_point(aes(x = Photic_Depth, y = `s(Photic_Depth)`),
             data = eg1, cex = 1.5, colour = "steelblue3") +
  geom_line(aes(x = Photic_Depth, y = est), lwd = 1.2) +
  labs(y = "Partial effect", title = "Vulcanococcus ~ s(Photic_Depth)") +
  geom_point(aes(x = Photic_Depth, y = `s(Photic_Depth)`,
                 colour = Season, shape = Limitation), # <-- map fac to colour aesthetic
             data = eg1, cex = 1.5)
plt_Vulcanococcus_Photic_Depth





#####DIN_DIP #####


model_Microcystis_DIN_DIP <- gam(Microcystis ~ s(DIN_DIP,k=10, bs = "cr") + s(Location,bs="re") + s(Outing,bs="re"),
                                      data = data, method = "REML", family= nb(link = "log"))
AIC(model_Microcystis_DIN_DIP)
k.check(model_Microcystis_DIN_DIP)
draw(model_Microcystis_DIN_DIP, residuals = TRUE, select = "s(DIN_DIP)")
summary(model_Microcystis_DIN_DIP)
appraise(model_Microcystis_DIN_DIP, method = 'simulate')
variance_comp(model_Microcystis_DIN_DIP)
resids = DHARMa::simulateResiduals(fittedModel = model_Microcystis_DIN_DIP, plot = FALSE)
plot(resids) #This is the better way to deteremine if a good model
#Abundance increses with incread DINDIP

sm <- smooth_estimates(model_Microcystis_DIN_DIP) %>%
  add_confint()
sm
eg1 <- data %>% filter(DIN_DIP >= 0) %>%
  add_partial_residuals(model_Microcystis_DIN_DIP)

#this one wont plot
plt_Microcystis_DIN_DIP <- sm %>%
  filter(smooth == "s(DIN_DIP)") %>%
  ggplot() +
  geom_rug(aes(x = DIN_DIP),
           data = eg1,
           sides = "b", length = grid::unit(0.02, "npc")) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, x = DIN_DIP),
              alpha = 0.2) +
  geom_point(aes(x = DIN_DIP, y = `s(DIN_DIP)`),
             data = eg1, cex = 1.5, colour = "steelblue3") +
  geom_line(aes(x = DIN_DIP, y = est), lwd = 1.2) +
  labs(y = "Partial effect", title = "Microcystis ~ s(DIN_DIP)") +
  geom_point(aes(x = DIN_DIP, y = `s(DIN_DIP)`,
                 colour = Season, shape = Limitation), # <-- map fac to colour aesthetic
             data = eg1, cex = 1.5)
plt_Microcystis_DIN_DIP



model_Raphidiopsis_DIN_DIP <- gam(Raphidiopsis ~ s(DIN_DIP,k=6, bs = "cr") + s(Location,bs="re") + s(Outing,bs="re"),
                                       data = data, method = "REML", family= nb(link = "log"))
AIC(model_Raphidiopsis_DIN_DIP)
k.check(model_Raphidiopsis_DIN_DIP)
draw(model_Raphidiopsis_DIN_DIP, residuals = TRUE, select = "s(DIN_DIP)")
summary(model_Raphidiopsis_DIN_DIP)
appraise(model_Raphidiopsis_DIN_DIP, method = 'simulate')
variance_comp(model_Raphidiopsis_DIN_DIP)
resids = DHARMa::simulateResiduals(fittedModel = model_Raphidiopsis_DIN_DIP, plot = FALSE)
plot(resids) #This is the better way to deteremine if a good model
#Abundance decreases with increasing DIN_DIP, it likes low photic depth??

sm <- smooth_estimates(model_Raphidiopsis_DIN_DIP) %>%
  add_confint()
sm
eg1 <- data %>% filter(DIN_DIP >= 0) %>%
  add_partial_residuals(model_Raphidiopsis_DIN_DIP)


plt_Raphidiopsis_DIN_DIP <- sm %>%
  filter(smooth == "s(DIN_DIP)") %>%
  ggplot() +
  geom_rug(aes(x = DIN_DIP),
           data = eg1,
           sides = "b", length = grid::unit(0.02, "npc")) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, x = DIN_DIP),
              alpha = 0.2) +
  geom_point(aes(x = DIN_DIP, y = `s(DIN_DIP)`),
             data = eg1, cex = 1.5, colour = "steelblue3") +
  geom_line(aes(x = DIN_DIP, y = est), lwd = 1.2) +
  labs(y = "Partial effect", title = "Raphidiopsis ~ s(DIN_DIP)") +
  geom_point(aes(x = DIN_DIP, y = `s(DIN_DIP)`,
                 colour = Season, shape = Limitation), # <-- map fac to colour aesthetic
             data = eg1, cex = 1.5)
plt_Raphidiopsis_DIN_DIP



model_Cuspidothrix_DIN_DIP <- gam(Cuspidothrix ~ s(DIN_DIP,k=6, bs = "cr") + s(Location,bs="re") + s(Outing,bs="re"),
                                  data = data, method = "REML", family= nb(link = "log"))
AIC(model_Cuspidothrix_DIN_DIP)
k.check(model_Cuspidothrix_DIN_DIP)
draw(model_Cuspidothrix_DIN_DIP, residuals = TRUE, select = "s(DIN_DIP)")
summary(model_Cuspidothrix_DIN_DIP)
appraise(model_Cuspidothrix_DIN_DIP, method = 'simulate')
variance_comp(model_Cuspidothrix_DIN_DIP)
resids = DHARMa::simulateResiduals(fittedModel = model_Cuspidothrix_DIN_DIP, plot = FALSE)
plot(resids) #This is the better way to deteremine if a good model
#Abundance decreases with increasing DIN_DIP, it likes low photic depth??

sm <- smooth_estimates(model_Cuspidothrix_DIN_DIP) %>%
  add_confint()
sm
eg1 <- data %>% filter(DIN_DIP >= 0) %>%
  add_partial_residuals(model_Cuspidothrix_DIN_DIP)


plt_Cuspidothrix_DIN_DIP <- sm %>%
  filter(smooth == "s(DIN_DIP)") %>%
  ggplot() +
  geom_rug(aes(x = DIN_DIP),
           data = eg1,
           sides = "b", length = grid::unit(0.02, "npc")) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, x = DIN_DIP),
              alpha = 0.2) +
  geom_point(aes(x = DIN_DIP, y = `s(DIN_DIP)`),
             data = eg1, cex = 1.5, colour = "steelblue3") +
  geom_line(aes(x = DIN_DIP, y = est), lwd = 1.2) +
  labs(y = "Partial effect", title = "Cuspidothrix ~ s(DIN_DIP)") +
  geom_point(aes(x = DIN_DIP, y = `s(DIN_DIP)`,
                 colour = Season, shape = Limitation), # <-- map fac to colour aesthetic
             data = eg1, cex = 1.5)
plt_Cuspidothrix_DIN_DIP


#
model_Dolichospermum_DIN_DIP <- gam(Dolichospermum ~ s(DIN_DIP,k = 4, bs="ps") + s(Location,bs="re") + s(Outing,bs="re"),
                                         data = data, method = "REML", family= nb(link = "log"))
AIC(model_Dolichospermum_DIN_DIP)
k.check(model_Dolichospermum_DIN_DIP)
draw(model_Dolichospermum_DIN_DIP, residuals = TRUE, select = "s(DIN_DIP)")
summary(model_Dolichospermum_DIN_DIP)
appraise(model_Dolichospermum_DIN_DIP, method = 'simulate')
variance_comp(model_Dolichospermum_DIN_DIP)
resids = DHARMa::simulateResiduals(fittedModel = model_Dolichospermum_DIN_DIP, plot = FALSE)
plot(resids) #This is the better way to deteremine if a good model
#how confusing lol

sm <- smooth_estimates(model_Dolichospermum_DIN_DIP) %>%
  add_confint()
sm
eg1 <- data %>% filter(DIN_DIP >= 0) %>%
  add_partial_residuals(model_Dolichospermum_DIN_DIP)


plt_Dolichospermum_DIN_DIP <- sm %>%
  filter(smooth == "s(DIN_DIP)") %>%
  ggplot() +
  geom_rug(aes(x = DIN_DIP),
           data = eg1,
           sides = "b", length = grid::unit(0.02, "npc")) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, x = DIN_DIP),
              alpha = 0.2) +
  geom_point(aes(x = DIN_DIP, y = `s(DIN_DIP)`),
             data = eg1, cex = 1.5, colour = "steelblue3") +
  geom_line(aes(x = DIN_DIP, y = est), lwd = 1.2) +
  labs(y = "Partial effect", title = "Dolichospermum ~ s(DIN_DIP)") +
  geom_point(aes(x = DIN_DIP, y = `s(DIN_DIP)`,
                 colour = Season, shape = Limitation), # <-- map fac to colour aesthetic
             data = eg1, cex = 1.5)
plt_Dolichospermum_DIN_DIP




model_Cyanobium_DIN_DIP <- gam(Cyanobium ~ s(DIN_DIP, bs = "cr") + s(Location,bs="re") + s(Outing,bs="re"),
                                    data = data, method = "REML", family= nb(link = "log"))
AIC(model_Cyanobium_DIN_DIP)
k.check(model_Cyanobium_DIN_DIP)
draw(model_Cyanobium_DIN_DIP, residuals = TRUE, select = "s(DIN_DIP)")
summary(model_Cyanobium_DIN_DIP)
appraise(model_Cyanobium_DIN_DIP, method = 'simulate')
variance_comp(model_Cyanobium_DIN_DIP)
resids = DHARMa::simulateResiduals(fittedModel = model_Cyanobium_DIN_DIP, plot = FALSE)
plot(resids) #This is the better way to deteremine if a good model
#DIN_DIPeratuer has no effect on cyanoboum 

sm <- smooth_estimates(model_Cyanobium_DIN_DIP) %>%
  add_confint()
sm
eg1 <- data %>% filter(DIN_DIP >= 0) %>%
  add_partial_residuals(model_Cyanobium_DIN_DIP)


plt_Cyanobium_DIN_DIP <- sm %>%
  filter(smooth == "s(DIN_DIP)") %>%
  ggplot() +
  geom_rug(aes(x = DIN_DIP),
           data = eg1,
           sides = "b", length = grid::unit(0.02, "npc")) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, x = DIN_DIP),
              alpha = 0.2) +
  geom_point(aes(x = DIN_DIP, y = `s(DIN_DIP)`),
             data = eg1, cex = 1.5, colour = "steelblue3") +
  geom_line(aes(x = DIN_DIP, y = est), lwd = 1.2) +
  labs(y = "Partial effect", title = "Cyanobium ~ s(DIN_DIP)") +
  geom_point(aes(x = DIN_DIP, y = `s(DIN_DIP)`,
                 colour = Season, shape = Limitation), # <-- map fac to colour aesthetic
             data = eg1, cex = 1.5)
plt_Cyanobium_DIN_DIP




#
model_Vulcanococcus_DIN_DIP <- gam(Vulcanococcus ~ s(DIN_DIP, k = 4) + s(Location,bs="re") + s(Outing,bs="re"),
                                        data = data, method = "REML", family= nb(link = "log"))
AIC(model_Vulcanococcus_DIN_DIP)
k.check(model_Vulcanococcus_DIN_DIP)
draw(model_Vulcanococcus_DIN_DIP, residuals = TRUE, select = "s(DIN_DIP)")
summary(model_Vulcanococcus_DIN_DIP)
appraise(model_Vulcanococcus_DIN_DIP, method = 'simulate')
variance_comp(model_Vulcanococcus_DIN_DIP)
resids = DHARMa::simulateResiduals(fittedModel = model_Vulcanococcus_DIN_DIP, plot = FALSE)
plot(resids) #This is the better way to deteremine if a good model
#Vulcanococcus likes it clear

sm <- smooth_estimates(model_Vulcanococcus_DIN_DIP) %>%
  add_confint()
sm
eg1 <- data %>% filter(DIN_DIP >= 0) %>%
  add_partial_residuals(model_Vulcanococcus_DIN_DIP)


plt_Vulcanococcus_DIN_DIP <- sm %>%
  filter(smooth == "s(DIN_DIP)") %>%
  ggplot() +
  geom_rug(aes(x = DIN_DIP),
           data = eg1,
           sides = "b", length = grid::unit(0.02, "npc")) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, x = DIN_DIP),
              alpha = 0.2) +
  geom_point(aes(x = DIN_DIP, y = `s(DIN_DIP)`),
             data = eg1, cex = 1.5, colour = "steelblue3") +
  geom_line(aes(x = DIN_DIP, y = est), lwd = 1.2) +
  labs(y = "Partial effect", title = "Vulcanococcus ~ s(DIN_DIP)") +
  geom_point(aes(x = DIN_DIP, y = `s(DIN_DIP)`,
                 colour = Season, shape = Limitation), # <-- map fac to colour aesthetic
             data = eg1, cex = 1.5)
plt_Vulcanococcus_DIN_DIP

#####Conductivity_SPC#####

model_Microcystis_Conductivity_SPC <- gam(Microcystis ~ s(Conductivity_SPC, k =4, bs = "cr") + s(Location,bs="re") + s(Outing,bs="re"),
                                 data = data, method = "REML", family= nb(link = "log"))
AIC(model_Microcystis_Conductivity_SPC)
k.check(model_Microcystis_Conductivity_SPC)
draw(model_Microcystis_Conductivity_SPC, residuals = TRUE, select = "s(Conductivity_SPC)")
summary(model_Microcystis_Conductivity_SPC)
appraise(model_Microcystis_Conductivity_SPC, method = 'simulate')
variance_comp(model_Microcystis_Conductivity_SPC)
resids = DHARMa::simulateResiduals(fittedModel = model_Microcystis_Conductivity_SPC, plot = FALSE)
plot(resids) #This is the better way to deteremine if a good model
#Abundance increses with incread DINDIP

sm <- smooth_estimates(model_Microcystis_Conductivity_SPC) %>%
  add_confint()
sm
eg1 <- data %>% filter(Conductivity_SPC >= 0) %>%
  add_partial_residuals(model_Microcystis_Conductivity_SPC)

#this one wont plot
plt_Microcystis_Conductivity_SPC <- sm %>%
  filter(smooth == "s(Conductivity_SPC)") %>%
  ggplot() +
  geom_rug(aes(x = Conductivity_SPC),
           data = eg1,
           sides = "b", length = grid::unit(0.02, "npc")) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, x = Conductivity_SPC),
              alpha = 0.2) +
  geom_point(aes(x = Conductivity_SPC, y = `s(Conductivity_SPC)`),
             data = eg1, cex = 1.5, colour = "steelblue3") +
  geom_line(aes(x = Conductivity_SPC, y = est), lwd = 1.2) +
  labs(y = "Partial effect", title = "Microcystis ~ s(Conductivity_SPC)") +
  geom_point(aes(x = Conductivity_SPC, y = `s(Conductivity_SPC)`,
                 colour = Season, shape = Limitation), # <-- map fac to colour aesthetic
             data = eg1, cex = 1.5)
plt_Microcystis_Conductivity_SPC



model_Raphidiopsis_Conductivity_SPC <- gam(Raphidiopsis ~ s(Conductivity_SPC,k=3, bs = "tp") + s(Location,bs="re") + s(Outing,bs="re"),
                                  data = data, method = "REML", family= nb(link = "log"))
AIC(model_Raphidiopsis_Conductivity_SPC)
k.check(model_Raphidiopsis_Conductivity_SPC)
draw(model_Raphidiopsis_Conductivity_SPC, residuals = TRUE, select = "s(Conductivity_SPC)")
summary(model_Raphidiopsis_Conductivity_SPC)
appraise(model_Raphidiopsis_Conductivity_SPC, method = 'simulate')
variance_comp(model_Raphidiopsis_Conductivity_SPC)
resids = DHARMa::simulateResiduals(fittedModel = model_Raphidiopsis_Conductivity_SPC, plot = FALSE)
plot(resids) #This is the better way to deteremine if a good model
#Abundance decreases with increasing Conductivity_SPC, it likes low photic depth??

sm <- smooth_estimates(model_Raphidiopsis_Conductivity_SPC) %>%
  add_confint()
sm
eg1 <- data %>% filter(Conductivity_SPC >= 0) %>%
  add_partial_residuals(model_Raphidiopsis_Conductivity_SPC)


plt_Raphidiopsis_Conductivity_SPC <- sm %>%
  filter(smooth == "s(Conductivity_SPC)") %>%
  ggplot() +
  geom_rug(aes(x = Conductivity_SPC),
           data = eg1,
           sides = "b", length = grid::unit(0.02, "npc")) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, x = Conductivity_SPC),
              alpha = 0.2) +
  geom_point(aes(x = Conductivity_SPC, y = `s(Conductivity_SPC)`),
             data = eg1, cex = 1.5, colour = "steelblue3") +
  geom_line(aes(x = Conductivity_SPC, y = est), lwd = 1.2) +
  labs(y = "Partial effect", title = "Raphdiopsis ~ s(Conductivity_SPC)") +
  geom_point(aes(x = Conductivity_SPC, y = `s(Conductivity_SPC)`,
                 colour = Season, shape = Limitation), # <-- map fac to colour aesthetic
             data = eg1, cex = 1.5)
plt_Raphidiopsis_Conductivity_SPC


model_Cuspidothrix_Conductivity_SPC <- gam(Cuspidothrix ~ s(Conductivity_SPC,k=5, bs = "tp") + s(Location,bs="re") + s(Outing,bs="re"),
                                           data = data, method = "REML", family= nb(link = "log"))
AIC(model_Cuspidothrix_Conductivity_SPC)
k.check(model_Cuspidothrix_Conductivity_SPC)
draw(model_Cuspidothrix_Conductivity_SPC, residuals = TRUE, select = "s(Conductivity_SPC)")
summary(model_Cuspidothrix_Conductivity_SPC)
appraise(model_Cuspidothrix_Conductivity_SPC, method = 'simulate')
variance_comp(model_Cuspidothrix_Conductivity_SPC)
resids = DHARMa::simulateResiduals(fittedModel = model_Cuspidothrix_Conductivity_SPC, plot = FALSE)
plot(resids) #This is the better way to deteremine if a good model
#Abundance decreases with increasing Conductivity_SPC, it likes low photic depth??

sm <- smooth_estimates(model_Cuspidothrix_Conductivity_SPC) %>%
  add_confint()
sm
eg1 <- data %>% filter(Conductivity_SPC >= 0) %>%
  add_partial_residuals(model_Cuspidothrix_Conductivity_SPC)


plt_Cuspidothrix_Conductivity_SPC <- sm %>%
  filter(smooth == "s(Conductivity_SPC)") %>%
  ggplot() +
  geom_rug(aes(x = Conductivity_SPC),
           data = eg1,
           sides = "b", length = grid::unit(0.02, "npc")) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, x = Conductivity_SPC),
              alpha = 0.2) +
  geom_point(aes(x = Conductivity_SPC, y = `s(Conductivity_SPC)`),
             data = eg1, cex = 1.5, colour = "steelblue3") +
  geom_line(aes(x = Conductivity_SPC, y = est), lwd = 1.2) +
  labs(y = "Partial effect", title = "Cuspidothrix ~ s(Conductivity_SPC)") +
  geom_point(aes(x = Conductivity_SPC, y = `s(Conductivity_SPC)`,
                 colour = Season, shape = Limitation), # <-- map fac to colour aesthetic
             data = eg1, cex = 1.5)
plt_Cuspidothrix_Conductivity_SPC

#
model_Dolichospermum_Conductivity_SPC <- gam(Dolichospermum ~ s(Conductivity_SPC,k = 4, bs="tp") + s(Location,bs="re") + s(Outing,bs="re"),
                                    data = data, method = "REML", family= nb(link = "log"))
AIC(model_Dolichospermum_Conductivity_SPC)
k.check(model_Dolichospermum_Conductivity_SPC)
draw(model_Dolichospermum_Conductivity_SPC, residuals = TRUE, select = "s(Conductivity_SPC)")
summary(model_Dolichospermum_Conductivity_SPC)
appraise(model_Dolichospermum_Conductivity_SPC, method = 'simulate')
variance_comp(model_Dolichospermum_Conductivity_SPC)
resids = DHARMa::simulateResiduals(fittedModel = model_Dolichospermum_Conductivity_SPC, plot = FALSE)
plot(resids) #This is the better way to deteremine if a good model
#how confusing lol

sm <- smooth_estimates(model_Dolichospermum_Conductivity_SPC) %>%
  add_confint()
sm
eg1 <- data %>% filter(Conductivity_SPC >= 0) %>%
  add_partial_residuals(model_Dolichospermum_Conductivity_SPC)


plt_Dolichospermum_Conductivity_SPC <- sm %>%
  filter(smooth == "s(Conductivity_SPC)") %>%
  ggplot() +
  geom_rug(aes(x = Conductivity_SPC),
           data = eg1,
           sides = "b", length = grid::unit(0.02, "npc")) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, x = Conductivity_SPC),
              alpha = 0.2) +
  geom_point(aes(x = Conductivity_SPC, y = `s(Conductivity_SPC)`),
             data = eg1, cex = 1.5, colour = "steelblue3") +
  geom_line(aes(x = Conductivity_SPC, y = est), lwd = 1.2) +
  labs(y = "Partial effect", title = "Dolichospermum ~ s(Conductivity_SPC)") +
  geom_point(aes(x = Conductivity_SPC, y = `s(Conductivity_SPC)`,
                 colour = Season, shape = Limitation), # <-- map fac to colour aesthetic
             data = eg1, cex = 1.5)
plt_Dolichospermum_Conductivity_SPC




model_Cyanobium_Conductivity_SPC <- gam(Cyanobium ~ s(Conductivity_SPC, k = 3, bs = "tp") + s(Location,bs="re") + s(Outing,bs="re"),
                               data = data, method = "REML", family= nb(link = "log"))
AIC(model_Cyanobium_Conductivity_SPC)
k.check(model_Cyanobium_Conductivity_SPC)
draw(model_Cyanobium_Conductivity_SPC, residuals = TRUE, select = "s(Conductivity_SPC)")
summary(model_Cyanobium_Conductivity_SPC)
appraise(model_Cyanobium_Conductivity_SPC, method = 'simulate')
variance_comp(model_Cyanobium_Conductivity_SPC)
resids = DHARMa::simulateResiduals(fittedModel = model_Cyanobium_Conductivity_SPC, plot = FALSE)
plot(resids) #This is the better way to deteremine if a good model
#Conductivity_SPCeratuer has no effect on cyanoboum 

sm <- smooth_estimates(model_Cyanobium_Conductivity_SPC) %>%
  add_confint()
sm
eg1 <- data %>% filter(Conductivity_SPC >= 0) %>%
  add_partial_residuals(model_Cyanobium_Conductivity_SPC)


plt_Cyanobium_Conductivity_SPC <- sm %>%
  filter(smooth == "s(Conductivity_SPC)") %>%
  ggplot() +
  geom_rug(aes(x = Conductivity_SPC),
           data = eg1,
           sides = "b", length = grid::unit(0.02, "npc")) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, x = Conductivity_SPC),
              alpha = 0.2) +
  geom_point(aes(x = Conductivity_SPC, y = `s(Conductivity_SPC)`),
             data = eg1, cex = 1.5, colour = "steelblue3") +
  geom_line(aes(x = Conductivity_SPC, y = est), lwd = 1.2) +
  labs(y = "Partial effect", title = "Cyanobium ~ s(Conductivity_SPC)") +
  geom_point(aes(x = Conductivity_SPC, y = `s(Conductivity_SPC)`,
                 colour = Season, shape = Limitation), # <-- map fac to colour aesthetic
             data = eg1, cex = 1.5)
plt_Cyanobium_Conductivity_SPC




#
model_Vulcanococcus_Conductivity_SPC <- gam(Vulcanococcus ~ s(Conductivity_SPC, k = 3) + s(Location,bs="re") + s(Outing,bs="re"),
                                   data = data, method = "REML", family= nb(link = "log"))
AIC(model_Vulcanococcus_Conductivity_SPC)
k.check(model_Vulcanococcus_Conductivity_SPC)
draw(model_Vulcanococcus_Conductivity_SPC, residuals = TRUE, select = "s(Conductivity_SPC)")
summary(model_Vulcanococcus_Conductivity_SPC)
appraise(model_Vulcanococcus_Conductivity_SPC, method = 'simulate')
variance_comp(model_Vulcanococcus_Conductivity_SPC)
resids = DHARMa::simulateResiduals(fittedModel = model_Vulcanococcus_Conductivity_SPC, plot = FALSE)
plot(resids) #This is the better way to deteremine if a good model
#Vulcanococcus likes it clear

sm <- smooth_estimates(model_Vulcanococcus_Conductivity_SPC) %>%
  add_confint()
sm
eg1 <- data %>% filter(Conductivity_SPC >= 0) %>%
  add_partial_residuals(model_Vulcanococcus_Conductivity_SPC)


plt_Vulcanococcus_Conductivity_SPC <- sm %>%
  filter(smooth == "s(Conductivity_SPC)") %>%
  ggplot() +
  geom_rug(aes(x = Conductivity_SPC),
           data = eg1,
           sides = "b", length = grid::unit(0.02, "npc")) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, x = Conductivity_SPC),
              alpha = 0.2) +
  geom_point(aes(x = Conductivity_SPC, y = `s(Conductivity_SPC)`),
             data = eg1, cex = 1.5, colour = "steelblue3") +
  geom_line(aes(x = Conductivity_SPC, y = est), lwd = 1.2) +
  labs(y = "Partial effect", title = "Vulcanococcus ~ s(Conductivity_SPC)") +
  geom_point(aes(x = Conductivity_SPC, y = `s(Conductivity_SPC)`,
                 colour = Season, shape = Limitation), # <-- map fac to colour aesthetic
             data = eg1, cex = 1.5)
plt_Vulcanococcus_Conductivity_SPC


#####TRP_mgL####



model_Microcystis_TRP_mgL <- gam(Microcystis ~ s(TRP_mgL,k=10, bs = "cr") + s(Location,bs="re") + s(Outing,bs="re"),
                                          data = data, method = "REML", family= nb(link = "log"))
AIC(model_Microcystis_TRP_mgL)
k.check(model_Microcystis_TRP_mgL)
draw(model_Microcystis_TRP_mgL, residuals = TRUE, select = "s(TRP_mgL)")
summary(model_Microcystis_TRP_mgL)
appraise(model_Microcystis_TRP_mgL, method = 'simulate')
variance_comp(model_Microcystis_TRP_mgL)
resids = DHARMa::simulateResiduals(fittedModel = model_Microcystis_TRP_mgL, plot = FALSE)
plot(resids) #This is the better way to deteremine if a good model
#Abundance increses with incread DINDIP

sm <- smooth_estimates(model_Microcystis_TRP_mgL) %>%
  add_confint()
sm
eg1 <- data %>% filter(TRP_mgL >= 0) %>%
  add_partial_residuals(model_Microcystis_TRP_mgL)

#this one wont plot
plt_Microcystis_TRP_mgL <- sm %>%
  filter(smooth == "s(TRP_mgL)") %>%
  ggplot() +
  geom_rug(aes(x = TRP_mgL),
           data = eg1,
           sides = "b", length = grid::unit(0.02, "npc")) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, x = TRP_mgL),
              alpha = 0.2) +
  geom_point(aes(x = TRP_mgL, y = `s(TRP_mgL)`),
             data = eg1, cex = 1.5, colour = "steelblue3") +
  geom_line(aes(x = TRP_mgL, y = est), lwd = 1.2) +
  labs(y = "Partial effect", title = "Microcystis ~ s(TRP_mgL)") +
  geom_point(aes(x = TRP_mgL, y = `s(TRP_mgL)`,
                 colour = Season, shape = Limitation), # <-- map fac to colour aesthetic
             data = eg1, cex = 1.5)
plt_Microcystis_TRP_mgL



model_Raphidiopsis_TRP_mgL <- gam(Raphidiopsis ~ s(TRP_mgL,k=3, bs = "cr") + s(Location,bs="re") + s(Outing,bs="re"),
                                           data = data, method = "REML", family= nb(link = "log"))
AIC(model_Raphidiopsis_TRP_mgL)
k.check(model_Raphidiopsis_TRP_mgL)
draw(model_Raphidiopsis_TRP_mgL, residuals = TRUE, select = "s(TRP_mgL)")
summary(model_Raphidiopsis_TRP_mgL)
appraise(model_Raphidiopsis_TRP_mgL, method = 'simulate')
variance_comp(model_Raphidiopsis_TRP_mgL)
resids = DHARMa::simulateResiduals(fittedModel = model_Raphidiopsis_TRP_mgL, plot = FALSE)
plot(resids) #This is the better way to deteremine if a good model
#Abundance decreases with increasing TRP_mgL, it likes low photic depth??

sm <- smooth_estimates(model_Raphidiopsis_TRP_mgL) %>%
  add_confint()
sm
eg1 <- data %>% filter(TRP_mgL >= 0) %>%
  add_partial_residuals(model_Raphidiopsis_TRP_mgL)


plt_Raphidiopsis_TRP_mgL <- sm %>%
  filter(smooth == "s(TRP_mgL)") %>%
  ggplot() +
  geom_rug(aes(x = TRP_mgL),
           data = eg1,
           sides = "b", length = grid::unit(0.02, "npc")) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, x = TRP_mgL),
              alpha = 0.2) +
  geom_point(aes(x = TRP_mgL, y = `s(TRP_mgL)`),
             data = eg1, cex = 1.5, colour = "steelblue3") +
  geom_line(aes(x = TRP_mgL, y = est), lwd = 1.2) +
  labs(y = "Partial effect", title = "Raphidiopsis ~ s(TRP_mgL)") +
  geom_point(aes(x = TRP_mgL, y = `s(TRP_mgL)`,
                 colour = Season, shape = Limitation), # <-- map fac to colour aesthetic
             data = eg1, cex = 1.5)
plt_Raphidiopsis_TRP_mgL


model_Cuspidothrix_TRP_mgL <- gam(Cuspidothrix ~ s(TRP_mgL,k=10, bs = "cr") + s(Location,bs="re") + s(Outing,bs="re"),
                                  data = data, method = "REML", family= nb(link = "log"))
AIC(model_Cuspidothrix_TRP_mgL)
k.check(model_Cuspidothrix_TRP_mgL)
draw(model_Cuspidothrix_TRP_mgL, residuals = TRUE, select = "s(TRP_mgL)")
summary(model_Cuspidothrix_TRP_mgL)
appraise(model_Cuspidothrix_TRP_mgL, method = 'simulate')
variance_comp(model_Cuspidothrix_TRP_mgL)
resids = DHARMa::simulateResiduals(fittedModel = model_Cuspidothrix_TRP_mgL, plot = FALSE)
plot(resids) #This is the better way to deteremine if a good model
#Abundance decreases with increasing TRP_mgL, it likes low photic depth??

sm <- smooth_estimates(model_Cuspidothrix_TRP_mgL) %>%
  add_confint()
sm
eg1 <- data %>% filter(TRP_mgL >= 0) %>%
  add_partial_residuals(model_Cuspidothrix_TRP_mgL)


plt_Cuspidothrix_TRP_mgL <- sm %>%
  filter(smooth == "s(TRP_mgL)") %>%
  ggplot() +
  geom_rug(aes(x = TRP_mgL),
           data = eg1,
           sides = "b", length = grid::unit(0.02, "npc")) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, x = TRP_mgL),
              alpha = 0.2) +
  geom_point(aes(x = TRP_mgL, y = `s(TRP_mgL)`),
             data = eg1, cex = 1.5, colour = "steelblue3") +
  geom_line(aes(x = TRP_mgL, y = est), lwd = 1.2) +
  labs(y = "Partial effect", title = "Cuspidothrix ~ s(TRP_mgL)") +
  geom_point(aes(x = TRP_mgL, y = `s(TRP_mgL)`,
                 colour = Season, shape = Limitation), # <-- map fac to colour aesthetic
             data = eg1, cex = 1.5)
plt_Cuspidothrix_TRP_mgL

#
model_Dolichospermum_TRP_mgL <- gam(Dolichospermum ~ s(TRP_mgL,k = 10, bs="cr") + s(Location,bs="re") + s(Outing,bs="re"),
                                             data = data, method = "REML", family= nb(link = "log"))
AIC(model_Dolichospermum_TRP_mgL)
k.check(model_Dolichospermum_TRP_mgL)
draw(model_Dolichospermum_TRP_mgL, residuals = TRUE, select = "s(TRP_mgL)")
summary(model_Dolichospermum_TRP_mgL)
appraise(model_Dolichospermum_TRP_mgL, method = 'simulate')
variance_comp(model_Dolichospermum_TRP_mgL)
resids = DHARMa::simulateResiduals(fittedModel = model_Dolichospermum_TRP_mgL, plot = FALSE)
plot(resids) #This is the better way to deteremine if a good model
#how confusing lol

sm <- smooth_estimates(model_Dolichospermum_TRP_mgL) %>%
  add_confint()
sm
eg1 <- data %>% filter(TRP_mgL >= 0) %>%
  add_partial_residuals(model_Dolichospermum_TRP_mgL)


plt_Dolichospermum_TRP_mgL <- sm %>%
  filter(smooth == "s(TRP_mgL)") %>%
  ggplot() +
  geom_rug(aes(x = TRP_mgL),
           data = eg1,
           sides = "b", length = grid::unit(0.02, "npc")) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, x = TRP_mgL),
              alpha = 0.2) +
  geom_point(aes(x = TRP_mgL, y = `s(TRP_mgL)`),
             data = eg1, cex = 1.5, colour = "steelblue3") +
  geom_line(aes(x = TRP_mgL, y = est), lwd = 1.2) +
  labs(y = "Partial effect", title = "Dolichospermum ~ s(TRP_mgL)") +
  geom_point(aes(x = TRP_mgL, y = `s(TRP_mgL)`,
                 colour = Season, shape = Limitation), # <-- map fac to colour aesthetic
             data = eg1, cex = 1.5)
plt_Dolichospermum_TRP_mgL




model_Cyanobium_TRP_mgL <- gam(Cyanobium ~ s(TRP_mgL, k = 3, bs = "tp") + s(Location,bs="re") + s(Outing,bs="re"),
                                        data = data, method = "REML", family= nb(link = "log"))
AIC(model_Cyanobium_TRP_mgL)
k.check(model_Cyanobium_TRP_mgL)
draw(model_Cyanobium_TRP_mgL, residuals = TRUE, select = "s(TRP_mgL)")
summary(model_Cyanobium_TRP_mgL)
appraise(model_Cyanobium_TRP_mgL, method = 'simulate')
variance_comp(model_Cyanobium_TRP_mgL)
resids = DHARMa::simulateResiduals(fittedModel = model_Cyanobium_TRP_mgL, plot = FALSE)
plot(resids) #This is the better way to deteremine if a good model
#TRP_mgLeratuer has no effect on cyanoboum 

sm <- smooth_estimates(model_Cyanobium_TRP_mgL) %>%
  add_confint()
sm
eg1 <- data %>% filter(TRP_mgL >= 0) %>%
  add_partial_residuals(model_Cyanobium_TRP_mgL)


plt_Cyanobium_TRP_mgL <- sm %>%
  filter(smooth == "s(TRP_mgL)") %>%
  ggplot() +
  geom_rug(aes(x = TRP_mgL),
           data = eg1,
           sides = "b", length = grid::unit(0.02, "npc")) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, x = TRP_mgL),
              alpha = 0.2) +
  geom_point(aes(x = TRP_mgL, y = `s(TRP_mgL)`),
             data = eg1, cex = 1.5, colour = "steelblue3") +
  geom_line(aes(x = TRP_mgL, y = est), lwd = 1.2) +
  labs(y = "Partial effect", title = "Cyanobium ~ s(TRP_mgL)") +
  geom_point(aes(x = TRP_mgL, y = `s(TRP_mgL)`,
                 colour = Season, shape = Limitation), # <-- map fac to colour aesthetic
             data = eg1, cex = 1.5)
plt_Cyanobium_TRP_mgL




#
model_Vulcanococcus_TRP_mgL <- gam(Vulcanococcus ~ s(TRP_mgL, k = 10, bs = "cr") + s(Location,bs="re") + s(Outing,bs="re"),
                                            data = data, method = "REML", family= nb(link = "log"))
AIC(model_Vulcanococcus_TRP_mgL)
k.check(model_Vulcanococcus_TRP_mgL)
draw(model_Vulcanococcus_TRP_mgL, residuals = TRUE, select = "s(TRP_mgL)")
summary(model_Vulcanococcus_TRP_mgL)
appraise(model_Vulcanococcus_TRP_mgL, method = 'simulate')
variance_comp(model_Vulcanococcus_TRP_mgL)
resids = DHARMa::simulateResiduals(fittedModel = model_Vulcanococcus_TRP_mgL, plot = FALSE)
plot(resids) #This is the better way to deteremine if a good model
#Vulcanococcus likes it clear

sm <- smooth_estimates(model_Vulcanococcus_TRP_mgL) %>%
  add_confint()
sm
eg1 <- data %>% filter(TRP_mgL >= 0) %>%
  add_partial_residuals(model_Vulcanococcus_TRP_mgL)


plt_Vulcanococcus_TRP_mgL <- sm %>%
  filter(smooth == "s(TRP_mgL)") %>%
  ggplot() +
  geom_rug(aes(x = TRP_mgL),
           data = eg1,
           sides = "b", length = grid::unit(0.02, "npc")) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, x = TRP_mgL),
              alpha = 0.2) +
  geom_point(aes(x = TRP_mgL, y = `s(TRP_mgL)`),
             data = eg1, cex = 1.5, colour = "steelblue3") +
  geom_line(aes(x = TRP_mgL, y = est), lwd = 1.2) +
  labs(y = "Partial effect", title = "Vulcanococcus ~ s(TRP_mgL)") +
  geom_point(aes(x = TRP_mgL, y = `s(TRP_mgL)`,
                 colour = Season, shape = Limitation), # <-- map fac to colour aesthetic
             data = eg1, cex = 1.5)
plt_Vulcanococcus_TRP_mgL


#####Depth_m####



model_Microcystis_Depth_m <- gam(Microcystis ~ s(Depth_m,k=9, bs = "cr") + s(Location,bs="re") + s(Outing,bs="re"),
                                 data = data, method = "REML", family= nb(link = "log"))
AIC(model_Microcystis_Depth_m)
k.check(model_Microcystis_Depth_m)
draw(model_Microcystis_Depth_m, residuals = TRUE, select = "s(Depth_m)")
summary(model_Microcystis_Depth_m)
appraise(model_Microcystis_Depth_m, method = 'simulate')
variance_comp(model_Microcystis_Depth_m)
resids = DHARMa::simulateResiduals(fittedModel = model_Microcystis_Depth_m, plot = FALSE)
plot(resids) #This is the better way to deteremine if a good model
#Abundance increses with incread DINDIP

sm <- smooth_estimates(model_Microcystis_Depth_m) %>%
  add_confint()
sm
eg1 <- data %>%
  add_partial_residuals(model_Microcystis_Depth_m)

#this one wont plot
plt_Microcystis_Depth_m <- sm %>%
  filter(smooth == "s(Depth_m)") %>%
  ggplot() +
  geom_rug(aes(x = Depth_m),
           data = eg1,
           sides = "b", length = grid::unit(0.02, "npc")) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, x = Depth_m),
              alpha = 0.2) +
  geom_point(aes(x = Depth_m, y = `s(Depth_m)`),
             data = eg1, cex = 1.5, colour = "steelblue3") +
  geom_line(aes(x = Depth_m, y = est), lwd = 1.2) +
  labs(y = "Partial effect", title = "Microcystis ~ s(Depth_m)") +
  geom_point(aes(x = Depth_m, y = `s(Depth_m)`,
                 colour = Season, shape = Limitation), # <-- map fac to colour aesthetic
             data = eg1, cex = 1.5)
plt_Microcystis_Depth_m



model_Raphidiopsis_Depth_m <- gam(Raphidiopsis ~ s(Depth_m,k=9, bs = "tp") + s(Location,bs="re") + s(Outing,bs="re"),
                                  data = data, method = "REML", family= nb(link = "log"))
AIC(model_Raphidiopsis_Depth_m)
k.check(model_Raphidiopsis_Depth_m)
draw(model_Raphidiopsis_Depth_m, residuals = TRUE, select = "s(Depth_m)")
summary(model_Raphidiopsis_Depth_m)
appraise(model_Raphidiopsis_Depth_m, method = 'simulate')
variance_comp(model_Raphidiopsis_Depth_m)
resids = DHARMa::simulateResiduals(fittedModel = model_Raphidiopsis_Depth_m, plot = FALSE)
plot(resids) #This is the better way to deteremine if a good model
#Abundance decreases with increasing Depth_m, it likes low photic depth??

sm <- smooth_estimates(model_Raphidiopsis_Depth_m) %>%
  add_confint()
sm
eg1 <- data %>%
  add_partial_residuals(model_Raphidiopsis_Depth_m)


plt_Raphidiopsis_Depth_m <- sm %>%
  filter(smooth == "s(Depth_m)") %>%
  ggplot() +
  geom_rug(aes(x = Depth_m),
           data = eg1,
           sides = "b", length = grid::unit(0.02, "npc")) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, x = Depth_m),
              alpha = 0.2) +
  geom_point(aes(x = Depth_m, y = `s(Depth_m)`),
             data = eg1, cex = 1.5, colour = "steelblue3") +
  geom_line(aes(x = Depth_m, y = est), lwd = 1.2) +
  labs(y = "Partial effect", title = "Raphidiopsis ~ s(Depth_m)") +
  geom_point(aes(x = Depth_m, y = `s(Depth_m)`,
                 colour = Season, shape = Limitation), # <-- map fac to colour aesthetic
             data = eg1, cex = 1.5)
plt_Raphidiopsis_Depth_m




model_Cuspidothrix_Depth_m <- gam(Cuspidothrix ~ s(Depth_m,k=9, bs = "cr") + s(Location,bs="re") + s(Outing,bs="re"),
                                  data = data, method = "REML", family= nb(link = "log"))
AIC(model_Cuspidothrix_Depth_m)
k.check(model_Cuspidothrix_Depth_m)
draw(model_Cuspidothrix_Depth_m, residuals = TRUE, select = "s(Depth_m)")
summary(model_Cuspidothrix_Depth_m)
appraise(model_Cuspidothrix_Depth_m, method = 'simulate')
variance_comp(model_Cuspidothrix_Depth_m)
resids = DHARMa::simulateResiduals(fittedModel = model_Cuspidothrix_Depth_m, plot = FALSE)
plot(resids) #This is the better way to deteremine if a good model
#Abundance decreases with increasing Depth_m, it likes low photic depth??

sm <- smooth_estimates(model_Cuspidothrix_Depth_m) %>%
  add_confint()
sm
eg1 <- data %>%
  add_partial_residuals(model_Cuspidothrix_Depth_m)


plt_Cuspidothrix_Depth_m <- sm %>%
  filter(smooth == "s(Depth_m)") %>%
  ggplot() +
  geom_rug(aes(x = Depth_m),
           data = eg1,
           sides = "b", length = grid::unit(0.02, "npc")) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, x = Depth_m),
              alpha = 0.2) +
  geom_point(aes(x = Depth_m, y = `s(Depth_m)`),
             data = eg1, cex = 1.5, colour = "steelblue3") +
  geom_line(aes(x = Depth_m, y = est), lwd = 1.2) +
  labs(y = "Partial effect", title = "Cuspidothrix ~ s(Depth_m)") +
  geom_point(aes(x = Depth_m, y = `s(Depth_m)`,
                 colour = Season, shape = Limitation), # <-- map fac to colour aesthetic
             data = eg1, cex = 1.5)
plt_Cuspidothrix_Depth_m





#
model_Dolichospermum_Depth_m <- gam(Dolichospermum ~ s(Depth_m,k = 7, bs="cr") + s(Location,bs="re") + s(Outing,bs="re"),
                                    data = data, method = "REML", family= nb(link = "log"))
AIC(model_Dolichospermum_Depth_m)
k.check(model_Dolichospermum_Depth_m)
draw(model_Dolichospermum_Depth_m, residuals = TRUE, select = "s(Depth_m)")
summary(model_Dolichospermum_Depth_m)
appraise(model_Dolichospermum_Depth_m, method = 'simulate')
variance_comp(model_Dolichospermum_Depth_m)
resids = DHARMa::simulateResiduals(fittedModel = model_Dolichospermum_Depth_m, plot = FALSE)
plot(resids) #This is the better way to deteremine if a good model
#how confusing lol

sm <- smooth_estimates(model_Dolichospermum_Depth_m) %>%
  add_confint()
sm
eg1 <- data %>%
  add_partial_residuals(model_Dolichospermum_Depth_m)


plt_Dolichospermum_Depth_m <- sm %>%
  filter(smooth == "s(Depth_m)") %>%
  ggplot() +
  geom_rug(aes(x = Depth_m),
           data = eg1,
           sides = "b", length = grid::unit(0.02, "npc")) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, x = Depth_m),
              alpha = 0.2) +
  geom_point(aes(x = Depth_m, y = `s(Depth_m)`),
             data = eg1, cex = 1.5, colour = "steelblue3") +
  geom_line(aes(x = Depth_m, y = est), lwd = 1.2) +
  labs(y = "Partial effect", title = "Dolichospermum ~ s(Depth_m)") +
  geom_point(aes(x = Depth_m, y = `s(Depth_m)`,
                 colour = Season, shape = Limitation), # <-- map fac to colour aesthetic
             data = eg1, cex = 1.5)
plt_Dolichospermum_Depth_m




model_Cyanobium_Depth_m <- gam(Cyanobium ~ s(Depth_m, k = 3, bs = "tp") + s(Location,bs="re") + s(Outing,bs="re"),
                               data = data, method = "REML", family= nb(link = "log"))
AIC(model_Cyanobium_Depth_m)
k.check(model_Cyanobium_Depth_m)
draw(model_Cyanobium_Depth_m, residuals = TRUE, select = "s(Depth_m)")
summary(model_Cyanobium_Depth_m)
appraise(model_Cyanobium_Depth_m, method = 'simulate')
variance_comp(model_Cyanobium_Depth_m)
resids = DHARMa::simulateResiduals(fittedModel = model_Cyanobium_Depth_m, plot = FALSE)
plot(resids) #This is the better way to deteremine if a good model
#Depth_meratuer has no effect on cyanoboum 

sm <- smooth_estimates(model_Cyanobium_Depth_m) %>%
  add_confint()
sm
eg1 <- data %>%
  add_partial_residuals(model_Cyanobium_Depth_m)


plt_Cyanobium_Depth_m <- sm %>%
  filter(smooth == "s(Depth_m)") %>%
  ggplot() +
  geom_rug(aes(x = Depth_m),
           data = eg1,
           sides = "b", length = grid::unit(0.02, "npc")) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, x = Depth_m),
              alpha = 0.2) +
  geom_point(aes(x = Depth_m, y = `s(Depth_m)`),
             data = eg1, cex = 1.5, colour = "steelblue3") +
  geom_line(aes(x = Depth_m, y = est), lwd = 1.2) +
  labs(y = "Partial effect", title = "Cyanobium ~ s(Depth_m)") +
  geom_point(aes(x = Depth_m, y = `s(Depth_m)`,
                 colour = Season, shape = Limitation), # <-- map fac to colour aesthetic
             data = eg1, cex = 1.5)
plt_Cyanobium_Depth_m




#
model_Vulcanococcus_Depth_m <- gam(Vulcanococcus ~ s(Depth_m, k = 5) + s(Location,bs="re") + s(Outing,bs="re"),
                                   data = data, method = "REML", family= nb(link = "log"))
AIC(model_Vulcanococcus_Depth_m)
k.check(model_Vulcanococcus_Depth_m)
draw(model_Vulcanococcus_Depth_m, residuals = TRUE, select = "s(Depth_m)")
summary(model_Vulcanococcus_Depth_m)
appraise(model_Vulcanococcus_Depth_m, method = 'simulate')
variance_comp(model_Vulcanococcus_Depth_m)
resids = DHARMa::simulateResiduals(fittedModel = model_Vulcanococcus_Depth_m, plot = FALSE)
plot(resids) #This is the better way to deteremine if a good model
#Vulcanococcus likes it clear

sm <- smooth_estimates(model_Vulcanococcus_Depth_m) %>%
  add_confint()
sm
eg1 <- data %>%
  add_partial_residuals(model_Vulcanococcus_Depth_m)


plt_Vulcanococcus_Depth_m <- sm %>%
  filter(smooth == "s(Depth_m)") %>%
  ggplot() +
  geom_rug(aes(x = Depth_m),
           data = eg1,
           sides = "b", length = grid::unit(0.02, "npc")) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, x = Depth_m),
              alpha = 0.2) +
  geom_point(aes(x = Depth_m, y = `s(Depth_m)`),
             data = eg1, cex = 1.5, colour = "steelblue3") +
  geom_line(aes(x = Depth_m, y = est), lwd = 1.2) +
  labs(y = "Partial effect", title = "Vulcanococcus ~ s(Depth_m)") +
  geom_point(aes(x = Depth_m, y = `s(Depth_m)`,
                 colour = Season, shape = Limitation), # <-- map fac to colour aesthetic
             data = eg1, cex = 1.5)
plt_Vulcanococcus_Depth_m














#####NOx_mgL####



model_Microcystis_NOx_mgL <- gam(Microcystis ~ s(NOx_mgL,k=9, bs = "cr") + s(Location,bs="re") + s(Outing,bs="re"),
                                 data = data, method = "REML", family= nb(link = "log"))
AIC(model_Microcystis_NOx_mgL)
k.check(model_Microcystis_NOx_mgL)
draw(model_Microcystis_NOx_mgL, residuals = TRUE, select = "s(NOx_mgL)")
summary(model_Microcystis_NOx_mgL)
appraise(model_Microcystis_NOx_mgL, method = 'simulate')
variance_comp(model_Microcystis_NOx_mgL)
resids = DHARMa::simulateResiduals(fittedModel = model_Microcystis_NOx_mgL, plot = FALSE)
plot(resids) #This is the better way to deteremine if a good model
#Abundance increses with incread DINDIP

sm <- smooth_estimates(model_Microcystis_NOx_mgL) %>%
  add_confint()
sm
eg1 <- data %>% 
  add_partial_residuals(model_Microcystis_NOx_mgL)

#this one wont plot
plt_Microcystis_NOx_mgL <- sm %>%
  filter(smooth == "s(NOx_mgL)") %>%
  ggplot() +
  geom_rug(aes(x = NOx_mgL),
           data = eg1,
           sides = "b", length = grid::unit(0.02, "npc")) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, x = NOx_mgL),
              alpha = 0.2) +
  geom_point(aes(x = NOx_mgL, y = `s(NOx_mgL)`),
             data = eg1, cex = 1.5, colour = "steelblue3") +
  geom_line(aes(x = NOx_mgL, y = est), lwd = 1.2) +
  labs(y = "Partial effect", title = "s(NOx_mgL)") +
  geom_point(aes(x = NOx_mgL, y = `s(NOx_mgL)`,
                 colour = Season, shape = Limitation), # <-- map fac to colour aesthetic
             data = eg1, cex = 1.5)
plt_Microcystis_NOx_mgL



model_Raphidiopsis_NOx_mgL <- gam(Cuspidothrix ~ s(NH3_mgL,k=9, bs = "cr") + s(Location,bs="re") + s(Outing,bs="re"),
                                  data = data, method = "REML", family= nb(link = "log"))
AIC(model_Raphidiopsis_NOx_mgL)
k.check(model_Raphidiopsis_NOx_mgL)
draw(model_Raphidiopsis_NOx_mgL, residuals = TRUE, select = "s(depth_thing)")
summary(model_Raphidiopsis_NOx_mgL)
appraise(model_Raphidiopsis_NOx_mgL, method = 'simulate')
variance_comp(model_Raphidiopsis_NOx_mgL)
resids = DHARMa::simulateResiduals(fittedModel = model_Raphidiopsis_NOx_mgL, plot = FALSE)
plot(resids) #This is the better way to deteremine if a good model
#Abundance decreases with increasing NOx_mgL, it likes low photic depth??

sm <- smooth_estimates(model_Raphidiopsis_NOx_mgL) %>%
  add_confint()
sm
eg1 <- data %>% filter(NH3_mgL >= 0) %>%
  add_partial_residuals(model_Raphidiopsis_NOx_mgL)


plt_Raphidiopsis_NOx_mgL <- sm %>%
  filter(smooth == "s(NH3_mgL)") %>%
  ggplot() +
  geom_rug(aes(x = NH3_mgL),
           data = eg1,
           sides = "b", length = grid::unit(0.02, "npc")) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, x = NH3_mgL),
              alpha = 0.2) +
  geom_point(aes(x = NH3_mgL, y = `s(NH3_mgL)`),
             data = eg1, cex = 1.5, colour = "steelblue3") +
  geom_line(aes(x = NH3_mgL, y = est), lwd = 1.2) +
  labs(y = "Partial effect", title = "s(NH4_NOX)") +
  geom_point(aes(x = NH3_mgL, y = `s(NH3_mgL)`,
                 colour = Season, shape = Limitation), # <-- map fac to colour aesthetic
             data = eg1, cex = 1.5)
plt_Raphidiopsis_NOx_mgL

#
model_Dolichospermum_NOx_mgL <- gam(Dolichospermum ~ s(NOx_mgL,k = 5, bs="cr") + s(Location,bs="re") + s(Outing,bs="re"),
                                    data = data, method = "REML", family= nb(link = "log"))
AIC(model_Dolichospermum_NOx_mgL)
k.check(model_Dolichospermum_NOx_mgL)
draw(model_Dolichospermum_NOx_mgL, residuals = TRUE, select = "s(NOx_mgL)")
summary(model_Dolichospermum_NOx_mgL)
appraise(model_Dolichospermum_NOx_mgL, method = 'simulate')
variance_comp(model_Dolichospermum_NOx_mgL)
resids = DHARMa::simulateResiduals(fittedModel = model_Dolichospermum_NOx_mgL, plot = FALSE)
plot(resids) #This is the better way to deteremine if a good model
#how confusing lol

sm <- smooth_estimates(model_Dolichospermum_NOx_mgL) %>%
  add_confint()
sm
eg1 <- data %>%
  add_partial_residuals(model_Dolichospermum_NOx_mgL)


plt_Dolichospermum_NOx_mgL <- sm %>%
  filter(smooth == "s(NOx_mgL)") %>%
  ggplot() +
  geom_rug(aes(x = NOx_mgL),
           data = eg1,
           sides = "b", length = grid::unit(0.02, "npc")) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, x = NOx_mgL),
              alpha = 0.2) +
  geom_point(aes(x = NOx_mgL, y = `s(NOx_mgL)`),
             data = eg1, cex = 1.5, colour = "steelblue3") +
  geom_line(aes(x = NOx_mgL, y = est), lwd = 1.2) +
  labs(y = "Partial effect", title = "s(NOx_mgL)") +
  geom_point(aes(x = NOx_mgL, y = `s(NOx_mgL)`,
                 colour = Season, shape = Limitation), # <-- map fac to colour aesthetic
             data = eg1, cex = 1.5)
plt_Dolichospermum_NOx_mgL




model_Cyanobium_NOx_mgL <- gam(Cyanobium ~ s(NOx_mgL, k = 3, bs = "tp") + s(Location,bs="re") + s(Outing,bs="re"),
                               data = data, method = "REML", family= nb(link = "log"))
AIC(model_Cyanobium_NOx_mgL)
k.check(model_Cyanobium_NOx_mgL)
draw(model_Cyanobium_NOx_mgL, residuals = TRUE, select = "s(NOx_mgL)")
summary(model_Cyanobium_NOx_mgL)
appraise(model_Cyanobium_NOx_mgL, method = 'simulate')
variance_comp(model_Cyanobium_NOx_mgL)
resids = DHARMa::simulateResiduals(fittedModel = model_Cyanobium_NOx_mgL, plot = FALSE)
plot(resids) #This is the better way to deteremine if a good model
#NOx_mgLeratuer has no effect on cyanoboum 

sm <- smooth_estimates(model_Cyanobium_NOx_mgL) %>%
  add_confint()
sm
eg1 <- data %>%
  add_partial_residuals(model_Cyanobium_NOx_mgL)


plt_Cyanobium_NOx_mgL <- sm %>%
  filter(smooth == "s(NOx_mgL)") %>%
  ggplot() +
  geom_rug(aes(x = NOx_mgL),
           data = eg1,
           sides = "b", length = grid::unit(0.02, "npc")) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, x = NOx_mgL),
              alpha = 0.2) +
  geom_point(aes(x = NOx_mgL, y = `s(NOx_mgL)`),
             data = eg1, cex = 1.5, colour = "steelblue3") +
  geom_line(aes(x = NOx_mgL, y = est), lwd = 1.2) +
  labs(y = "Partial effect", title = "s(NOx_mgL)") +
  geom_point(aes(x = NOx_mgL, y = `s(NOx_mgL)`,
                 colour = Season, shape = Limitation), # <-- map fac to colour aesthetic
             data = eg1, cex = 1.5)
plt_Cyanobium_NOx_mgL




#
model_Vulcanococcus_NOx_mgL <- gam(Vulcanococcus ~ s(NOx_mgL, k = 5) + s(Location,bs="re") + s(Outing,bs="re"),
                                   data = data, method = "REML", family= nb(link = "log"))
AIC(model_Vulcanococcus_NOx_mgL)
k.check(model_Vulcanococcus_NOx_mgL)
draw(model_Vulcanococcus_NOx_mgL, residuals = TRUE, select = "s(NOx_mgL)")
summary(model_Vulcanococcus_NOx_mgL)
appraise(model_Vulcanococcus_NOx_mgL, method = 'simulate')
variance_comp(model_Vulcanococcus_NOx_mgL)
resids = DHARMa::simulateResiduals(fittedModel = model_Vulcanococcus_NOx_mgL, plot = FALSE)
plot(resids) #This is the better way to deteremine if a good model
#Vulcanococcus likes it clear

sm <- smooth_estimates(model_Vulcanococcus_NOx_mgL) %>%
  add_confint()
sm
eg1 <- data %>%
  add_partial_residuals(model_Vulcanococcus_NOx_mgL)


plt_Vulcanococcus_NOx_mgL <- sm %>%
  filter(smooth == "s(NOx_mgL)") %>%
  ggplot() +
  geom_rug(aes(x = NOx_mgL),
           data = eg1,
           sides = "b", length = grid::unit(0.02, "npc")) +
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci, x = NOx_mgL),
              alpha = 0.2) +
  geom_point(aes(x = NOx_mgL, y = `s(NOx_mgL)`),
             data = eg1, cex = 1.5, colour = "steelblue3") +
  geom_line(aes(x = NOx_mgL, y = est), lwd = 1.2) +
  labs(y = "Partial effect", title = "s(NOx_mgL)") +
  geom_point(aes(x = NOx_mgL, y = `s(NOx_mgL)`,
                 colour = Season, shape = Limitation), # <-- map fac to colour aesthetic
             data = eg1, cex = 1.5)
plt_Vulcanococcus_NOx_mgL










#####
library(patchwork)


plt_Cuspidothrix_temp
plt_Cuspidothrix_Photic_Depth
plt_Cuspidothrix_Depth_m
plt_Cuspidothrix_Conductivity_SPC
plt_Cuspidothrix_DIN
plt_Cuspidothrix_TRP_mgL
plt_Cuspidothrix_DIN_DIP

plt_Microcystis_temp
plt_Microcystis_Photic_Depth
plt_Microcystis_Depth_m
plt_Microcystis_Conductivity_SPC
plt_Microcystis_DIN
plt_Microcystis_TRP_mgL
plt_Microcystis_DIN_DIP

plt_Dolichospermum_temp
plt_Dolichospermum_Photic_Depth
plt_Dolichospermum_Depth_m
plt_Dolichospermum_Conductivity_SPC
plt_Dolichospermum_DIN
plt_Dolichospermum_TRP_mgL
plt_Dolichospermum_DIN_DIP


plt_Raphidiopsis_temp
plt_Raphidiopsis_Photic_Depth
plt_Raphidiopsis_Depth_m
plt_Raphidiopsis_Conductivity_SPC
plt_Raphidiopsis_DIN
plt_Raphidiopsis_TRP_mgL
plt_Raphidiopsis_DIN_DIP

plt_Cyanobium_temp
plt_Cyanobium_Photic_Depth
plt_Cyanobium_Depth_m
plt_Cyanobium_Conductivity_SPC
plt_Cyanobium_DIN
plt_Cyanobium_TRP_mgL
plt_Cyanobium_DIN_DIP

plt_Vulcanococcus_temp
plt_Vulcanococcus_Photic_Depth
plt_Vulcanococcus_Depth_m
plt_Vulcanococcus_Conductivity_SPC
plt_Vulcanococcus_DIN
plt_Vulcanococcus_TRP_mgL
plt_Vulcanococcus_DIN_DIP




plt_Cuspidothrix_DIN








(plt_Cuspidothrix_temp + theme(plot.title = element_text(size = 9)) + labs(x = NULL) | plt_Microcystis_temp + theme(plot.title = element_text(size = 9)) + labs(y = NULL, x = NULL) | plt_Dolichospermum_temp + theme(plot.title = element_text(size = 9)) + labs(y = NULL, x = NULL) | plt_Raphidiopsis_temp + theme(plot.title = element_text(size = 9)) + labs(y = NULL, x = NULL)) /
(plt_Cuspidothrix_Photic_Depth + theme(plot.title = element_text(size = 9))+ labs(x = NULL) | plt_Microcystis_Photic_Depth + theme(plot.title = element_text(size = 9)) + labs(y = NULL, x = NULL) | plt_Dolichospermum_Photic_Depth + theme(plot.title = element_text(size = 9)) + labs(y = NULL, x = NULL) | plt_Raphidiopsis_Photic_Depth + theme(plot.title = element_text(size = 9)) + labs(y = NULL, x = NULL)) /
(plt_Cuspidothrix_Depth_m + theme(plot.title = element_text(size = 9))+ labs(x = NULL) | plt_Microcystis_Depth_m + theme(plot.title = element_text(size = 9)) + labs(y = NULL, x = NULL) | plt_Dolichospermum_Depth_m + theme(plot.title = element_text(size = 9)) + labs(y = NULL, x = NULL) | plt_Raphidiopsis_Depth_m + theme(plot.title = element_text(size = 9)) + labs(y = NULL, x = NULL))  + plot_layout(guides = 'collect')  + plot_annotation(tag_levels = 'A')
#(plt_Cuspidothrix_Conductivity_SPC + theme(plot.title = element_text(size = 9))+ labs(x = NULL) | plt_Microcystis_Conductivity_SPC + theme(plot.title = element_text(size = 9)) + labs(y = NULL, x = NULL) | plt_Dolichospermum_Conductivity_SPC + theme(plot.title = element_text(size = 9)) + labs(y = NULL, x = NULL) | plt_Raphidiopsis_Conductivity_SPC + theme(plot.title = element_text(size = 9)) + labs(y = NULL, x = NULL)) + plot_layout(guides = 'collect')  + plot_annotation(tag_levels = 'A')

  
(plt_Cuspidothrix_DIN + theme(plot.title = element_text(size = 9))+ labs(x = NULL) | plt_Microcystis_DIN + theme(plot.title = element_text(size = 9)) + labs(y = NULL, x = NULL) | plt_Dolichospermum_DIN + theme(plot.title = element_text(size = 9)) + labs(y = NULL, x = NULL) | plt_Raphidiopsis_DIN + theme(plot.title = element_text(size = 9)) + labs(y = NULL, x = NULL)) /
(plt_Cuspidothrix_TRP_mgL + theme(plot.title = element_text(size = 9))+ labs(x = NULL) | plt_Microcystis_TRP_mgL + theme(plot.title = element_text(size = 9)) + labs(y = NULL, x = NULL) | plt_Dolichospermum_TRP_mgL + theme(plot.title = element_text(size = 9)) + labs(y = NULL, x = NULL) | plt_Raphidiopsis_TRP_mgL + theme(plot.title = element_text(size = 9)) + labs(y = NULL, x = NULL)) /
(plt_Cuspidothrix_DIN_DIP + theme(plot.title = element_text(size = 9))+ labs(x = NULL) | plt_Microcystis_DIN_DIP + theme(plot.title = element_text(size = 9)) + labs(y = NULL, x = NULL) |  plt_Dolichospermum_DIN_DIP + theme(plot.title = element_text(size = 9)) + labs(y = NULL, x = NULL) | plt_Raphidiopsis_DIN_DIP + theme(plot.title = element_text(size = 9)) + labs(y = NULL, x = NULL)) + plot_layout(guides = 'collect')  + plot_annotation(tag_levels = 'A')



(plt_Cyanobium_temp + theme(plot.title = element_text(size = 9)) + labs(x = NULL) | plt_Cyanobium_Photic_Depth + theme(plot.title = element_text(size = 9)) + labs(y = NULL, x = NULL) | plt_Cyanobium_Depth_m + theme(plot.title = element_text(size = 9)) + labs(y = NULL, x = NULL)) /
(plt_Vulcanococcus_temp + theme(plot.title = element_text(size = 9))+ labs(x = NULL) | plt_Vulcanococcus_Photic_Depth + theme(plot.title = element_text(size = 9)) + labs(y = NULL, x = NULL) | plt_Vulcanococcus_Depth_m + theme(plot.title = element_text(size = 9)) + labs(y = NULL, x = NULL)) + plot_layout(guides = 'collect')  + plot_annotation(tag_levels = 'A')


(plt_Cyanobium_DIN + theme(plot.title = element_text(size = 9))+ labs(x = NULL) | plt_Cyanobium_TRP_mgL + theme(plot.title = element_text(size = 9)) + labs(y = NULL, x = NULL) | plt_Cyanobium_DIN_DIP + theme(plot.title = element_text(size = 9)) + labs(y = NULL, x = NULL)) /
(plt_Vulcanococcus_DIN + theme(plot.title = element_text(size = 9))+ labs(x = NULL) | plt_Vulcanococcus_TRP_mgL + theme(plot.title = element_text(size = 9)) + labs(y = NULL, x = NULL) |  plt_Vulcanococcus_DIN_DIP + theme(plot.title = element_text(size = 9)) + labs(y = NULL, x = NULL)) + plot_layout(guides = 'collect')  + plot_annotation(tag_levels = 'A')











(plt_Cyanobium_temp + theme(plot.title = element_text(size = 9)) | plt_Cyanobium_DIN_DIP + theme(plot.title = element_text(size = 9)) + labs(y = NULL) | plt_Cyanobium_DIN + theme(plot.title = element_text(size = 9)) + labs(y = NULL) | plt_Cyanobium_Photic_Depth + theme(plot.title = element_text(size = 9)) + labs(y = NULL)) / 
(plt_Cyanobium_Conductivity_SPC + theme(plot.title = element_text(size = 9)) | plt_Vulcanococcus_DIN + theme(plot.title = element_text(size = 9)) + labs(y = NULL) | plt_Vulcanococcus_Photic_Depth + theme(plot.title = element_text(size = 9)) + labs(y = NULL) | plt_Vulcanococcus_DIN_DIP + theme(plot.title = element_text(size = 9)) + labs(y = NULL)) / 
(plt_Vulcanococcus_TRP_mgL + theme(plot.title = element_text(size = 9)) | plt_Vulcanococcus_Conductivity_SPC + theme(plot.title = element_text(size = 9)) + labs(y = NULL) | plt_Vulcanococcus_Depth_m+ theme(plot.title = element_text(size = 9)) + labs(y = NULL)) + plot_layout(guides = 'collect')  + plot_annotation(tag_levels = 'A')


(plt_Dolichospermum_TRP_mgL + theme(plot.title = element_text(size = 9)) | plt_Dolichospermum_Photic_Depth + theme(plot.title = element_text(size = 9)) + labs(y = NULL) | plt_Dolichospermum_Depth_m + theme(plot.title = element_text(size = 9)) + labs(y = NULL) | plt_Dolichospermum_Conductivity_SPC + theme(plot.title = element_text(size = 9))) /
(plt_Dolichospermum_DIN_DIP + theme(plot.title = element_text(size = 9)) | plt_Microcystis_DIN_DIP + theme(plot.title = element_text(size = 9)) + labs(y = NULL) | plt_Microcystis_temp + theme(plot.title = element_text(size = 9)) + labs(y = NULL) | plt_Microcystis_Photic_Depth + theme(plot.title = element_text(size = 9)) + labs(y = NULL)) /
(plt_Raphidiopsis_DIN + theme(plot.title = element_text(size = 9)) | plt_Raphidiopsis_DIN_DIP + theme(plot.title = element_text(size = 9)) + labs(y = NULL) | plt_Raphidiopsis_Photic_Depth + theme(plot.title = element_text(size = 9)) + labs(y = NULL) | plt_Raphidiopsis_Depth_m + theme(plot.title = element_text(size = 9)) + labs(y = NULL)) /
(plt_Raphidiopsis_temp + theme(plot.title = element_text(size = 9)) | plt_Cyanobium_Depth_m + theme(plot.title = element_text(size = 9)) + labs(y = NULL)  | plt_Cyanobium_TRP_mgL + theme(plot.title = element_text(size = 9)) + labs(y = NULL) | plt_Vulcanococcus_temp + theme(plot.title = element_text(size = 9)) + labs(y = NULL)) + plot_layout(guides = 'collect')  + plot_annotation(tag_levels = 'A')


#####


model <- gam(Dolichospermum ~ s(Depth_m,k = 7, bs="cr"),
                                    data = data, method = "REML", family= nb(link = "log"))



new_df <- with(data,
               tibble(Depth_m = gratia:::seq_min_max(Depth_m, n = 200)), tibble(Location = gratia:::seq_min_max(Location, n = 200)))

fs <- fitted_samples(model, n = 10000, newdata = new_df,
                     seed = 42)

## a function to find the env value where abundance is maximum
max_loc <- function(row, abund, env) {
  r <- row[which.max(abund)]
  env[r]
}

max_post <- fs %>%
  group_by(draw) %>%
  summarise(env = max_loc(row, fitted, new_df$Depth_m))

summ <- max_post %>%
  summarise(mean = mean(env), median = median(env),
            "q2.5" = quantile(env, prob = 0.025),
            "q97.5" = quantile(env, prob = 0.975))
summ

summ <- summ %>%
  mutate(abundance = 0)

## plot
data %>%
  ggplot(aes(x = Depth_m, y = Dolichospermum)) +
  geom_point() +
  geom_pointrange(data = summ,
                  aes(y = abundance, x = median,
                      xmax = q97.5, xmin = q2.5), col = "red")




#####
model_Microcystis_xxxx <- gam(Microcystis ~ s(Temp) + s(Photic_Depth) +
                                s(Location,bs="re") + s(Outing,bs="re"),
                              data = data, method = "REML", family= nb(link = "log"))
AIC(model_Microcystis_xxxx)
k.check(model_Microcystis_xxxx)
draw(model_Microcystis_xxxx, residuals = TRUE)
summary(model_Microcystis_xxxx)
appraise(model_Microcystis_xxxx, method = 'simulate')
variance_comp(model_Microcystis_xxxx)
resids = DHARMa::simulateResiduals(fittedModel = model_Microcystis_xxxx, plot = FALSE)
plot(resids) #This is the better way to deteremine if a good model
