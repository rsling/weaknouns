##############################################
#
# I did this when I was young and needed the money ;)
# Seriously, although I think the results are valid,
# I would run and report my GLMMs differently now.
# - Roland, 2015-11-09
#
##############################################


rm(list=ls(all=T))

require(aplpack)
require(fmsb)
require(car)
#require(safeBinaryRegression)
library(MASS)



# ---------- FUNCTIONS ----------

# Sums all strong counts for one stem.
sum.countstr <- function(stem) {
  sum(wk[which(wk$Stem==stem),]$CountStr)
}

sum.countwk <- function(stem) {
  sum(wk[which(wk$Stem==stem),]$CountWk);
}

# Sums all weak counts for one stem.
sum.weakstr <- function(stem) {
  sum(wk[which(wk$Stem==stem),]$CountWk)
}

sum.all <- function(stem) {
  sum(wk[which(wk$Stem==stem),]$CountWk) +
    sum(wk[which(wk$Stem==stem),]$CountStr)
}

collo <- function(arow){ 
  a1 <- as.numeric(arow["CountStr"]);
  a2 <- as.numeric(arow["CountWk"]);
  b1 <- as.numeric(arow["Field3"]);
  b2 <- as.numeric(arow["Field4"]);
  fisher.test(matrix(c(a1,a2,b1,b2),nrow=2,byrow=T))$p.value;
}

siglev <- function(p)
{
  if (p <= 0.1)
    if (p <= 0.05)
      if (p <= 0.01)
        if (p <= 0.001)
          sig <- "***"
        else
          sig <- "**"
      else
        sig <- "*"
    else
      sig <- "."
  else
    sig <- ""
}

plot.coef <- function(glm, alpha = 0.1, cex=1.5, pcol="Pr(>|t|)")
{
  # Extract p values and coefs.
  p <- coef(summary(glm))[,pcol]
  c <- coef(glm)

  # Separate intercept from other coefs.
  p.int <- as.numeric(p[1])
  c.int <- as.numeric(c[1])
  p <- p[-1]
  c <- c[-1]

  # Find the order in which to put coefficients.
  c <- c[which(p <= alpha)]
  p <- p[which(p <= alpha)]
  
  # Find the order according to coef.
  order <- order(c)
  c <- c[order]
  p <- p[order]
  names <- names(c)

  # Put intercept and coefs back together.
  names <- c("(Intercept)", names)
  vals <- c(round(c.int, 3), round(c, 3))

  c <- c(as.numeric(c.int), as.numeric(c))
  p <- c(as.numeric(p.int), as.numeric(p))
  vals <- paste(vals, lapply(p, siglev), sep=" ")
  names <- paste(names, vals, sep="\n")
  
  range <- max(c)-min(c)
  ylim <- c(min(c)-0.25*range, max(c)+0.25*range)
  bars <- barplot(c, ylim=ylim, col=c("white", rep("lightgray", length(c)-1)))
  pos <- ifelse(c > 0, c+0.4*cex, c-0.4*cex)
  text(bars, pos, names, cex=cex)
}

lemma.pprint <- function(row)
{
  #print(row)
  cat(row[1], " (", row[2], "/", round(as.numeric(row[3]), 3), "±", round(as.numeric(row[4]), 3), "), ", sep="")
}

ci.prop <- function(p, n, cl = 0.05, round = 3)
{
  z <- qnorm(cl/2, lower.tail=F);
  se <- sqrt(p*(1-p)/n);
  ci <- z*se;
  c(round(ci, round), round(p-ci, round), round(p+ci,round));
}

create.rows.w <- function(row)
{
 as.data.frame(row)[, rep(1, row[13])]
}

create.rows.s <- function(row)
{
  as.data.frame(row)[, rep(1, row[14])]
}

# ---------- LEMMA FREQS ----------

f00 <- read.delim("~/workingcopies/weaknouns/data/Rweaknouns/00.csv", header=F)
f01 <- read.delim("~/workingcopies/weaknouns/data/Rweaknouns/01.csv", header=F)
f02 <- read.delim("~/workingcopies/weaknouns/data/Rweaknouns/02.csv", header=F)
f03 <- read.delim("~/workingcopies/weaknouns/data/Rweaknouns/03.csv", header=F)
f04 <- read.delim("~/workingcopies/weaknouns/data/Rweaknouns/04.csv", header=F)
f05 <- read.delim("~/workingcopies/weaknouns/data/Rweaknouns/05.csv", header=F)
f06 <- read.delim("~/workingcopies/weaknouns/data/Rweaknouns/06.csv", header=F)
f07 <- read.delim("~/workingcopies/weaknouns/data/Rweaknouns/07.csv", header=F)

f <- data.frame(cbind(f00[,1], f01[,1], f02[,1], f03[,1], f04[,1], f05[,1], f06[,1], f07[,1]))
f$Stem <- f00[,2]
f$Freq <- apply(f[,1:8], 1, sum)

# "Friese" and "Pole" aren't properly tagged by TT.
# Add adequate Frequencies.
f[which(f$Stem=="Friese"),]$Freq <- 10000
f[which(f$Stem=="Pole"),]$Freq <- 10000

f$Fpm <- f$Freq/9100
f$Lfpm <- log10(f$Fpm)

minF <- abs(min(f$Lfpm))
f$Lfpm <- f$Lfpm+minF
f <- f[,-c(1:8)]

# ---------- DATA ----------

wk <- read.table("wk.csv", sep=",", quote="", header=T)
wk$Freq <- f[match(wk$Stem, f$Stem), "Freq"]
#wk <- wk[-which(wk$Freq==0),]
wk$LogFreq <- f[match(wk$Stem, f$Stem), "Lfpm"]
wk$Loan <- factor(ifelse(wk$Loan==1,1,0))
#wk$Case <- factor(wk$Case)
wk$Case <- factor(wk$Case, levels=c("Gen", "Dat", "Acc"))
wk$Monosyll <- factor(wk$Monosyll)
wk$Loan <- factor(wk$Loan, levels=c("1", "0", "2"))
wk$Schwa <- factor(ifelse(wk$Rhyme=="E", 1, 0))
wk$Accent <- factor(wk$Accent, levels=c("Ult", "Ini", "Pen"))
wk$Sem <- factor(wk$Sem, levels=c("Hum", "Ina", "Ani"))

# All phonotactics in one factor:
# Mono, PolySchwa, PolyNult, PolyUlt
wk$Pt <- factor(ifelse(wk$Monosyll=="1", "Mono",
                      ifelse(wk$Schwa=="1", "PolySchwa",
                             ifelse(wk$Accent=="Ult", "PolyUlt", "PolyNult")
                             )
                      ),
                levels=c("PolySchwa", "PolyUlt", "Mono", "PolyNult"))


wk.lemmas <- as.character(unique(wk$Stem))
wk.lemmas.freqs <- unlist(lapply(wk.lemmas, sum.all))
wk.lemmas.wk <- unlist(lapply(wk.lemmas, sum.countwk))
wk.lemmas.str <- unlist(lapply(wk.lemmas, sum.countstr))
wk.lemfreq <- cbind(wk.lemmas, wk.lemmas.freqs, wk.lemmas.str, wk.lemmas.wk, wk.lemmas.str/(wk.lemmas.str+wk.lemmas.wk))

wk.ci <- function(row)
{
  ci <- ci.prop(as.numeric(row[5]), as.numeric(row[2]), 0.05, 4)
  ci[1]
}
wk.lemmas.ci <- unlist(apply(wk.lemfreq, 1, wk.ci))

wk.lemfreq <- cbind(wk.lemfreq, wk.lemmas.ci)

colnames(wk.lemfreq) <- c("Lemma", "Freq", "FStr", "FWk", "Prop", "95%CI")

# ---------- DUMP LIST (DESCRIPTIVE) ----------

sink("lemmas_alpha.txt")
apply(wk.lemfreq, 1, lemma.pprint)
sink()

sink("lemmas_freq.txt")
apply(wk.lemfreq[order(as.numeric(wk.lemfreq[,2]), decreasing=T),], 1, lemma.pprint)
sink()

# ---------- FILTER STEMS ----------

stems.filter <- c("Bauer", "Buchstabe", "Friede", "Gedanke", "Glaube", "Götze", "Junge", "Mensch", "Name", "Pole", "Same", "Wille", "Page", "Artist", "Resident", "Steinmetz", "Titan", "Herr")
wk.reduced <- wk[which(!(wk$Stem %in% stems.filter)),]
wk.lemmas.reduced <- wk.lemmas[which(!(wk.lemmas %in% stems.filter))]
wk.lemfreq.reduced <- wk.lemfreq[which(!(wk.lemfreq[,1] %in% stems.filter)),]

# ---------- DATA AND GLM ON SINGLE OBS ----------

wk.single.wk <- t(as.data.frame(apply(wk.reduced, 1, create.rows.w)))
wk.single.wk <- as.data.frame(wk.single.wk[sample(1:nrow(wk.single.wk), nrow(wk.single.wk)),])
wk.single.wk$Strong <- "0"

wk.single.st <- t(as.data.frame(apply(wk.reduced, 1, create.rows.s)))
wk.single.st <- as.data.frame(wk.single.st[sample(1:nrow(wk.single.st), nrow(wk.single.wk)),])
wk.single.st$Strong <- "1"

wk.single <- rbind(wk.single.st, wk.single.wk)
wk.single <- wk.single[sample(1:nrow(wk.single), nrow(wk.single)),]

wk.single$Stem <- factor(wk.single$Stem)
wk.single$Freq <- as.numeric(f[match(wk.single$Stem, f$Stem), "Freq"])
wk.single$LogFreq <- as.numeric(f[match(wk.single$Stem, f$Stem), "Lfpm"])
#wk.single$Freq <- as.numeric(wk.single$Freq)
#wk.single$LogFreq <- as.numeric(wk.single$LogFreq)
wk.single$Strong <- factor(wk.single$Strong, levels=c(1,0))
wk.single$Sem <- factor(wk.single$Sem, levels=c("Hum", "Ani", "Ina"))
wk.single$Accent <- factor(wk.single$Accent, levels=c("Ult", "Ini", "Pen"))
wk.single$Class <- factor(wk.single$Class)
wk.single$Rhyme <- factor(wk.single$Rhyme)
wk.single$Coda <- factor(wk.single$Coda)
wk.single$Syll <- as.numeric(wk.single$Syll)
wk.single$Monosyll <- factor(wk.single$Monosyll, levels=c("0", "1"))
wk.single$Case <- factor(wk.single$Case, levels=c("Gen", "Acc", "Dat"))
wk.single$CountStr <- as.numeric(wk.single$CountStr)
wk.single$CountWk <- as.numeric(wk.single$CountWk)
wk.single$PropStr <- as.numeric(wk.single$PropStr)
wk.single$Schwa <- factor(wk.single$Schwa, levels=c("0", "1"))
wk.single$Pt <- factor(wk.single$Pt, levels=c("PolySchwa", "PolyUlt", "Mono", "PolyNult"))

wk.single.red <- wk.single[sample(1:nrow(wk.single), 10000),]

wk.single.glm <- glm(Strong~Case+Sem+Pt+LogFreq, family=binomial, data=wk.single.red)
wk.single.glm.0 <- glm(Strong~1, family=binomial, data=wk.single.red)

sink("glm.txt")
print(summary(wk.single.glm))
cat("Rounded coefs\n")
print(round(coef(wk.single.glm), 3))
cat("\nOdds ratios\n")
print(round(exp(coef(wk.single.glm)), 3))
cat("\nRounded Pr(>|t|)\n")
print(round(summary(wk.single.glm)$coefficients[,"Pr(>|z|)"], 5))
cat("\n\n")
cat("            n=", nrow(wk.single.red), "\n\n")
cat("           R²=", NagelkerkeR2(wk.single.glm)$R2, "\n\n")
cat("  deviance/df=", round(summary(wk.single.glm)$deviance/summary(wk.single.glm)$df.residual, 2), "\n\n")
cat("Dispersion  φ=", summary(glm(Strong~Case+Sem+Pt+LogFreq, family=quasibinomial, data=wk.single.red))$dispersion, "\n\n")
print(vif(wk.single.glm))
cat("\n\n")
print(anova(wk.single.glm, wk.single.glm.0, test="Chisq"))
cat("\n\n")
cat("Filtered stems: ", stems.filter[order(stems.filter)], "\n")
sink()

svg("coefs.svg")
plot.coef(wk.single.glm, cex=0.8, alpha=0.11, pcol="Pr(>|z|)")
dev.off()


# ---------- PLOT ----------

svg("distribution.svg")
bagplot(as.numeric(wk.lemfreq.reduced[,3]), as.numeric(wk.lemfreq.reduced[,4]), xlim=c(0,200), ylim=c(0,8000), show.outlier=T, approx.limit=600,
        col.loophull="white", col.baghull="white", col.looppoints="black", col.bagpoints="black", transparency=F, cex.axis=1.4, cex.lab=1.4, cex.main=1.4,
        pch=19, cex=0.6, lwd=2, show.whiskers=F, xlab="raw frequency of strong form", ylab="raw frequency of weak form",
        main="Distribution of weak and strong forms per lemma")
dev.off()

# ---------- CASE DISTRIBUTION (DESCRIPTIVE) ----------

acc <- wk.reduced[which(wk.reduced$Case=="Acc"),]
dat <- wk.reduced[which(wk.reduced$Case=="Dat"),]
gen <- wk.reduced[which(wk.reduced$Case=="Gen"),]

acc.rel <- sum(na.omit(acc$CountStr)) / ( sum(na.omit(acc$CountWk)) + sum(na.omit(acc$CountStr)) )
dat.rel <- sum(na.omit(dat$CountStr)) / ( sum(na.omit(dat$CountWk)) + sum(na.omit(dat$CountStr)) )
gen.rel <- sum(na.omit(gen$CountStr)) / ( sum(na.omit(gen$CountWk)) + sum(na.omit(gen$CountStr)) )

case.dist <- matrix(
  c(sum(na.omit(acc$CountStr)), sum(na.omit(acc$CountWk)),
    sum(na.omit(dat$CountStr)), sum(na.omit(dat$CountWk)),
    sum(na.omit(gen$CountStr)), sum(na.omit(gen$CountWk))), 2, 3)
rownames(case.dist) <- c("Str", "Wk")
colnames(case.dist) <- c("Acc", "Dat", "Gen")

sink("casedist.txt")
print(case.dist)
acc.p <- case.dist[1,1]/case.dist[2,1]
dat.p <- case.dist[1,2]/case.dist[2,2]
gen.p <- case.dist[1,3]/case.dist[2,3]
acc.n <- case.dist[1,1]+case.dist[2,1]
dat.n <- case.dist[1,2]+case.dist[2,2]
gen.n <- case.dist[1,3]+case.dist[2,3]
acc.ci <- ci.prop(acc.p, acc.n, 0.01, 4)
dat.ci <- ci.prop(dat.p, dat.n, 0.01, 4)
gen.ci <- ci.prop(gen.p, gen.n, 0.01, 4)
cat("\n\nProportion Acc (99% ci): ", acc.p, " (±", acc.ci[1], ": ", acc.ci[2], " .. ", acc.ci[3], ", n=", acc.n, ")", sep="")
cat("\nProportion Dat (99% ci): ", dat.p, " (±", dat.ci[1], ": ", dat.ci[2], " .. ", dat.ci[3], ", n=", dat.n, ")", sep="")
cat("\nProportion Gen (99% ci): ", gen.p, " (±", gen.ci[1], ": ", gen.ci[2], " .. ", gen.ci[3], ", n=", gen.n, ")", sep="")
sink()

