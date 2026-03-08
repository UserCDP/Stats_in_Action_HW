library(readr)
NHANES <- read_csv("NHANES_age_prediction 3.csv")

data = NHANES[,c("DIQ010","age_group","RIDAGEYR","RIAGENDR","PAQ605","BMXBMI","LBXGLU","LBXGLT","LBXIN")]
#colnames(data) = c("Diabete","age_group","Age","Sex","Phys_activ","BMI","Glu","Glu2h","BIL")

#rmarkdown::paged_table(data)

load("liver_data.rda")

gene_cor <- cor(liver[ , colnames(liver) != "cholesterol"],
                liver$cholesterol)
gene_cor

# the gene with the highest positive correlation
max(gene_cor, na.rm = TRUE)

# the gene with the highest overall correlation
max(abs(gene_cor), na.rm = TRUE)

# and its name
max_index <- which.max(abs(gene_cor))
gene_name <- rownames(gene_cor)[max_index]

gene_expr <- liver[[gene_name]]

image(gene_cor)

plot(gene_expr,
     liver$cholesterol,
     xlab = gene_name,
     ylab = "Cholesterol",
     main = paste("Cholesterol vs", gene_name),
     pch = 19,
     col = "steelblue")

model <- lm(cholesterol ~ gene_expr, data = liver)
summary(model)

abline(model, col = "red", lwd = 2)
coef(model)

## Question 2
gene_matrix <- liver[, colnames(liver) != "cholesterol"]
p_values <- apply(gene_matrix, 2, function(gene) {
  summary(lm(liver$cholesterol ~ gene))$coefficients[2,4]
})
length(p_values)

## Question 3
m <- length(p_values)   # should be 3116

p_sorted <- sort(p_values)

plot(p_sorted,
     ylab = "Ordered p-values",
     xlab = "Rank",
     main = "Ordered p-values vs Rank",
     pch = 19,
     cex = 0.5)

# Add line y = x/m
abline(a = 0, b = 1/m, col = "red", lwd = 2)

## Question 4
# Adjust the p-values for multiple testing to limit the expected proportion of
# false positives among significant results
p_adj_fdr <- p.adjust(p_values, method = "BH")

# Filter the significant discoveries
discoveries_fdr <- names(p_adj_fdr)[p_adj_fdr < 0.05]

# Count the discoveries after filtering
length(discoveries_fdr)

## Question 5
p_adj_bonf <- p.adjust(p_values, method = "bonferroni")

discoveries_bonf <- names(p_adj_bonf)[p_adj_bonf < 0.05]

length(discoveries_bonf)
