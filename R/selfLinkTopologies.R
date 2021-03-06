GenerateSelfLinks <- function(){
df1 <- data.frame(Source = c("A"), Target = c("A"), Type = (1))
df2 <- data.frame(Source = c("A"), Target = c("A"), Type = c(2))
df3 <- data.frame(Source = c("B"), Target = c("B"), Type = c(1))
df4 <- data.frame(Source = c("B"), Target = c("B"), Type = c(2))
df5 <- data.frame(Source = c("C"), Target = c("C"), Type = c(1))
df6 <- data.frame(Source = c("C"), Target = c("C"), Type = c(2))
df7 <- data.frame(Source = c("A", "B"), Target = c("A", "B"), Type = c(1, 1))
df8 <- data.frame(Source = c("A", "B"), Target = c("A", "B"), Type = c(1, 2))
df9 <- data.frame(Source = c("A", "B"), Target = c("A", "B"), Type = c(2, 1))
df10 <- data.frame(Source = c("A", "B"), Target = c("A", "B"), Type = c(2, 2))
df11 <- data.frame(Source = c("B", "C"), Target = c("B", "C"), Type = c(1, 1))
df12 <- data.frame(Source = c("B", "C"), Target = c("B", "C"), Type = c(1, 2))
df13 <- data.frame(Source = c("B", "C"), Target = c("B", "C"), Type = c(2, 1))
df14 <- data.frame(Source = c("B", "C"), Target = c("B", "C"), Type = c(2, 2))
df15 <- data.frame(Source = c("A", "C"), Target = c("A", "C"), Type = c(1, 1))
df16 <- data.frame(Source = c("A", "C"), Target = c("A", "C"), Type = c(1, 2))
df17 <- data.frame(Source = c("A", "C"), Target = c("A", "C"), Type = c(2, 1))
df18 <- data.frame(Source = c("A", "C"), Target = c("A", "C"), Type = c(2, 2))
df19 <- data.frame(Source = c("A", "B", "C"), Target = c("A", "B", "C"), Type = c(1, 1, 1))
df20 <- data.frame(Source = c("A", "B", "C"), Target = c("A", "B", "C"), Type = c(1, 1, 2))
df21 <- data.frame(Source = c("A", "B", "C"), Target = c("A", "B", "C"), Type = c(1, 2, 1))
df22 <- data.frame(Source = c("A", "B", "C"), Target = c("A", "B", "C"), Type = c(1, 2, 2))
df23 <- data.frame(Source = c("A", "B", "C"), Target = c("A", "B", "C"), Type = c(2, 1, 1))
df24 <- data.frame(Source = c("A", "B", "C"), Target = c("A", "B", "C"), Type = c(2, 1, 2))
df25 <- data.frame(Source = c("A", "B", "C"), Target = c("A", "B", "C"), Type = c(2, 2, 1))
df26 <- data.frame(Source = c("A", "B", "C"), Target = c("A", "B", "C"), Type = c(2, 2, 2))
self.links <- list(df1, df2, df3, df4, df5, df6, df7, df8, df9, df10, df11, df12, df13, df14, df15, df16, df17, df18, df19, df20, df21, df22, df23, df24, df25, df26)

names(self.links)<-c("Source", "Target", "Type")
return (self.links)
}

