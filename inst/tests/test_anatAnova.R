context("anatAnova")

gf <- read.csv("/tmp/rminctestdata/CIVET_TEST.csv")
gf <- civet.getAllFilenames(gf,"ID","TEST","/tmp/rminctestdata/CIVET","TRUE","1.1.12")
gf <- civet.readAllCivetFiles("/tmp/rminctestdata/AAL.csv",gf)

rmincAnova = verboseRun("anatAnova(~ Sex,gf,gf$lobeThickness)",getOption("verbose"))

lobeThickness = gf$lobeThickness[,1]
Age = gf$Age
Sex = gf$Sex
rAnova = anova(lm(lobeThickness~Sex))

test_that("anatAnova Two Factors",{
	expect_that(rmincAnova[1,1],is_equivalent_to(rAnova$F[1]))
	expect_that(attr(rmincAnova,"df")[[1]][2],is_equivalent_to(rAnova$Df[2]))
	expect_that(attr(rmincAnova,"df")[[1]][1],is_equivalent_to(rAnova$Df[1]))
})

rmincAnova = verboseRun("anatAnova(~ Age*Sex,gf,gf$lobeThickness)",getOption("verbose"))

lobeThickness = gf$lobeThickness[,1]
Age = gf$Age
Sex = gf$Sex
rAnova = anova(lm(lobeThickness~Age*Sex))

test_that("anatAnova Interaction",{
	expect_that(rmincAnova[1,1],is_equivalent_to(rAnova$F[1]))
	expect_that(rmincAnova[1,2],is_equivalent_to(rAnova$F[2]))
	expect_that(rmincAnova[1,3],is_equivalent_to(rAnova$F[3]))
	expect_that(attr(rmincAnova,"df")[[1]][1],is_equivalent_to(rAnova$Df[1]))
	expect_that(attr(rmincAnova,"df")[[1]][2],is_equivalent_to(rAnova$Df[4]))
	expect_that(attr(rmincAnova,"df")[[2]][1],is_equivalent_to(rAnova$Df[2]))
	expect_that(attr(rmincAnova,"df")[[2]][2],is_equivalent_to(rAnova$Df[4]))
	expect_that(attr(rmincAnova,"df")[[3]][1],is_equivalent_to(rAnova$Df[3]))
	expect_that(attr(rmincAnova,"df")[[3]][2],is_equivalent_to(rAnova$Df[4]))
})
rmincAnova = verboseRun("anatAnova(~ Primary.Diagnosis,gf,gf$lobeThickness)",getOption("verbose"))

lobeThickness = gf$lobeThickness[,1]
Primary.Diagnosis = gf$Primary.Diagnosis
rAnova = anova(lm(lobeThickness~Primary.Diagnosis))


test_that("anatAnova Three Factors",{
	expect_that(attr(rmincAnova,"df")[[1]][2],is_equivalent_to(rAnova$Df[2]))
	expect_that(attr(rmincAnova,"df")[[1]][1],is_equivalent_to(rAnova$Df[1]))
})


rmincAnova = verboseRun("anatAnova(~Age*Primary.Diagnosis,gf,gf$lobeThickness)",getOption("verbose"))

lobeThickness = gf$lobeThickness[,1]
Primary.Diagnosis = as.factor(gf$Primary.Diagnosis)
rAnova = anova(lm(lobeThickness~Age*Primary.Diagnosis))

test_that("anatAnova Three Factors Interaction",{
         expect_that(rmincAnova[1,1],is_equivalent_to(rAnova$F[1]))
         expect_that(rmincAnova[1,2],is_equivalent_to(rAnova$F[2]))
         expect_that(rmincAnova[1,3],is_equivalent_to(rAnova$F[3]))
         expect_that(attr(rmincAnova,"df")[[1]][1],is_equivalent_to(rAnova$Df[1]))
         expect_that(attr(rmincAnova,"df")[[1]][2],is_equivalent_to(rAnova$Df[4]))
         expect_that(attr(rmincAnova,"df")[[2]][1],is_equivalent_to(rAnova$Df[2]))
         expect_that(attr(rmincAnova,"df")[[2]][2],is_equivalent_to(rAnova$Df[4]))
         expect_that(attr(rmincAnova,"df")[[3]][1],is_equivalent_to(rAnova$Df[3]))
         expect_that(attr(rmincAnova,"df")[[3]][2],is_equivalent_to(rAnova$Df[4])) 
})
