
gf = read.csv("/tmp/rminctestdata/CIVET_TEST.csv")
gf = civet.getAllFilenames(gf,"ID","POND","/tmp/rminctestdata/CIVET","TRUE","1.1.12")
gf = civet.readAllCivetFiles("/tmp/rminctestdata/AAL.csv",gf)


context("anatMean")

#Calculate mean
sink("/dev/null"); vm <- anatMean(gf$lobeThickness); sink();

test_that("anatMean", {
    for (j in 1:dim(gf$lobeThickness)[2]) {
    		expect_that(mean(gf$lobeThickness[,j]),is_equivalent_to(vm[j]))
	}
})

context("anatSum")

#Calculate sum
sink("/dev/null"); vs <- anatSum(gf$lobeThickness); sink();

test_that("anatSum", {
    for (j in 1:dim(gf$lobeThickness)[2]) {
    		expect_that(sum(gf$lobeThickness[,j]),is_equivalent_to(vs[j]))
	}	
})

context("anatVar")

#Calculate variance
sink("/dev/null"); vv <- anatVar(gf$lobeThickness); sink();

test_that("anatVar", {
    for (j in 1:dim(gf$lobeThickness)[2]) {
    		expect_that(var(gf$lobeThickness[,j]),is_equivalent_to(vv[j]))
	}	
})


context("anatSd")

#Calculate standard deviation
sink("/dev/null"); vsd <- anatSd(gf$lobeThickness); sink();

test_that("anatSd", {
    for (j in 1:dim(gf$lobeThickness)[2]) {
    		expect_that(sd(gf$lobeThickness[,j]),is_equivalent_to(vsd[j]))
	}	
})

