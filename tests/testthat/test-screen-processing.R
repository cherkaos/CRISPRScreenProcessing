context("test-screen-processing")

test_that("Basic test", {

   resultfile <- tempfile()
   inputfile="APK-1-and-2-final.txt"
   controlStart="sample35"
   controlEnd="sample19-II"
   maxsgRNA=15
   minReadCount=100
   zscore=TRUE
   screenProcessing(inputfile,controlStart,controlEnd, resultfile, 15,100,TRUE)
   expect_true(file.exists(resultfile))
   output=read.table(resultfile, check.name=FALSE,sep='\t')

   expect_equal(dim(output),c(865,151),1)

   # TODO: load result file, check dimensions
   #      compute some stats on columns
   #      use: expect_equal(a, b, tolerance)
   #      Implement Coverage  https://github.com/r-lib/covr

   # modify input file to have some strange column names (spaces in it  or starting with numbers),
   #   use such names also for controlEnd, ...
   #   finally check the created result file if these names are there as expected.

})

test_that("Errors test", {
  expect_error(screenProcessing("definitively not a file", 0, 0, "out", 15,100,TRUE),"No such file .*")
  expect_error(screenProcessing('APK-1-and-2-final.txt',"XXX", "XXX", "out", 15,100,TRUE),"The clone/organoid .*")

})


