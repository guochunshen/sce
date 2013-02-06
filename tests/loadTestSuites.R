#run all the test function after buid and reload the package

#define test a testSuite for basic object creation
testsuite.createObject <- defineTestSuite("createObject",
                    dirs="./tests/objects")


runTestSuite(testsuite.createObject)