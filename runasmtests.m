function runasmtests()
%runasmtests Runs all test cases in the tests folder
p = mfilename('fullpath');
asmdir = fileparts(p);
testsuite = matlab.unittest.TestSuite.fromFolder(fullfile(asmdir, 'tests'));
run(testsuite)
end