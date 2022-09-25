function tests = wsOFMM_main_test
    tests = functiontests(localfunctions);
end

function setupOnce(testCase)
    % Input directory
    testCase.TestData.in_dir = "C:/Users/Lang/Documents/Harvard/Research/Briana/supRPC/wsOFMM/Toy_Example/";   
    % Output directory 
    testCase.TestData.out_dir = "C:/Users/Lang/Documents/Harvard/Research/Briana/supRPC/wsOFMM/Toy_Example/";  
    testCase.TestData.pop_data = importdata(strcat(in_dir, 'simdata_scen2_iter1.mat'));
    testCase.TestData.samp_data = importdata(strcat(in_dir, 'simdata_scen14_iter1_samp1.mat'));

end

function test_wtd_get_data_vars(testCase)
    
end

function teardownOnce(testCase)
    delete(testCase.TestData.in_dir);
    delete(testCase.TestData.out_dir);
end