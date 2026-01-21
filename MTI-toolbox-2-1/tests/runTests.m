function runTests(createWithArtifacts,htmlCoverageReport)
% runTests - Run the tests for all files in the project and generate optional a test report.
%
%   Input parameters:
%   - createWithArtifacts: [logical] Indicates whether to create test artifacts or not.
%
%   Output:
%   - None
%
%   Example:
%   runTests(true) % Run tests and generate test report with artifacts
%   runTests(false) % Run tests without generating test artifacts
%   Other:
%   Documentation: https://de.mathworks.com/help/matlab/matlab_prog/generate-artifacts-using-matlab-unit-test-plugins.html
%   this file need an update after enable the CICD pipeline
arguments
    createWithArtifacts logical = false
    htmlCoverageReport logical = false
end
    import matlab.unittest.TestSuite;
    import matlab.unittest.TestRunner
    import matlab.unittest.plugins.TestReportPlugin
    import matlab.unittest.Verbosity
    import matlab.unittest.plugins.CodeCoveragePlugin
    import matlab.unittest.plugins.XMLPlugin
    import matlab.unittest.plugins.codecoverage.CoberturaFormat
    import matlab.unittest.plugins.codecoverage.CoverageReport



    cd(erase(fileparts(which(mfilename)),'tests')); % set current path to parent folder
    %Run all tests in this project
    if (batchStartupOptionUsed)
        suite = TestSuite.fromFolder(strcat(pwd,'/tests'), 'IncludingSubfolders', true); %Create testsuite from testfolder
    else
        suite = TestSuite.fromFolder(strcat(pwd,'\tests'), 'IncludingSubfolders', true); %Create testsuite from testfolder
    end

    if createWithArtifacts
        mkdir("tests/artifacts")
        runner = TestRunner.withTextOutput("OutputDetail",Verbosity.Detailed);
        
        if htmlCoverageReport
            runner.addPlugin(CodeCoveragePlugin.forFolder(fileparts(which('mlgreyest.m')),'IncludingSubfolders',true, ...
            'Producing',CoverageReport('tests/artifacts')))
        else
            runner.addPlugin(CodeCoveragePlugin.forFolder(fileparts(which('mlgreyest.m')),'IncludingSubfolders',true, ...
            'Producing',CoberturaFormat('tests/artifacts/cobertura.xml')))
        end
        
        runner.addPlugin(XMLPlugin.producingJUnitFormat('tests/artifacts/results.xml'))
        

        %if (batchStartupOptionUsed)
        %    pdfFile = "tests/artifacts/TestReport.pdf";
        %else
        %    pdfFile = "tests\artifacts\TestReport.pdf";
        %end
        %p1 = TestReportPlugin.producingPDF(pdfFile,"IncludingCommandWindowText",true);
        %runner.addPlugin(p1)

        result = runner.run(suite);
    else

        result = run(suite);
    end

    if (any([result.Failed]) )
        disp("At least one test failed")
        if (batchStartupOptionUsed)
            exit(1)
        end
    else
        disp(strcat("All (",string(length(result)),") Tests passed"))
    end

    if (batchStartupOptionUsed)
        exit(0)
    end
end
