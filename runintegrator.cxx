// ***************************************************************
// This file was created using the bat-project script
// for project integrator.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#include <BAT/BCAux.h>
#include <BAT/BCLog.h>
#include <BAT/BCMath.h>

#include "Integrator.h"

struct TestModel : public BCModel
{
    TestModel(const std::string name, unsigned n) : BCModel(name)
        {
            for (unsigned i = 1; i <= n; ++i)
                AddParameter("x_" + std::to_string(i), -5, 5.);
            SetPriorConstantAll();
        }

    double LogLikelihood(const std::vector<double>& par)
        {
            double LL = 0;
            for (auto p : par)
                LL += BCMath::LogGaus(p, 0., 1., true);
            return LL;
        }
};

int main()
{
    // // open log file
    // BCLog::OpenLog("log.txt", BCLog::detail, BCLog::detail);
    BCLog::SetLogLevel(BCLog::nothing);

    TestModel test("test_model", 10);
    test.GetParameter(0).Fix(0);
    test.GetParameter(5).Fix(0);
    test.SetPrecision(BCEngineMCMC::kMedium);
    test.SetNIterationsRun(1e5);
    test.SetFlagFillHistograms(false);
    test.MarginalizeAll(BCIntegrate::kMargMetropolis);
    // test.PrintSummary();

    Integrator int_test(test);
    int_test.SetPrecision(BCEngineMCMC::kMedium);
    int_test.SetNIterationsRun(1e5);
    int_test.MarginalizeAll(BCIntegrate::kMargMetropolis);
    // int_test.PrintSummary();

    double sampling_integral = approx_posterior_integral(int_test);

    std::cout << "approx posterior integral = " << sampling_integral << std::endl;

    std::cout << "ratio = " << int_test.GetIntegral()
              << " +- " << int_test.GetIntegralUncertainty()
              << std::endl;
    
    std::cout << "integral = " << int_test.GetIntegral(sampling_integral)
              << " +- " << int_test.GetIntegralUncertainty(sampling_integral)
              << " (" << (100. * int_test.GetIntegralUncertainty() / int_test.GetIntegral()) << " %)"
              << "\tin [" << int_test.GetIntegral(sampling_integral) - 3 * int_test.GetIntegralUncertainty(sampling_integral)
              << ", " << int_test.GetIntegral(sampling_integral) + 3 * int_test.GetIntegralUncertainty(sampling_integral)
              << "] with 95% credibility"
              << std::endl;

    // BCLog::SetLogLevel(BCLog::detail);
    // BCLog::OutSummary("Exiting");
    // BCLog::CloseLog();

    return 0;
}
