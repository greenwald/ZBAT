// ***************************************************************
// This file was created using the bat-project script.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#ifndef __BAT__INTEGRATOR__H
#define __BAT__INTEGRATOR__H

#include <BAT/BCModel.h>

#include <TDecompChol.h>
#include <TMath.h>
#include <TMatrixDSym.h>
#include <TVectorD.h>

#include <algorithm>
#include <string>
#include <vector>

// This is a Integrator header file.
// Model source code is located in file integrator/Integrator.cxx

// ---------------------------------------------------------
class Integrator : public BCModel
{

public:

    /// Constructor
    Integrator(BCModel& model) : BCModel(model.GetName() + "_integrator"),
                                 fModel(model),
                                 fApproxPosteriorInverseCovariance(fModel.GetNParameters()),
                                 fApproxPosteriorCoverianceDeterminant(0)
    {
        // get covariance as ROOT matrix
        TMatrixDSym C(fModel.GetNParameters());
        for (size_t i = 0; i < fModel.GetNParameters(); ++i)
            for (size_t j = 0; j < fModel.GetNParameters(); ++j)
                C[i][j] = fModel.GetStatistics().covariance[i][j];
        // invert it, storing the determinant
        fApproxPosteriorInverseCovariance = C.Invert(&fApproxPosteriorCoverianceDeterminant);

        // Copy model's paramters
        fParameters = fModel.GetParameters();
        SetPriorConstantAll();
        
        // Create Observable for model's posterior / approx posterior
        AddObservable("posterior_to_approx_posterior", 0, 100);

        SetFlagFillHistograms(false);
    }
    
    // \return flat
    double LogLikelihood(const std::vector<double>& pars)
    {
        // x[i] := pars[i] - mean[i]
        std::vector<double> x;
        x.reserve(pars.size());
        std::transform(pars.begin(), pars.end(), fModel.GetStatistics().mean.begin(), std::back_inserter(x), std::minus<double>());

        double LL = 0;
        for (size_t i = 0; i < x.size(); ++i) {
            for (size_t j = 0; j < x.size(); ++j)
                LL += x[i] * fApproxPosteriorInverseCovariance[i][j] * x[j];
        }
        return -0.5 * LL - 0.5 * x.size() * log(2. * TMath::Pi()) - 0.5 * log(fApproxPosteriorCoverianceDeterminant);
    }

    /// Store (posterior of model - pedestal) / (approx posterior)
    void CalculateObservables(const std::vector<double>& pars)
    //    { GetObservable(0) = exp(fModel.LogEval(pars) - fModel.GetStatistics().probability_mean - GetLogProbx(GetCurrentChain())); }
    { GetObservable(0) = exp(fModel.LogEval(pars) - fModel.GetStatistics().probability_mean - LogEval(pars)); }

    const BCModel& model() const
    { return fModel; }
    
    /// \return inverse covariance of approx posterior
    const TMatrixDSym& GetApproxPosteriorInverseCovariance() const
    { return fApproxPosteriorInverseCovariance; }

    /// \return determinant of covariance of approx posterior
    double GetApproxPosteriorCovarianceDeterminant() const
    { return fApproxPosteriorCoverianceDeterminant; }

    // \return integral
    double GetIntegral(double prior_integral = 1) const
    { return GetStatistics().mean[GetNParameters()] * exp(fModel.GetStatistics().probability_mean) * prior_integral; }

    // \return uncertainty on integral
    double GetIntegralUncertainty(double prior_integral = 1) const
    { return sqrt(GetStatistics().variance[GetNParameters()] / GetStatistics().n_samples) * exp(fModel.GetStatistics().probability_mean) * prior_integral; }
    
protected:

    BCModel& fModel;

    TMatrixDSym fApproxPosteriorInverseCovariance;

    double fApproxPosteriorCoverianceDeterminant;
};

/// \return integral of a multivariate Gaussian
/// \param mean mean of Gaussian
/// \param Cinv inverse of covariance matrix of Gaussian
/// \param Cdet determinant of covariance matrix of Gaussian
/// \param xmin lower limit of integral range
/// \param xmax upper limit of integral range
inline double multivariate_gaussian_integral(const std::vector<double>& mean,
                                             const TMatrixDSym& Cinv, double Cdet,
                                             const std::vector<double>& xmin,
                                             const std::vector<double>& xmax)
{
    // get Cholesky decomposition
    //   Cinv = Lc Lc^T
    TDecompChol chol_dec;
    chol_dec.SetMatrix(Cinv);
    chol_dec.Decompose();
    // ROOT gives upper-triangular matrix
    TMatrixD L = chol_dec.GetU();
    // We need the lower-triangular one
    L.Transpose(L);

    //   Cinv = L D L^T
    // from cholesky decomposition
    //   S := diag(L)
    //   D = S^2
    //   L -> L S^-1

    // store vector of diag(S) = sqrt(D)
    std::vector<double> S;
    S.reserve(L.GetNrows());
    for (int i = 0; i < L.GetNrows(); ++i) {
        S.push_back(L[i][i]);
        for (int j = 0; j < L.GetNcols(); ++j)
            L[i][j] /= S.back();
    }
    // we need L^T
    L = L.Transpose(L);

    // transform into diagonal basis and divide by sqrt(2 / D_ii)
    std::vector<double> diag_min;
    std::vector<double> diag_max;
    std::vector<double> diag_mean;

    double integral = pow(4., -0.5 * L.GetNrows()) / sqrt(Cdet);
    for (int i = 0; i < L.GetNrows(); ++i) {

        // transform into diagonal basis
        double amin = 0;
        double amax = 0;
        for (int j = 0; j < L.GetNcols(); ++j) {
            amin += L[i][j] * (xmin[j] - mean[j]);
            amax += L[i][j] * (xmax[j] - mean[j]);
        }

        integral *= (TMath::Erf(amax * S[i] / sqrt(2)) - TMath::Erf(amin * S[i] / sqrt(2))) / S[i];
    }

    return integral;
    
}

/// \return integral of approx posterior
/// \param m Integrator model to integrate approx posterior of
/// \param xmin lower limit of integration
/// \param xmax upper limit of integration
inline double approx_posterior_integral(const Integrator& m,
                                        std::vector<double> xmin = std::vector<double>(),
                                        std::vector<double> xmax = std::vector<double>())
{
    if (xmin.empty()) {
        xmin.reserve(m.GetNParameters());
        for (size_t i = 0; i < m.GetNParameters(); ++i)
            xmin.push_back(m.GetParameter(i).GetLowerLimit());
    }

    if (xmax.empty()) {
        xmax.reserve(m.GetNParameters());
        for (size_t i = 0; i < m.GetNParameters(); ++i)
            xmax.push_back(m.GetParameter(i).GetUpperLimit());
    }

    return multivariate_gaussian_integral(m.model().GetStatistics().mean, m.GetApproxPosteriorInverseCovariance(),
                                          m.GetApproxPosteriorCovarianceDeterminant(), xmin, xmax);
}



#endif
