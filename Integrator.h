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
                                 fApproxPosteriorInverseCovariance(fModel.GetNFreeParameters())
    {
        // get covariance as ROOT matrix
        TMatrixDSym C(fModel.GetNFreeParameters());
        size_t I = 0;
        for (size_t i = 0; i < fModel.GetNParameters(); ++i) {
            if (fModel.GetParameter(i).Fixed()) continue;
            size_t J = 0;
            for (size_t j = 0; j < fModel.GetNParameters(); ++j) {
                if (fModel.GetParameter(j).Fixed()) continue;
                C[I][J] = fModel.GetStatistics().covariance[i][j];
                C[J][I] = fModel.GetStatistics().covariance[i][j];
                ++J;
            }
            ++I;
        }
        
        // invert it, storing the determinant
        double det;
        fApproxPosteriorInverseCovariance = C.Invert(&det);

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
        // x[i] := pars[i] - mean[i], if par(i) not fixed
        std::vector<double> x;
        x.reserve(GetNFreeParameters());
        for (size_t i = 0; i < GetNParameters(); ++i)
            if (!GetParameter(i).Fixed())
                x.push_back(pars[i] - fModel.GetStatistics().mean[i]);

        double LL = 0;
        for (size_t i = 0; i < x.size(); ++i) {
            LL += pow(x[i], 2) * fApproxPosteriorInverseCovariance[i][i];
            for (size_t j = i + 1; j < x.size(); ++j)
                LL += 2 * x[i] * fApproxPosteriorInverseCovariance[i][j] * x[j];
        }
        return -0.5 * LL;
    }

    /// Store (posterior of model - pedestal) / (approx posterior)
    void CalculateObservables(const std::vector<double>& pars)
    {
        fModel.UpdateChainIndex(GetCurrentChain());
        GetObservable(0) = exp(fModel.LogEval(pars) - fModel.GetStatistics().probability_mean - LogEval(pars));
    }

    const BCModel& model() const
    { return fModel; }
    
    /// \return inverse covariance of approx posterior
    const TMatrixDSym& GetApproxPosteriorInverseCovariance() const
    { return fApproxPosteriorInverseCovariance; }

    // \return integral
    double GetIntegral(double prior_integral = 1) const
    { return GetStatistics().mean[GetNParameters()] * exp(fModel.GetStatistics().probability_mean) * prior_integral; }
    
    // \return uncertainty on integral
    double GetIntegralUncertainty(double prior_integral = 1) const
    { return sqrt(GetStatistics().variance[GetNParameters()] / GetStatistics().n_samples) * exp(fModel.GetStatistics().probability_mean) * prior_integral; }
    
protected:

    BCModel& fModel;

    TMatrixDSym fApproxPosteriorInverseCovariance;

};

/// \return integral of a multivariate Gaussian
/// \param mean mean of Gaussian
/// \param Cinv inverse of covariance matrix of Gaussian
/// \param xmin lower limit of integral range
/// \param xmax upper limit of integral range
inline double multivariate_gaussian_integral(const std::vector<double>& mean, const TMatrixDSym& Cinv,
                                             const std::vector<double>& xmin, const std::vector<double>& xmax,
                                             double epsilon = 1.e-3)
{
    if (static_cast<int>(mean.size()) < Cinv.GetNrows())
        throw;
    
    if (static_cast<int>(xmin.size()) < Cinv.GetNrows())
        throw;

    if (static_cast<int>(xmax.size()) < Cinv.GetNrows())
        throw;
    
    // get Cholesky decomposition
    //   Cinv = Lc Lc^T
    TDecompChol chol_dec;
    chol_dec.SetMatrix(Cinv);

    // try to decompose
    if(!chol_dec.Decompose()) {
        // if failed, nudge diagonal elements
        TMatrixDSym Cinv_nudged = Cinv;
        for (int i = 0; i < Cinv_nudged.GetNrows(); ++i)
            Cinv_nudged[i][i] *= 1 + epsilon;
        // try again
        if (!chol_dec.Decompose())
            // if fails still, throw
            throw;
    }
        
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

    double integral = pow(TMath::Pi() / 2., 0.5 * L.GetNrows());
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
        xmin.reserve(m.GetNFreeParameters());
        for (size_t i = 0; i < m.GetNParameters(); ++i)
            if (!m.GetParameter(i).Fixed())
                xmin.push_back(m.GetParameter(i).GetLowerLimit());
    }

    if (xmax.empty()) {
        xmax.reserve(m.GetNFreeParameters());
        for (size_t i = 0; i < m.GetNParameters(); ++i)
            if (!m.GetParameter(i).Fixed())
                xmax.push_back(m.GetParameter(i).GetUpperLimit());
    }

    std::vector<double> mean;
    mean.reserve(m.GetNFreeParameters());
    for (size_t i = 0; i < m.GetNParameters(); ++i)
        if (!m.GetParameter(i).Fixed())
            mean.push_back(m.model().GetStatistics().mean[i]);
    
    return multivariate_gaussian_integral(mean, m.GetApproxPosteriorInverseCovariance(), xmin, xmax);
}



#endif
