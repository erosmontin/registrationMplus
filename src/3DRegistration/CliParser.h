#ifndef __CliParser_h
#define __CliParser_h

#include <string>
#include <vector>
#include <sstream>
#include <algorithm>
#include <stdexcept>

/**
 * @class CliParser
 * @brief Unified CLI parser supporting both old (verbose) and new (array-based) formats.
 *
 * Supports automatic format detection:
 *   Old: --alpha 1.0 --lambda 0.5 --nu 0.8
 *   New: --weights "1.0,0.5,0.8,0,0,0,0"  or --preset multimodal
 *
 * Derivatives default to weight values if not specified separately.
 * Sampling percentages default to 10% for all metrics if not specified.
 */
class CliParser
{
public:
    struct MetricWeights
    {
        double alpha = 1.0;
        double lambda = 0.0;
        double nu = 0.0;
        double rho = 0.0;
        double yota = 0.0;
        double kappa = 0.0;
        double sigma = 0.0;
    };

    struct MetricDerivatives
    {
        double alphaDerivative = 1.0;
        double lambdaDerivative = 0.0;
        double nuDerivative = 0.0;
        double rhoDerivative = 0.0;
        double yotaDerivative = 0.0;
        double kappaDerivative = 0.0;
        double sigmaDerivative = 0.0;
    };

    struct SamplingPercentages
    {
        double ma = 0.1;
        double ngf = 0.1;
        double mse = 0.1;
        double gd = 0.1;
        double nc = 0.1;
        double nmi = 0.1;
        double label = 0.1;
    };

    /**
     * Parse comma-separated array of metric weights.
     * Format: "alpha,lambda,nu,rho,yota,kappa,sigma"
     * Example: "1.0,0.5,0.2,0,0,0,0"
     */
    static MetricWeights ParseWeightsArray(const std::string& weightsStr)
    {
        MetricWeights w;
        std::vector<double> values = ParseDoubleArray(weightsStr);
        
        if (values.size() != 7)
            throw std::runtime_error("Weights array must have exactly 7 values: "
                                   "alpha,lambda,nu,rho,yota,kappa,sigma");
        
        w.alpha = values[0];
        w.lambda = values[1];
        w.nu = values[2];
        w.rho = values[3];
        w.yota = values[4];
        w.kappa = values[5];
        w.sigma = values[6];
        
        return w;
    }

    /**
     * Parse comma-separated array of sampling percentages.
     * Format: "ma,ngf,mse,gd,nc,nmi,label"
     * Example: "0.1,0.1,0.1,0.1,0.1,0.1,0.1"
     */
    static SamplingPercentages ParseSamplingArray(const std::string& samplingStr)
    {
        SamplingPercentages s;
        std::vector<double> values = ParseDoubleArray(samplingStr);
        
        if (values.size() != 7)
            throw std::runtime_error("Sampling array must have exactly 7 values: "
                                   "ma,ngf,mse,gd,nc,nmi,label");
        
        s.ma = values[0];
        s.ngf = values[1];
        s.mse = values[2];
        s.gd = values[3];
        s.nc = values[4];
        s.nmi = values[5];
        s.label = values[6];
        
        return s;
    }

    /**
     * Get preset metric weights by name.
     * 'multimodal': alpha=1.0, lambda=0.5 (MI + NGF for different modalities)
     * 'singlemodal': nu=1.0, yota=0.5 (MSE + NC for same modality)
     * 'rigid': alpha=1.0, lambda=0.5, nu=0.25 (balanced, structural)
     */
    static MetricWeights GetPreset(const std::string& presetName)
    {
        MetricWeights w;
        
        std::string lower = presetName;
        std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
        
        if (lower == "multimodal")
        {
            w.alpha = 1.0;
            w.lambda = 0.5;
            w.nu = 0.0;
        }
        else if (lower == "singlemodal")
        {
            w.alpha = 0.0;
            w.lambda = 0.0;
            w.nu = 1.0;
            w.yota = 0.5;
        }
        else if (lower == "rigid")
        {
            w.alpha = 1.0;
            w.lambda = 0.5;
            w.nu = 0.25;
        }
        else if (lower == "none" || lower == "custom")
        {
            // All zeros, user will specify manually
        }
        else
        {
            throw std::runtime_error("Unknown preset: '" + presetName + 
                                   "'. Valid: multimodal, singlemodal, rigid, custom");
        }
        
        return w;
    }

    /**
     * Derive derivative weights from regular weights.
     * If a weight is non-zero, its derivative is non-zero (default to weight value).
     * This keeps value and gradient consistent and simplifies user input.
     */
    static MetricDerivatives DeriveFromWeights(const MetricWeights& w)
    {
        MetricDerivatives d;
        
        d.alphaDerivative = (w.alpha != 0.0) ? w.alpha : 0.0;
        d.lambdaDerivative = (w.lambda != 0.0) ? w.lambda : 0.0;
        d.nuDerivative = (w.nu != 0.0) ? w.nu : 0.0;
        d.rhoDerivative = (w.rho != 0.0) ? w.rho : 0.0;
        d.yotaDerivative = (w.yota != 0.0) ? w.yota : 0.0;
        d.kappaDerivative = (w.kappa != 0.0) ? w.kappa : 0.0;
        d.sigmaDerivative = (w.sigma != 0.0) ? w.sigma : 0.0;
        
        return d;
    }

    /**
     * Check if new (simplified) CLI format is being used.
     * Returns true if --weights or --preset parameters are present.
     */
    static bool IsNewFormat(const std::string& weightsParam,
                            const std::string& presetParam)
    {
        return (!weightsParam.empty() && weightsParam != "") ||
               (!presetParam.empty() && presetParam != "");
    }

private:
    /**
     * Parse comma-separated list of doubles.
     * Example: "1.0,0.5,0.2" -> {1.0, 0.5, 0.2}
     */
    static std::vector<double> ParseDoubleArray(const std::string& str)
    {
        std::vector<double> result;
        std::stringstream ss(str);
        std::string token;
        
        while (std::getline(ss, token, ','))
        {
            // Trim whitespace
            token.erase(0, token.find_first_not_of(" \t"));
            token.erase(token.find_last_not_of(" \t") + 1);
            
            if (token.empty())
                throw std::runtime_error("Empty value in comma-separated array");
            
            try
            {
                result.push_back(std::stod(token));
            }
            catch (const std::exception& e)
            {
                throw std::runtime_error(std::string("Failed to parse '") + token + 
                                       "' as double: " + e.what());
            }
        }
        
        return result;
    }
};

#endif // __CliParser_h
