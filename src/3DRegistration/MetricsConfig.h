#ifndef __MetricsConfig_h
#define __MetricsConfig_h

#include <string>
#include <vector>
#include <sstream>
#include <array>
#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <iomanip>

/**
 * @class MetricsConfig
 * @brief Unified configuration for all metrics (main + label).
 *
 * Supports three submission modes (flexible):
 *   1. Arrays (all-at-once): --metrics "1.0,0.5,0,0,0,0"
 *   2. Individual: --alpha 1.0 --alphaderivative 1.0 --mapercentage 0.1
 *   3. Mixed: --preset multimodal --metric-sampling "0.1,0.1,0.1,0.1,0.1,0.1"
 *
 * Main metrics: MI, NGF, MSE, GD, NC, NMI (6 total)
 * Label metric: Separate vector of per-label weights
 */
class MetricsConfig
{
public:
    /**
     * Single metric specification: weight, derivative, sampling %.
     */
    struct MetricSpec
    {
        double weight = 0.0;
        double derivative = 0.0;
        double samplingPercent = 0.1;

        void Print(const std::string& name = "", const std::string& indent = "  ") const
        {
            if (!name.empty())
                std::cout << indent << name << ": ";
            else
                std::cout << indent;
            
            std::cout << "w=" << std::fixed << std::setprecision(3) << weight
                      << ", d=" << derivative
                      << ", samp=" << samplingPercent << std::endl;
        }
    };

    /**
     * All 6 main metrics.
     * Order matters: [MI, NGF, MSE, GD, NC, NMI]
     */
    struct MainMetricsConfig
    {
        MetricSpec mi;   // alpha
        MetricSpec ngf;  // lambda
        MetricSpec mse;  // nu
        MetricSpec gd;   // rho
        MetricSpec nc;   // yota
        MetricSpec nmi;  // sigma

        /**
         * Convert to array for easier indexing.
         * Index: 0=MI, 1=NGF, 2=MSE, 3=GD, 4=NC, 5=NMI
         */
        std::array<MetricSpec, 6> ToArray() const
        {
            return {mi, ngf, mse, gd, nc, nmi};
        }

        /**
         * Convert from array.
         */
        static MainMetricsConfig FromArray(const std::array<MetricSpec, 6>& arr)
        {
            MainMetricsConfig cfg;
            cfg.mi = arr[0];
            cfg.ngf = arr[1];
            cfg.mse = arr[2];
            cfg.gd = arr[3];
            cfg.nc = arr[4];
            cfg.nmi = arr[5];
            return cfg;
        }

        void Print(const std::string& indent = "  ") const
        {
            std::cout << indent << "Main Metrics:" << std::endl;
            mi.Print("MI (alpha)", indent + "  ");
            ngf.Print("NGF (lambda)", indent + "  ");
            mse.Print("MSE (nu)", indent + "  ");
            gd.Print("GD (rho)", indent + "  ");
            nc.Print("NC (yota)", indent + "  ");
            nmi.Print("NMI (sigma)", indent + "  ");
        }
    };

    /**
     * Label metric configuration (separate from main metrics).
     */
    struct LabelMetricConfig
    {
        std::vector<double> weights;       // Per-label weights
        std::vector<double> derivatives;   // Per-label derivatives
        double samplingPercent = 0.1;      // Same for all labels

        bool IsEnabled() const
        {
            for (double w : weights)
                if (w > 1e-6)
                    return true;
            return false;
        }

        void Print(const std::string& indent = "  ") const
        {
            if (!IsEnabled())
            {
                std::cout << indent << "Label Metric: disabled" << std::endl;
                return;
            }

            std::cout << indent << "Label Metric (" << weights.size() << " structure(s)):" << std::endl;
            std::cout << indent << "  Weights: ";
            for (size_t i = 0; i < weights.size(); ++i)
            {
                if (i > 0) std::cout << ", ";
                std::cout << std::fixed << std::setprecision(3) << weights[i];
            }
            std::cout << std::endl;

            std::cout << indent << "  Sampling: " << samplingPercent << std::endl;
        }
    };

    /**
     * **TOTAL** configuration: main + label metrics.
     */
    struct Config
    {
        MainMetricsConfig main;
        LabelMetricConfig label;

        void Print() const
        {
            main.Print();
            label.Print();
        }
    };

    // ========================================================================
    // PARSERS: Main Metrics Arrays
    // ========================================================================

    /**
     * Parse main metrics weight array.
     * Format: "alpha,lambda,nu,rho,yota,sigma"
     * Example: "1.0,0.5,0.8,0,0,0"
     */
    static MainMetricsConfig ParseMainWeights(const std::string& weightsStr)
    {
        auto arr = ParseDoubleArray(weightsStr, 6, "main metrics weights");
        MainMetricsConfig cfg;
        cfg.mi.weight = arr[0];
        cfg.ngf.weight = arr[1];
        cfg.mse.weight = arr[2];
        cfg.gd.weight = arr[3];
        cfg.nc.weight = arr[4];
        cfg.nmi.weight = arr[5];
        // Auto-derive derivatives
        cfg.mi.derivative = (arr[0] > 1e-6) ? arr[0] : 0.0;
        cfg.ngf.derivative = (arr[1] > 1e-6) ? arr[1] : 0.0;
        cfg.mse.derivative = (arr[2] > 1e-6) ? arr[2] : 0.0;
        cfg.gd.derivative = (arr[3] > 1e-6) ? arr[3] : 0.0;
        cfg.nc.derivative = (arr[4] > 1e-6) ? arr[4] : 0.0;
        cfg.nmi.derivative = (arr[5] > 1e-6) ? arr[5] : 0.0;
        return cfg;
    }

    /**
     * Parse main metrics derivative array.
     * Format: "alpha_d,lambda_d,nu_d,rho_d,yota_d,sigma_d"
     */
    static MainMetricsConfig ParseMainDerivatives(const MainMetricsConfig& base,
                                                   const std::string& derivStr)
    {
        auto arr = ParseDoubleArray(derivStr, 6, "main metrics derivatives");
        MainMetricsConfig cfg = base;
        cfg.mi.derivative = arr[0];
        cfg.ngf.derivative = arr[1];
        cfg.mse.derivative = arr[2];
        cfg.gd.derivative = arr[3];
        cfg.nc.derivative = arr[4];
        cfg.nmi.derivative = arr[5];
        return cfg;
    }

    /**
     * Parse main metrics sampling array.
     * Format: "ma%,ngf%,mse%,gd%,nc%,nmi%"
     * Example: "0.1,0.1,0.1,0.1,0.1,0.1"
     */
    static MainMetricsConfig ParseMainSampling(const MainMetricsConfig& base,
                                                const std::string& samplingStr)
    {
        auto arr = ParseDoubleArray(samplingStr, 6, "main metrics sampling");
        MainMetricsConfig cfg = base;
        cfg.mi.samplingPercent = arr[0];
        cfg.ngf.samplingPercent = arr[1];
        cfg.mse.samplingPercent = arr[2];
        cfg.gd.samplingPercent = arr[3];
        cfg.nc.samplingPercent = arr[4];
        cfg.nmi.samplingPercent = arr[5];
        return cfg;
    }

    // ========================================================================
    // PARSERS: Label Metrics
    // ========================================================================

    /**
     * Parse label weights array.
     * Format: "w1,w2,w3,..."
     * Example: "0.5,0.3,0.2"
     */
    static LabelMetricConfig ParseLabelWeights(const std::string& weightsStr,
                                                double defaultSampling = 0.1)
    {
        auto weights = ParseDoubleArray(weightsStr, -1, "label weights");  // -1 = any length
        
        LabelMetricConfig cfg;
        cfg.weights = weights;
        cfg.samplingPercent = defaultSampling;
        
        // Auto-derive derivatives
        for (double w : weights)
            cfg.derivatives.push_back((w > 1e-6) ? w : 0.0);
        
        return cfg;
    }

    /**
     * Parse label derivatives array.
     * Format: "d1,d2,d3,..."
     */
    static LabelMetricConfig ParseLabelDerivatives(const LabelMetricConfig& base,
                                                    const std::string& derivStr)
    {
        auto derivs = ParseDoubleArray(derivStr, -1, "label derivatives");
        
        LabelMetricConfig cfg = base;
        if (derivs.size() != base.weights.size())
            throw std::runtime_error("Label derivatives count (" + std::to_string(derivs.size()) +
                                   ") must match weights count (" + std::to_string(base.weights.size()) + ")");
        cfg.derivatives = derivs;
        return cfg;
    }

    // ========================================================================
    // PRESETS (convenience)
    // ========================================================================

    /**
     * Get preset main metrics config.
     * 'multimodal': MI + NGF for different modalities
     * 'singlemodal': MSE + NC for same modality
     * 'rigid': Balanced combination
     */
    static MainMetricsConfig GetMainPreset(const std::string& presetName)
    {
        MainMetricsConfig cfg;
        
        std::string lower = presetName;
        std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);

        if (lower == "multimodal")
        {
            cfg.mi.weight = 1.0;
            cfg.ngf.weight = 0.5;
            cfg.mi.derivative = 1.0;
            cfg.ngf.derivative = 0.5;
        }
        else if (lower == "singlemodal")
        {
            cfg.mse.weight = 1.0;
            cfg.nc.weight = 0.5;
            cfg.mse.derivative = 1.0;
            cfg.nc.derivative = 0.5;
        }
        else if (lower == "rigid")
        {
            cfg.mi.weight = 1.0;
            cfg.ngf.weight = 0.5;
            cfg.mse.weight = 0.25;
            cfg.mi.derivative = 1.0;
            cfg.ngf.derivative = 0.5;
            cfg.mse.derivative = 0.25;
        }
        // else: all zeros (custom)

        return cfg;
    }

    /**
     * Create scalar label preset (single weight for all labels, to be expanded later).
     */
    static LabelMetricConfig GetLabelPreset(double scalarWeight, unsigned int numLabels = 0)
    {
        LabelMetricConfig cfg;
        
        if (numLabels > 0)
        {
            // Expand to vector for known label count
            for (unsigned int i = 0; i < numLabels; ++i)
            {
                cfg.weights.push_back(scalarWeight);
                cfg.derivatives.push_back((scalarWeight > 1e-6) ? scalarWeight : 0.0);
            }
        }
        else
        {
            // Scalar form (to be expanded later)
            cfg.weights.push_back(scalarWeight);
            cfg.derivatives.push_back((scalarWeight > 1e-6) ? scalarWeight : 0.0);
        }
        
        return cfg;
    }

    // ========================================================================
    // HELPER: Merge individual parameters into config
    // ========================================================================

    /**
     * Detect conflicts between array and individual parameter submission.
     * Warns user if both --metrics and --alpha (or other individual params) are provided.
     *
     * @param hasMetricsArray  True if --metrics or --metric-derivatives provided
     * @param alpha, etc.      Individual parameter values (-1 = not provided)
     * @param verbose          If true, print conflict info
     * @return String describing conflicts (empty if none)
     */
    static std::string DetectConflicts(bool hasMetricsArray,
                                        double alpha = -1, double lambda = -1,
                                        double nu = -1, double rho = -1,
                                        double yota = -1, double sigma = -1,
                                        bool verbose = true)
    {
        if (!hasMetricsArray)
            return "";  // No conflict if using individual params only

        std::vector<std::string> conflicts;
        if (alpha >= 0) conflicts.push_back("alpha");
        if (lambda >= 0) conflicts.push_back("lambda");
        if (nu >= 0) conflicts.push_back("nu");
        if (rho >= 0) conflicts.push_back("rho");
        if (yota >= 0) conflicts.push_back("yota");
        if (sigma >= 0) conflicts.push_back("sigma");

        if (conflicts.empty())
            return "";  // No conflicts

        std::string msg = "WARNING: Both --metrics array and individual parameters provided: ";
        for (size_t i = 0; i < conflicts.size(); ++i)
        {
            if (i > 0) msg += ", ";
            msg += "--" + conflicts[i];
        }
        msg += ". Individual parameters take precedence.";

        if (verbose)
            std::cout << "\n" << msg << std::endl << std::endl;

        return msg;
    }

    /**
     * Override/merge individual metric parameters into config.
     * Used when both array (--metrics) and individual (--alpha) are provided.
     * Individual takes precedence.
     *
     * NOTE: Call DetectConflicts() first to warn user if desired.
     */
    static MainMetricsConfig MergeIndividual(const MainMetricsConfig& base,
                                             double alpha = -1, double alphaDeriv = -1,
                                             double lambda = -1, double lambdaDeriv = -1,
                                             double nu = -1, double nuDeriv = -1,
                                             double rho = -1, double rhoDeriv = -1,
                                             double yota = -1, double yotaDeriv = -1,
                                             double sigma = -1, double sigmaDeriv = -1)
    {
        MainMetricsConfig cfg = base;
        
        if (alpha >= 0) cfg.mi.weight = alpha;
        if (alphaDeriv >= 0) cfg.mi.derivative = alphaDeriv;
        
        if (lambda >= 0) cfg.ngf.weight = lambda;
        if (lambdaDeriv >= 0) cfg.ngf.derivative = lambdaDeriv;
        
        if (nu >= 0) cfg.mse.weight = nu;
        if (nuDeriv >= 0) cfg.mse.derivative = nuDeriv;
        
        if (rho >= 0) cfg.gd.weight = rho;
        if (rhoDeriv >= 0) cfg.gd.derivative = rhoDeriv;
        
        if (yota >= 0) cfg.nc.weight = yota;
        if (yotaDeriv >= 0) cfg.nc.derivative = yotaDeriv;
        
        if (sigma >= 0) cfg.nmi.weight = sigma;
        if (sigmaDeriv >= 0) cfg.nmi.derivative = sigmaDeriv;
        
        return cfg;
    }

private:
    /**
     * Helper: Parse comma-separated double array.
     * @param str Input string
     * @param expectedLength Expected length (-1 for any length)
     * @param context Error message context
     */
    static std::vector<double> ParseDoubleArray(const std::string& str,
                                                 int expectedLength = -1,
                                                 const std::string& context = "array")
    {
        std::vector<double> result;
        
        if (str.empty())
            throw std::runtime_error("Empty " + context);

        std::stringstream ss(str);
        std::string token;

        while (std::getline(ss, token, ','))
        {
            // Trim whitespace
            token.erase(0, token.find_first_not_of(" \t"));
            token.erase(token.find_last_not_of(" \t") + 1);

            if (token.empty())
                throw std::runtime_error("Empty value in " + context);

            try
            {
                result.push_back(std::stod(token));
            }
            catch (const std::exception& e)
            {
                throw std::runtime_error("Failed to parse '" + token + "' in " + context + ": " + e.what());
            }
        }

        if (expectedLength > 0 && static_cast<int>(result.size()) != expectedLength)
            throw std::runtime_error(context + " has " + std::to_string(result.size()) +
                                   " values, expected " + std::to_string(expectedLength));

        return result;
    }
};

#endif  // __MetricsConfig_h
