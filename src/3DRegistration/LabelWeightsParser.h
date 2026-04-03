#ifndef __LabelWeightsParser_h
#define __LabelWeightsParser_h

#include <string>
#include <vector>
#include <sstream>
#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <iomanip>

/**
 * @class LabelWeightsParser
 * @brief Handle scalar and vector label weights for multi-label registration.
 *
 * Labels in registration can be weighted in two ways:
 *
 *   1. SCALAR (simple): One `labelkappa` applied equally to ALL labels
 *      Pros: Simple, no need to know label count upfront
 *      Usage: --labelkappa 0.5
 *      Result: All labels weighted 0.5
 *
 *   2. VECTOR (advanced): Per-label weights for different importance
 *      Pros: Fine-grained control
 *      Usage: --label-weights "0.5,0.3,0.2"
 *      Result: Label 1 → 0.5, Label 2 → 0.3, Label 3 → 0.2
 *
 * This class auto-expands scalar to vector when a labelmap is loaded.
 */
class LabelWeightsParser
{
public:
    struct LabelWeights
    {
        std::vector<double> kappaValues;       ///< Per-label kappa weights
        std::vector<double> kappaDerivatives;  ///< Per-label derivatives
        bool isVector = false;                  ///< True if multi-label
        unsigned int numLabels = 0;             ///< Total labels in map

        /**
         * Get scalar summary (mean of all label weights).
         * Useful for UI display when user wants one number.
         */
        double GetScalarKappa() const
        {
            if (kappaValues.empty())
                return 0.0;
            if (kappaValues.size() == 1)
                return kappaValues[0];
            
            double sum = 0.0;
            for (double v : kappaValues)
                sum += v;
            return sum / kappaValues.size();
        }

        /**
         * Check if labels are enabled (any non-zero weight).
         */
        bool IsEnabled() const
        {
            for (double v : kappaValues)
                if (v > 1e-6)
                    return true;
            return false;
        }

        /**
         * Print weights to console.
         */
        void Print(const std::string& indent = "  ") const
        {
            if (kappaValues.empty())
            {
                std::cout << indent << "LabelWeights: disabled" << std::endl;
                return;
            }

            std::cout << indent << "LabelWeights (" << kappaValues.size() << " label(s)): ";
            for (size_t i = 0; i < kappaValues.size(); ++i)
            {
                if (i > 0) std::cout << ", ";
                std::cout << std::fixed << std::setprecision(3) 
                          << "κ[" << i << "]=" << kappaValues[i];
            }
            std::cout << std::endl;
        }
    };

    /**
     * Parse comma-separated label weights from string.
     *
     * @param weightsStr  Comma-separated values: "0.5,0.3,0.2"
     * @return LabelWeights with isVector=true, derivatives auto-derived
     *
     * Example:
     *   "0.5,0.3,0.2" → 3 labels with weights [0.5, 0.3, 0.2]
     */
    static LabelWeights ParseVector(const std::string& weightsStr)
    {
        LabelWeights lw;
        std::vector<double> values = ParseDoubleArray(weightsStr);

        if (values.empty())
            throw std::runtime_error("Label weights array is empty");

        lw.kappaValues = values;
        lw.numLabels = values.size();
        lw.isVector = (values.size() > 1);

        // Auto-derive: if weight > 1e-6, use it; else 0
        for (double kappa : values)
        {
            lw.kappaDerivatives.push_back((kappa > 1e-6) ? kappa : 0.0);
        }

        return lw;
    }

    /**
     * Create scalar label weight (applied to ALL labels).
     * Will be expanded to vector when labelmap is loaded.
     *
     * @param kappaValue       Weight to apply to all labels
     * @param kappaDerivative  Derivative (default: same as weight, or 0 if disabledweight is 0)
     * @return LabelWeights with isVector=false, will expand later
     *
     * Example:
     *   ScalarLabelWeights(0.5) → {kappaValues: [0.5], isVector: false}
     *   Later when 3 labels detected:
     *     → ExpandToVector() → {kappaValues: [0.5, 0.5, 0.5], isVector: true}
     */
    static LabelWeights ScalarLabelWeights(double kappaValue,
                                            double kappaDerivative = -1.0)
    {
        LabelWeights lw;
        lw.kappaValues.push_back(kappaValue);

        double deriv = (kappaDerivative < 0.0) ?
            ((kappaValue > 1e-6) ? kappaValue : 0.0) :
            kappaDerivative;
        lw.kappaDerivatives.push_back(deriv);

        lw.isVector = false;
        lw.numLabels = 0;  // Unknown until labelmap is read

        return lw;
    }

    /**
     * Expand scalar label weight to vector (one per detected label).
     *
     * Call this after you've detected the number of labels in the labelmap.
     * If lw is already a vector, this does nothing.
     *
     * @param lw              LabelWeights (scalar form)
     * @param numLabelsInMap  Number of labels detected in labelmap
     * @return Expanded LabelWeights
     *
     * Example:
     *   lw = ScalarLabelWeights(0.5)           → [0.5], isVector=false
     *   lw = ExpandToVector(lw, 3)             → [0.5, 0.5, 0.5], isVector=true
     */
    static LabelWeights ExpandToVector(const LabelWeights& lw,
                                        unsigned int numLabelsInMap)
    {
        // If already vector, return as-is
        if (lw.isVector)
            return lw;

        // If scalar with 0 labels, expand it
        if (lw.kappaValues.size() == 1)
        {
            LabelWeights expanded;
            double scalarKappa = lw.kappaValues[0];
            double scalarDeriv = lw.kappaDerivatives[0];

            for (unsigned int i = 0; i < numLabelsInMap; ++i)
            {
                expanded.kappaValues.push_back(scalarKappa);
                expanded.kappaDerivatives.push_back(scalarDeriv);
            }

            expanded.isVector = true;
            expanded.numLabels = numLabelsInMap;

            return expanded;
        }

        // Already expanded or malformed
        return lw;
    }

    /**
     * Disable all label weights (set all to 0).
     * Useful for conditional registration (no label penalty).
     */
    static LabelWeights Disabled()
    {
        return ScalarLabelWeights(0.0);
    }

    /**
     * Detect number of unique non-zero labels in a labelmap file.
     *
     * NOTE: This is a STUB. For real implementation, you need ITK includes
     *       and your existing image reading code from imageUtils.h.
     *
     * @param labelmapPath  Path to label map file (NIfTI, DICOM, etc.)
     * @return Number of unique labels (or 0 if file not found)
     *
     * Implementation pattern (use your imageUtils):
     * ```cpp
     * auto labelmap = ReadImage<LabelImageType>(labelmapPath);
     * itk::ImageRegionConstIterator<LabelImageType> it(labelmap, labelmap->GetLargestPossibleRegion());
     * std::set<LabelPixelType> uniqueLabels;
     * for (it.GoToBegin(); !it.IsAtEnd(); ++it)
     *     if (it.Get() > 0) uniqueLabels.insert(it.Get());
     * return uniqueLabels.size();
     * ```
     */
    static unsigned int DetectNumberOfLabels(const std::string& labelmapPath)
    {
        // Stub: in production, actually read the image and count labels
        if (labelmapPath.empty() || labelmapPath == "N" || labelmapPath == "None")
            return 0;

        std::cout << "  [WARNING] DetectNumberOfLabels() stub called." << std::endl
                  << "           Requires ITK image reading to be implemented." << std::endl
                  << "           Using default: 5 labels (may be incorrect)" << std::endl;

        return 5;  // Placeholder
    }

private:
    /**
     * Helper: Parse comma-separated list of doubles.
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

#endif  // __LabelWeightsParser_h
