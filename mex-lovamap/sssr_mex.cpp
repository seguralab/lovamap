#include "mex.hpp"
#include "mexAdapter.hpp"

#include "SSSR.h"
#include <sstream>
#include <vector>


using namespace matlab::data;
using namespace vspc;

class MexFunction : public matlab::mex::Function
{
public:
    using MatlabPtr = std::shared_ptr<matlab::engine::MATLABEngine>;

    /// Matlab-facing interface.
    void operator()(matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs);

private:
    /// Helper function for checking valid input and output arguments from Matlab.
    void checkArgs(matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs);

    /// @brief Helper function for converting all cycle data into a 2D array.
    /// @details Each row represents the nodes in a cycle, and rows are padded with 0s on
    ///          the right in order for the data to fit into a rectangular array.
    TypedArray<int> outputAsArray(const std::vector<UndirectedGraph>& cycles);

    /// Helper function to print things to the Matlab console.
    void printToConsole(std::ostringstream& oss);

    MatlabPtr           mEngine = getEngine();
    ArrayFactory        mFactory;
    std::ostringstream  mOSS;
};

// -----------------------------------------------------------------------------

void
MexFunction::operator()(
    matlab::mex::ArgumentList outputs,
    matlab::mex::ArgumentList inputs)
{
    checkArgs(outputs, inputs);
    const TypedArray<int>& inEdges = inputs[0];
    const int maxCycleLength = inputs[1][0];

    // Initialize SSSR settings
    SSSR::Settings settings;
    settings.setMaxCycleLength(maxCycleLength);

    // Initialize edge data
    std::vector<UndirectedGraph::Edge> edgeData;

    // Extract edge data from input
    ArrayDimensions dim = inEdges.getDimensions();
    for (size_t i = 0; i < dim[0]; ++i) {
        edgeData.emplace_back(inEdges[i][0], inEdges[i][1]);
    }

    // Initialize SSSR functor and solve
    SSSR op(UndirectedGraph{edgeData}, settings);
    const std::vector<UndirectedGraph>& cycles = op.run();

    outputs[0] = cycles.size() > 0 ? outputAsArray(cycles) : mFactory.createEmptyArray();
}

void
MexFunction::checkArgs(
    matlab::mex::ArgumentList outputs,
    matlab::mex::ArgumentList inputs)
{
    if (inputs.size() != 2)
    {
        mEngine->feval(u"error", 0,
            std::vector<matlab::data::Array>(
                {
                    mFactory.createScalar(
                        "Error using sssr\n"
                        "Two input arguments required: edge data of an undirected graph,"
                        " and the maximum allowable length of loops.")
                }));
    }

    if (outputs.size() > 1)
    {
        mEngine->feval(u"error", 0,
            std::vector<matlab::data::Array>(
                {
                    mFactory.createScalar(
                        "Error using sssr\n"
                        "At most one output argument allowed.")
                }));
    }

    if (inputs[0].getType() != matlab::data::ArrayType::INT32 ||
        inputs[0].getDimensions()[1] != 2)
    {
        mEngine->feval(u"error", 0,
            std::vector<matlab::data::Array>(
                {
                    mFactory.createScalar(
                            "Error using sssr\n"
                            "First argument must be an N x 2 matrix of type int32.")
                }));
    }

    if (inputs[1].getNumberOfElements() != 1)
    {
        mEngine->feval(u"error", 0,
            std::vector<matlab::data::Array>(
                {
                    mFactory.createScalar(
                            "Error using sssr\n"
                            "Second argument must be a scalar.")
                }));
    }
}

TypedArray<int>
MexFunction::outputAsArray(const std::vector<UndirectedGraph>& cycles)
{
    const size_t nCycles = cycles.size();
    const size_t maxLength = convertGraphToPath(cycles.back()).length();
    const ArrayDimensions outDim{nCycles, maxLength + 1};

    TypedArray<int> result = mFactory.createArray<int>(outDim);

    for (size_t i = 0; i < nCycles; ++i) {
        SSSR::Path path = convertGraphToPath(cycles[i]);
        size_t j = 0;
        for (const auto& node : path) {
            result[i][j++] = node;
        }
    }
    return result;
}

void
MexFunction::printToConsole(std::ostringstream& oss)
{
    mEngine->feval(u"fprintf", 0,
        std::vector<matlab::data::Array>({ mFactory.createScalar(oss.str()) }));
    oss.str("");
    oss.clear();
}
