// Include the defined classes that are to be exported to python
#include "BondFlipUpdater.h"

#include <hoomd/extern/pybind/include/pybind11/pybind11.h>

// specify the python module. Note that the name must expliclty match the PROJECT() name provided in CMakeLists
// (with an underscore in front)
PYBIND11_PLUGIN(_bondflip_plugin)
    {
    pybind11::module m("_bondflip_plugin");
    export_BondFlipUpdater(m);

    #ifdef ENABLE_CUDA
    export_BondFlipUpdaterGPU(m);
    #endif

    return m.ptr();
    }
