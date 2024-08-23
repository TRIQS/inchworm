#ifdef USE_NVTX
#include "nvtx3/nvToolsExt.h"
#endif

class NVTXTracer {
public:
    NVTXTracer(const char* name, int color_id = 0);
    ~NVTXTracer();
};

#ifdef USE_NVTX
#define NVTX_RANGE(name, cid) NVTXTracer uniq_name_using_macros(name, cid);
#else
#define NVTX_RANGE(name, cid) // Empty definition
#endif