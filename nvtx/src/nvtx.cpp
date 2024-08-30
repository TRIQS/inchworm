#include "nvtx.hpp"

#ifdef USE_NVTX
NVTXTracer::NVTXTracer(const char *name, int color_id) {
  const uint32_t colors[] = {
     0xff00ff00, // Green
     0xff0000ff, // Blue
     0xffffff00, // Yellow
     0xffff00ff, // Magenta
     0xff00ffff, // Cyan
     0xffff0000, // Red
     0xffffffff, // White
     0xff808080, // Gray
     0xff800000, // Maroon
     0xff008000, // Dark Green
     0xff000080, // Navy
     0xff800080, // Purple
     0xff808000, // Olive
     0xff008080, // Teal
     0xffc0c0c0  // Silver
  };
  const int num_colors = sizeof(colors) / sizeof(uint32_t);

  color_id = color_id % num_colors;

  nvtxEventAttributes_t eventAttrib = {0};
  eventAttrib.version               = NVTX_VERSION;
  eventAttrib.size                  = NVTX_EVENT_ATTRIB_STRUCT_SIZE;
  eventAttrib.colorType             = NVTX_COLOR_ARGB;
  eventAttrib.color                 = colors[color_id];
  eventAttrib.messageType           = NVTX_MESSAGE_TYPE_ASCII;
  eventAttrib.message.ascii         = name;
  nvtxRangePushEx(&eventAttrib);
}

NVTXTracer::~NVTXTracer() { nvtxRangePop(); }
#endif