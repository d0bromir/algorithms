// Used when CERTA is built without CUDA.
#include <stdexcept>

#include "mapper.h"

namespace certa {

struct GpuMapper::Impl {};

GpuMapper::GpuMapper(const Reference&, const Index&, int) {
  throw std::runtime_error("this build has no GPU support; rebuild with -DCERTA_CUDA=ON");
}
GpuMapper::~GpuMapper() = default;
void GpuMapper::map(const Params&, const ReadBatch&, std::vector<Result>&) {}
std::string GpuMapper::device_name() const { return ""; }
bool GpuMapper::compiled_in() { return false; }

}  // namespace certa
