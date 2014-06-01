namespace principia {
namespace integrators {

template<typename Position, typename Momentum>
SymplecticIntegrator::Parameters<Position, Momentum>::Parameters()
    : p_error(nullptr),
      q_error(nullptr),
      t_error(0 * Time::SIUnit()) {}

}  // namespace integrators
}  // namespace principia
