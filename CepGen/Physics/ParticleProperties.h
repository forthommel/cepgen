#ifndef CepGen_Physics_ParticleProperties_h
#define CepGen_Physics_ParticleProperties_h

#include <string>

namespace cepgen
{
  /// A collection of physics constants associated to a single particle
  struct ParticleProperties
  {
    std::string name; ///< Particle name
    std::string description; ///< Human-readable name
    short colours; ///< Colour factor
    double mass; ///< Mass, in GeV/c\f$^2\f$
    double width; ///< Decay width, in GeV/c\f$^2\f$
    short charge; ///< Electric charge, in \f$e\f$/3
    bool fermion; ///< Is the particle a fermion?
  };
}

#endif
