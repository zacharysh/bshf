#include "wavefunction.hpp"

Wavefunction::Wavefunction(int n_, int l_, int m_)
: n(n_), l(l_), m(m_)
{
    std::cout << "> Generating wavefunction for atomic state (n, l, m) = (" << n_ << ", " << l_ << ", " << m_ << ").\n";
}