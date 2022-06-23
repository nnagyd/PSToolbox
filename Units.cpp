#include "Units.h"
//! Unit conversion constants
/*! Class filled up with conversion coefficients between various Imperial
	and SI measurements. */
Units::Units() {
	ft_to_m = 304.8 / 1000.;
	m_to_ft = 1 / ft_to_m;
	inch_to_m = 2.54 / 100.;
	m_to_inch = 1. / inch_to_m;
	lbm_to_kg = 0.45359237;
	kg_to_lbm = 1. / lbm_to_kg;
	lbf_to_N = 4.44822;
	N_to_lbf = 1. / lbf_to_N;
	psi_to_Pa = 6894.75;
	Pa_to_psi = 1. / psi_to_Pa;
	psi_to_bar = psi_to_Pa / 1e5;
	bar_to_psi = 1. / psi_to_bar;
	inch3_to_m3 = inch_to_m * inch_to_m * inch_to_m;
	m3_to_inch3 = 1. / inch3_to_m3;
}

Units::~Units() {}
