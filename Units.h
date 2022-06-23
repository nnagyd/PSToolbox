#pragma once

//! Class to contain unit conversion values
/*! A class that contains the unit conversion values used in the various simulations
    Serves as a base class for the other classes s it easily provides conversion values that way.
  \sa LWP, SCP, Valve, Reservoir
*/
class Units
{
public:
  //! Constructor
  /*! Create an instance of the Units class. Values of the conversion
      are assigned during this phase. (and as such are not constants.)
  */
  Units();
  //! Destuctor
  ~Units();
  double ft_to_m; //!< Converts foot to meters
  double m_to_ft; //!< Convert meters to foot
  double inch_to_m; //!< Convert inches to meters
  double m_to_inch; //!< Convert metres to inches
  double lbm_to_kg; //!< Convert pounds to kilograms
  double kg_to_lbm; //!< Multiplier for conversion from kilogram to pound
  double lbf_to_N; //!< Multiplier for conversion of poundforce to Newtons
  double N_to_lbf; //!< Mulitplier for conversion of Newtons to poundforce
  double psi_to_Pa; //!<Mulitplier ofr conversion of pound per square inch to Pa
  double Pa_to_psi; //!< Multiplier for conversion of pascals to pound per square inch
  double psi_to_bar; //!< Mulitplier for conversion of psi to bar
  double bar_to_psi; //!< Multiplier for conversion of bar to psi
  double m3_to_inch3; //!< Multiplier for conversion of cubic meter to cumbic inch
  double inch3_to_m3; //!< Multiplier for conversion of cubic inches to cubic metres
};
