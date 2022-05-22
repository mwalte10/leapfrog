#ifndef LEAPFROG_H
#define LEAPFROG_H

#include <unsupported/Eigen/CXX11/Tensor>

using Eigen::TensorBase;
using Eigen::TensorMap;
using Eigen::Tensor;
using Eigen::TensorFixedSize;
using Eigen::Sizes;

//' The output memory is passed as an argument rather than constructed
//' by the function.
//'
//' Template parameters
//' @param Type data type variables involved in calculations (typically 'double'
//'   but required to be a templated parameter for autodiff integration).
//' @param NG number of genders (= 2)
//' @param pAG number of population age groups (e.g. 81 for ages 0 to 80+)
//' @param pIDX_FERT first index eligible for fertility
//' @param pAG_FERT number of ages eligible for fertility
//'
//' @details
//' State space dimensions are specified as template parameters so that
//' these are specified at compile time, allowing stack allocation of
//' working arrays.
//'

template <typename Type, int NG, int pAG, int pIDX_FERT, int pAG_FERT>
  void leapfrog_sim(const Type *p_basepop,
                    const Type *p_sx,
                    const Type *p_netmigr,
                    const Type *p_asfr,
                    const Type *p_births_sex_prop,
                    //
                    //settings
                    const int sim_years,
		    //
                    //outputs
                    Type *p_totpop1,
		    Type *p_births,
                    Type *p_natdeaths) {
  
  // macros
  // TODO: unsure if these should be defined here or elsewhere. Maybe in a future class.

  typedef Eigen::TensorMap<Eigen::Tensor<const Type, 2>> TensorMapX2cT;
  typedef Eigen::TensorMap<Eigen::Tensor<const Type, 3>> TensorMapX3cT;

  typedef Eigen::TensorMap<Eigen::Tensor<Type, 1>> TensorMapX1T;
  typedef Eigen::TensorMap<Eigen::Tensor<Type, 3>> TensorMapX3T;

  const int MALE = 0;
  const int FEMALE = 1;

  // // inputs

  // demography
  const TensorMapX2cT basepop(p_basepop, pAG, NG);
  const TensorMapX3cT sx(p_sx, pAG+1, NG, sim_years);
  const TensorMapX3cT netmigr(p_netmigr, pAG, NG, sim_years);
  const TensorMapX2cT asfr(p_asfr, pAG_FERT, sim_years);
  const TensorMapX2cT births_sex_prop(p_births_sex_prop, NG, sim_years);

  // outputs
  TensorMapX3T totpop1(p_totpop1, pAG, NG, sim_years);
  TensorMapX1T births(p_births, sim_years);  
  TensorMapX3T natdeaths(p_natdeaths, pAG, NG, sim_years);

  // initialise population

  for(int g = 0; g < NG; g++) {
    for(int a = 0; a < pAG; a++) {
      totpop1(a, g, 0) = basepop(a, g);
    }
  }

  ////////////////////////////////////
  ////  do population projection  ////
  ////////////////////////////////////

  for(int t = 1; t < sim_years; t++){

    TensorFixedSize<Type, Sizes<pAG, NG>> migrate_ag;

    // ageing and non-HIV mortality
    for(int g = 0; g < NG; g++){

      for(int a = 1; a < pAG; a++) {
        natdeaths(a, g, t) = totpop1(a-1, g, t-1) * (1.0 - sx(a, g, t));
        totpop1(a, g, t) = totpop1(a-1, g, t-1) - natdeaths(a, g, t);
      }

      // open age group
      Type natdeaths_open_age = totpop1(pAG-1, g, t-1) * (1.0 - sx(pAG, g, t));
      natdeaths(pAG-1, g, t) += natdeaths_open_age;
      totpop1(pAG-1, g, t) += totpop1(pAG-1, g, t-1) - natdeaths_open_age;

      // net migration
      for(int a = 1; a < pAG; a++) {
        migrate_ag(a, g) = netmigr(a, g, t) * (1.0 + sx(a, g, t)) * 0.5 / totpop1(a, g, t);
        totpop1(a, g, t) *= 1.0 + migrate_ag(a, g);
      }
    }

    // fertility

    births(t) = 0.0;
    for(int af = 0; af < pAG_FERT; af++) {
      births(t) += (totpop1(pIDX_FERT + af, FEMALE, t-1) + totpop1(pIDX_FERT + af, FEMALE, t)) * 0.5 * asfr(af, t);
    }

    // add births
    for(int g = 0; g < NG; g++) {
      Type births_sex = births(t) * births_sex_prop(g, t);
      natdeaths(0, g, t) = births_sex * (1.0 - sx(0, g, t));
      totpop1(0, g, t) =  births_sex * sx(0, g, t);

      // Assume 2/3 survival rate since mortality in first six months higher than
      // second 6 months (Spectrum manual, section 6.2.7.4)
      Type migrate_a0 = netmigr(0, g, t) * (1.0 + 2.0 * sx(0, g, t)) / 3.0 / totpop1(0, g, t);
      totpop1(0, g, t) *= 1.0 + migrate_a0;
    }

  }

  return;
}


#endif // LEAPFROG_H
