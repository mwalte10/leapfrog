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
//' Notes:
//' * Consider adding half the migrants at the start and half at the end?
//' * Unsure whether to have this function accept raw pointers or TensorMap
//'   - If passing Tensor objects, it would be ideal to pass TensorBase, but
//'     this is not possible because dimensions are not known.

template <typename Type, size_t NG, size_t pAG, size_t pIDX_FERT, size_t pAG_FERT>
void leapfrog_sim(const Type *p_basepop,
		  const Type *p_sx,
		  const Type *p_netmigr,
		  const Type *p_asfr,
		  const Type *p_births_sex_prop,
		  //
		  // natural history
		  //
		  //settings
		  const int sim_years,
		  //
		  //outputs
		  Type *p_totpop1) {

  // macros
  // TODO: unsure if these should be defined here or elsewhere. Maybe in a future class.
  
  typedef Eigen::TensorMap<Eigen::Tensor<const Type, 2>> TensorMapXXcT;
  typedef Eigen::TensorMap<Eigen::Tensor<const Type, 3>> TensorMapXXXcT;
  typedef Eigen::TensorMap<Eigen::Tensor<Type, 3>> TensorMapXXXT;

  const size_t MALE = 0;
  const size_t FEMALE = 1;

  // // inputs

  // demography
  const TensorMapXXcT basepop(p_basepop, pAG, NG);
  const TensorMapXXXcT sx(p_sx, pAG+1, NG, sim_years);
  const TensorMapXXXcT netmigr(p_netmigr, pAG, NG, sim_years);
  const TensorMapXXcT asfr(p_asfr, pAG_FERT, sim_years);
  const TensorMapXXcT births_sex_prop(p_births_sex_prop, NG, sim_years);

  // HIV natural history

  // const TensorMapXXXcT cd4_initdist(REAL(
  //   multi_array_ref<double, 3> cd4_initdist(REAL(getListElement(s_fp, "cd4_initdist")), extents[NG][hAG][hDS]);
  //   multi_array_ref<double, 3> cd4_prog(REAL(getListElement(s_fp, "cd4_prog")), extents[NG][hAG][hDS-1]);
  //   multi_array_ref<double, 3> cd4_mort(REAL(getListElement(s_fp, "cd4_mort")), extents[NG][hAG][hDS]);
  //   multi_array_ref<double, 4> art_mort(REAL(getListElement(s_fp, "art_mort")), extents[NG][hAG][hDS][hTS]);
  //   multi_array_ref<double, 2> artmx_timerr(REAL(getListElement(s_fp, "artmx_timerr")), extents[PROJ_YEARS][hTS]);

  
  // outputs
  TensorMapXXXT totpop1(p_totpop1, pAG, NG, sim_years);
  Tensor<Type, 3> natdeaths(pAG, NG, sim_years);

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
      TensorFixedSize<Type, Sizes<pAG>> migrate_a;
      for(int a = 1; a < pAG; a++) {
	migrate_a(a) = netmigr(a, g, t) * (1.0 + sx(a, g, t)) * 0.5 / totpop1(a, g, t);
	totpop1(a, g, t) *= 1.0 + migrate_a(a);
      }
    }

    // fertility

    Type births = 0.0;
    for(int af = 0; af < pAG_FERT; af++) {
      births += (totpop1(pIDX_FERT + af, FEMALE, t-1) + totpop1(pIDX_FERT + af, FEMALE, t)) * 0.5 * asfr(af, t);
    }

    // add births
    for(int g = 0; g < NG; g++) {
      totpop1(0, g, t) = births * births_sex_prop(g, t) * sx(0, g, t);

      // Assume 2/3 survival rate since mortality in first six months higher than
      // second 6 months (Spectrum manual, section 6.2.7.4)
      Type migrate_a0 = netmigr(0, g, t) * (1.0 + 2.0 * sx(0, g, t)) / 3.0 / totpop1(0, g, t);
      totpop1(0, g, t) *= 1.0 + migrate_a0;
    }

    // demographic projection of the HIV population
    
      
    //   }
    // 	pop(pAG-1, g, t) += pop_t(pAG-1, g, t-1); // open age group
    //   }
    
    // TensorFixedSize<double, Sizes<hAG, NG>> hiv_ag_prob;
    //   // Tensor<double, 2> hiv_ag_prob(hAG, NG);      
    //   hiv_ag_prob.setZero();
      
    //   for(int g = 0; g < NG; g++){
    //     int a = 0;
    //     for(int ha = 0; ha < (hAG-1); ha++){
    //       for(int i = 0; i < hAG_SPAN[ha]; i++){
    //         hiv_ag_prob(ha, g) += pop_t(a, g, HIVP, t-1);
    //         a++;
    //       }
    //       hiv_ag_prob(ha, g) = (hiv_ag_prob(ha, g) > 0) ? pop_t(a-1, g, HIVP, t-1) / hiv_ag_prob(ha, g) : 0;
    //     }
    // 	// Note: loop stops at hAG-1; no one ages out of the open-ended age group
    //   }

    //   for(int g = 0; g < NG; g++)
    //     for(int ha = 1; ha < hAG; ha++)
    //       for(int hm = 0; hm < hDS; hm++){
    //         hivpop_t(hm, ha, g, t) = (1-hiv_ag_prob(ha, g)) * hivpop_t(hm, ha, g, t-1) + hiv_ag_prob(ha-1, g)*hivpop_t(hm, ha-1, g, t-1);
    //         if(t > t_ART_start)
    //           for(int hu = 0; hu < hTS; hu++)
    //             artpop_t(hu, hm, ha, g, t) = (1-hiv_ag_prob(ha, g)) * artpop_t(hu, hm, ha, g, t-1) + hiv_ag_prob(ha-1, g)*artpop_t(hu, hm, ha-1, g, t-1);
    //       }


    //   // non-HIV mortality and netmigration
    //   for(int g = 0; g < NG; g++){
    //     int a = 0;
    //     for(int ha = 0; ha < hAG; ha++){
    //       double deathsmig_ha = 0, hivpop_ha = 0;
    //       for(int i = 0; i < hAG_SPAN[ha]; i++){

    //         hivpop_ha += pop_t(a, g, HIVP, t);

    //         // non-HIV mortality
    //         double qx = 1.0 - Sx(a, g, t);
    //         double ndeaths_a = pop_t(a, g, HIVN, t) * qx;
    //         pop_t(a, g, HIVN, t) -= ndeaths_a; // survival HIV- population
    //         double hdeaths_a = pop_t(a, g, HIVP, t) * qx;
    //         deathsmig_ha -= hdeaths_a;
    //         pop_t(a, g, HIVP, t) -= hdeaths_a;   // survival HIV+ population
    //         natdeaths(a, g, t) = ndeaths_a + hdeaths_a;

    //         // net migration
    //         double migrate_a = netmigr(a, g, t) * (1+Sx(a, g, t))/2.0 / (pop_t(a, g, HIVN, t) + pop_t(a, g, HIVP, t));
    //         pop_t(a, g, HIVN, t) *= 1+migrate_a;
    //         double hmig_a = migrate_a * pop_t(a, g, HIVP, t);
    //         deathsmig_ha += hmig_a;
    //         pop_t(a, g, HIVP, t) += hmig_a;

    //         a++;
    //       }

    //       // migration and deaths for hivpop
    //       double deathmigrate_ha = hivpop_ha > 0 ? deathsmig_ha / hivpop_ha : 0.0;
    //       for(int hm = 0; hm < hDS; hm++){
    //         hivpop_t(hm, ha, g, t) *= 1+deathmigrate_ha;
    //         if(t > t_ART_start)
    //           for(int hu = 0; hu < hTS; hu++)
    //             artpop_t(hu, hm, ha, g , t) *= 1+deathmigrate_ha;
    //       } // loop over hm
    //     } // loop over ha
    //   } // loop over g


    //   // fertility
    //   double births = 0.0, births_by_ha[hAG_FERT];
    //   memset(births_by_ha, 0, hAG_FERT*sizeof(double));
    //   for(int m = 0; m < pDS; m++){
    //     int a = pIDX_FERT;
    //     for(int ha = hIDX_FERT; ha < hIDX_FERT+hAG_FERT; ha++){
    //       for(int i = 0; i < hAG_SPAN[ha]; i++){
    //         births_by_ha[ha-hIDX_FERT] += (pop_t(a, FEMALE, m, t-1) + pop_t(a, FEMALE, m, t))/2 * asfr(a, t);
    //         a++;
    //       }
    //     }
    //   }
    //   for(int ha = hIDX_FERT; ha < hAG_FERT; ha++)
    //     births += births_by_ha[ha-hIDX_FERT];

    //   if(t + AGE_START < PROJ_YEARS)
    //     for(int g = 0; g < NG; g++)
    //       birthslag(g, t + AGE_START-1) = srb(g, t) * births;

  }
  
  return;
}

	     




#endif // LEAPFROG_H
