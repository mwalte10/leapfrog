# Spectrum test files

These are a series of test files constructed to 




Files were constructed using default data for Botswana
in Spectrum v5.13 on 12 February 2022.


## `bwa_aim-adult-no

This file includes 

ART is removed from the projecting

Paediatric HIV is removed by setting mother-to-child transmission probabilities to zero.

* Change all Adult ART percentages to zero.
* Change Adult ART "Percent lost to follow-up each year" to zero.
* Change all Child ART percentages to zero.
* Change child Percent receiving cotrimoxazole to zero.
* Set [Advanced Options] -> [MTCT transition probabilities] all equal zero.


## `v6.13/bwa_aim-adult-art-no-special-elig_v6.13_2022-04-18.PJNZ`

Paediatric HIV is removed by setting mother-to-child transmission probabilities to zero.

All ART inputs are percentages in the demo data file. Adult ART manually converted from 
percentages to counts for years 2000 to 2020 for consistency with actual usage.

The default file had no 'special population' ART eligibility. This is retained.

* Set [Advanced Options] -> [MTCT transition probabilities] all equal zero.
* Manually calculate number on ART 15+ and input for years 2000 to 2020.
