/*======================================================================================================*/
/*                                     Beam Power Configurations                                        */	
/*                            (Use this chart to detemine @power in units of pot/yr)                    */
/*  MI prot/pulse    Energy (GeV)    Cycle time    Beam power (MW)    Uptime&efficiency    pot/year     */
/*     7.50E+13          120            1.2           1.20E+00              0.56           1.10E+21     */
/*     7.50E+13           80            0.9           1.07E+00              0.56           1.47E+21     */
/*     7.50E+13           60            0.7           1.03E+00              0.56           1.89E+21     */
/*                                                                                                      */
/*======================================================================================================*/


nuflux(#DplusNuMode)<
@flux_file="./Fluxes_smeft/AllFluxesDUNE/FluxesDplusNuModeND.txt"
@time = 3.5 /* years */
@power = 11.0  /*e20 POT per year ****MATCH FLUX POWER TO TABLE ABOVE**** */
@eft_coeff_file = "./Fluxes_smeft/ProductionCoefficients_DUNE/ProdALLDpDecay.txt"
@norm = 1e13 
@quark_flavors = 1 : 0
>

nuflux(#DplusNuBarMode)<
@flux_file="./Fluxes_smeft/AllFluxesDUNE/FluxesDplusNuModeND.txt"
@time = 3.5 /* years */
@power = 11.0  /*e20 POT per year ****MATCH FLUX POWER TO TABLE ABOVE**** */
@eft_coeff_file = "./Fluxes_smeft/ProductionCoefficients_DUNE/ProdALLDpDecay.txt"
@norm = 1e13
@quark_flavors = 1 : 0
>

nuflux(#DSNuMode)<
@flux_file="./Fluxes_smeft/AllFluxesDUNE/FluxesDSNuModeND.txt"
@time = 3.5 /* years */
@power = 11.0  /*e20 POT per year ****MATCH FLUX POWER TO TABLE ABOVE**** */
@eft_coeff_file = "./Fluxes_smeft/ProductionCoefficients_DUNE/ProdALLDSDecay.txt"
@norm = 1e13
@quark_flavors = 1 : 1
>

nuflux(#DSNuBarMode)<
@flux_file="./Fluxes_smeft/AllFluxesDUNE/FluxesDSNuModeND.txt"
@time = 3.5 /* years */
@power = 11.0  /*e20 POT per year ****MATCH FLUX POWER TO TABLE ABOVE**** */
@eft_coeff_file = "./Fluxes_smeft/ProductionCoefficients_DUNE/ProdALLDSDecay.txt"
@norm = 1e13
@quark_flavors = 1 : 1
>


nuflux(#KaonNuMode)<
@flux_file="./Fluxes_smeft/AllFluxesDUNE/FluxesKaonNuModeND.txt"
@time = 3.5 /* years */
@power = 11.0  /*e20 POT per year ****MATCH FLUX POWER TO TABLE ABOVE**** */
@eft_coeff_file = "./Fluxes_smeft/ProductionCoefficients_DUNE/ProdALLKaonDecayCombinedNuModeND.txt"
@norm = 4e14
@quark_flavors = 0 : 1
>

nuflux(#KaonNuBarMode)<
@flux_file="./Fluxes_smeft/AllFluxesDUNE/FluxesKaonNuBarModeND.txt"
@time = 3.5 /* years */
@power = 11.0  /*e20 POT per year ****MATCH FLUX POWER TO TABLE ABOVE**** */
@eft_coeff_file = "./Fluxes_smeft/ProductionCoefficients_DUNE/ProdALLKaonDecayCombinedNuBarModeND.txt"
@norm = 4e14
@quark_flavors = 0 : 1
>

nuflux(#MuonNuMode)<
@flux_file="./Fluxes_smeft/AllFluxesDUNE/FluxesMuonNuModeND.txt"
@time = 3.5 /* years */
@power = 11.0  /*e20 POT per year ****MATCH FLUX POWER TO TABLE ABOVE**** */
@norm = 4e14
>

nuflux(#MuonNuBarMode)<
@flux_file="./Fluxes_smeft/AllFluxesDUNE/FluxesMuonNuBarModeND.txt"
@time = 3.5 /* years */
@power = 11.0  /*e20 POT per year ****MATCH FLUX POWER TO TABLE ABOVE**** */
@norm = 4e14
>

nuflux(#PionNuMode)<
@flux_file="./Fluxes_smeft/AllFluxesDUNE/FluxesPionNuModeND.txt"
@time = 3.5 /* years */
@power = 11.0  /*e20 POT per year ****MATCH FLUX POWER TO TABLE ABOVE**** */
@eft_coeff_file = "./Fluxes_smeft/ProductionCoefficients_DUNE/ProdALLPionDecay.txt"
@norm = 4e14
@quark_flavors = 0 : 0
>

nuflux(#PionNuBarMode)<
@flux_file="./Fluxes_smeft/AllFluxesDUNE/FluxesPionNuBarModeND.txt"
@time = 3.5 /* years */
@power = 11.0  /*e20 POT per year ****MATCH FLUX POWER TO TABLE ABOVE**** */
@eft_coeff_file = "./ProductionCoefficients_DUNE/ProdALLPionDecay.txt"
@norm = 4e14
@quark_flavors = 0 : 0
>


