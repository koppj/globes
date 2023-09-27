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
@flux_file="./Fluxes_smeft/AllFluxesDUNE/FluxesDplusNuModeFD.txt"
@time = 3.5 /* years */
@power = 11.0  /*e20 POT per year ****MATCH FLUX POWER TO TABLE ABOVE**** */
@eft_coeff_file = "./Fluxes_smeft/ProductionCoefficients_DUNE_separated/plus_ProdALLDpDecay.txt"
@norm = 8.956972e+15
@quark_flavors = 1 : 0
>

nuflux(#DplusNuBarMode)<
@flux_file="./Fluxes_smeft/AllFluxesDUNE/FluxesDplusNuModeFD.txt"
@time = 3.5 /* years */
@power = 11.0  /*e20 POT per year ****MATCH FLUX POWER TO TABLE ABOVE**** */
@eft_coeff_file = "./Fluxes_smeft/ProductionCoefficients_DUNE_separated/plus_ProdALLDpDecay.txt"
@norm = 8.956972e+15
@quark_flavors = 1 : 0
>

nuflux(#DSNuMode)<
@flux_file="./Fluxes_smeft/AllFluxesDUNE/FluxesDSNuModeFD.txt"
@time = 3.5 /* years */
@power = 11.0  /*e20 POT per year ****MATCH FLUX POWER TO TABLE ABOVE**** */
@eft_coeff_file = "./Fluxes_smeft/ProductionCoefficients_DUNE_separated/plus_ProdALLDSDecay.txt"
@norm = 8.956972e+15
@quark_flavors = 1 : 1
>

nuflux(#DSNuBarMode)<
@flux_file="./Fluxes_smeft/AllFluxesDUNE/FluxesDSNuModeFD.txt"
@time = 3.5 /* years */
@power = 11.0  /*e20 POT per year ****MATCH FLUX POWER TO TABLE ABOVE**** */
@eft_coeff_file = "./Fluxes_smeft/ProductionCoefficients_DUNE_separated/plus_ProdALLDSDecay.txt"
@norm = 8.956972e+15
@quark_flavors = 1 : 1
>


nuflux(#KaonNuMode)<
@flux_file="./Fluxes_smeft/AllFluxesDUNE/FluxesKaonNuModeFD.txt"
@time = 3.5 /* years */
@power = 11.0  /*e20 POT per year ****MATCH FLUX POWER TO TABLE ABOVE**** */
@eft_coeff_file = "./Fluxes_smeft/ProductionCoefficients_DUNE_separated/plus_ProdALLKaonDecayCombinedNuModeFD.txt"
@norm =8.956972e+15
@quark_flavors = 0 : 1
>

nuflux(#KaonNuBarMode)<
@flux_file="./Fluxes_smeft/AllFluxesDUNE/FluxesKaonNuBarModeFD.txt"
@time = 3.5 /* years */
@power = 11.0  /*e20 POT per year ****MATCH FLUX POWER TO TABLE ABOVE**** */
@eft_coeff_file = "./Fluxes_smeft/ProductionCoefficients_DUNE_separated/plus_ProdALLKaonDecayCombinedNuBarModeFD.txt"
@norm =8.956972e+15
@quark_flavors = 0 : 1
>

nuflux(#MuonNuMode)<
@flux_file="./Fluxes_smeft/AllFluxesDUNE/FluxesMuonNuModeFD.txt"
@time = 3.5 /* years */
@power = 11.0  /*e20 POT per year ****MATCH FLUX POWER TO TABLE ABOVE**** */
@norm = 8.956972e+15
>

nuflux(#MuonNuBarMode)<
@flux_file="./Fluxes_smeft/AllFluxesDUNE/FluxesMuonNuBarModeFD.txt"
@time = 3.5 /* years */
@power = 11.0  /*e20 POT per year ****MATCH FLUX POWER TO TABLE ABOVE**** */
@norm =8.956972e+15
>

nuflux(#PionNuMode)<
@flux_file="./Fluxes_smeft/AllFluxesDUNE/FluxesPionNuModeFD.txt"
@time = 3.5 /* years */
@power = 11.0  /*e20 POT per year ****MATCH FLUX POWER TO TABLE ABOVE**** */
@eft_coeff_file = "./Fluxes_smeft/ProductionCoefficients_DUNE_separated/plus_ProdALLPionDecay.txt"
@norm = 8.956972e+15
@quark_flavors = 0 : 0
>

nuflux(#PionNuBarMode)<
@flux_file="./Fluxes_smeft/AllFluxesDUNE/FluxesPionNuBarModeFD.txt"
@time = 3.5 /* years */
@power = 11.0  /*e20 POT per year ****MATCH FLUX POWER TO TABLE ABOVE**** */
@eft_coeff_file = "./Fluxes_smeft/ProductionCoefficients_DUNE_separated/plus_ProdALLPionDecay.txt"
@norm = 8.956972e+15
@quark_flavors = 0 : 0
>

nuflux(#DplusNuMode_minus)<
@flux_file="./Fluxes_smeft/AllFluxesDUNE/FluxesDplusNuModeFD.txt"
@time = 3.5 /* years */
@power = 11.0  /*e20 POT per year ****MATCH FLUX POWER TO TABLE ABOVE**** */
@eft_coeff_file = "./Fluxes_smeft/ProductionCoefficients_DUNE_separated/minus_ProdALLDpDecay.txt"
@norm = 8.956972e+15
@quark_flavors = 1 : 0
>

nuflux(#DplusNuBarMode_minus)<
@flux_file="./Fluxes_smeft/AllFluxesDUNE/FluxesDplusNuModeFD.txt"
@time = 3.5 /* years */
@power = 11.0  /*e20 POT per year ****MATCH FLUX POWER TO TABLE ABOVE**** */
@eft_coeff_file = "./Fluxes_smeft/ProductionCoefficients_DUNE_separated/minus_ProdALLDpDecay.txt"
@norm = 8.956972e+15
@quark_flavors = 1 : 0
>

nuflux(#DSNuMode_minus)<
@flux_file="./Fluxes_smeft/AllFluxesDUNE/FluxesDSNuModeFD.txt"
@time = 3.5 /* years */
@power = 11.0  /*e20 POT per year ****MATCH FLUX POWER TO TABLE ABOVE**** */
@eft_coeff_file = "./Fluxes_smeft/ProductionCoefficients_DUNE_separated/minus_ProdALLDSDecay.txt"
@norm = 8.956972e+15
@quark_flavors = 1 : 1
>

nuflux(#DSNuBarMode_minus)<
@flux_file="./Fluxes_smeft/AllFluxesDUNE/FluxesDSNuModeFD.txt"
@time = 3.5 /* years */
@power = 11.0  /*e20 POT per year ****MATCH FLUX POWER TO TABLE ABOVE**** */
@eft_coeff_file = "./Fluxes_smeft/ProductionCoefficients_DUNE_separated/minus_ProdALLDSDecay.txt"
@norm =8.956972e+15
@quark_flavors = 1 : 1
>


nuflux(#KaonNuMode_minus)<
@flux_file="./Fluxes_smeft/AllFluxesDUNE/FluxesKaonNuModeFD.txt"
@time = 3.5 /* years */
@power = 11.0  /*e20 POT per year ****MATCH FLUX POWER TO TABLE ABOVE**** */
@eft_coeff_file = "./Fluxes_smeft/ProductionCoefficients_DUNE_separated/minus_ProdALLKaonDecayCombinedNuModeND.txt"
@norm = 8.956972e+15
@quark_flavors = 0 : 1
>

nuflux(#KaonNuBarMode_minus)<
@flux_file="./Fluxes_smeft/AllFluxesDUNE/FluxesKaonNuBarModeFD.txt"
@time = 3.5 /* years */
@power = 11.0  /*e20 POT per year ****MATCH FLUX POWER TO TABLE ABOVE**** */
@eft_coeff_file = "./Fluxes_smeft/ProductionCoefficients_DUNE_separated/minus_ProdALLKaonDecayCombinedNuBarModeFD.txt"
@norm = 8.956972e+15
@quark_flavors = 0 : 1
>

nuflux(#MuonNuMode_minus)<
@flux_file="./Fluxes_smeft/AllFluxesDUNE/FluxesMuonNuModeFD.txt"
@time = 3.5 /* years */
@power = 11.0  /*e20 POT per year ****MATCH FLUX POWER TO TABLE ABOVE**** */
@norm = 8.956972e+15
>

nuflux(#MuonNuBarMode_minus)<
@flux_file="./Fluxes_smeft/AllFluxesDUNE/FluxesMuonNuBarModeFD.txt"
@time = 3.5 /* years */
@power = 11.0  /*e20 POT per year ****MATCH FLUX POWER TO TABLE ABOVE**** */
@norm = 8.956972e+15
>

nuflux(#PionNuMode_minus)<
@flux_file="./Fluxes_smeft/AllFluxesDUNE/FluxesPionNuModeFD.txt"
@time = 3.5 /* years */
@power = 11.0  /*e20 POT per year ****MATCH FLUX POWER TO TABLE ABOVE**** */
@eft_coeff_file = "./Fluxes_smeft/ProductionCoefficients_DUNE_separated/minus_ProdALLPionDecay.txt"
@norm =  8.956972e+15
@quark_flavors = 0 : 0
>

nuflux(#PionNuBarMode_minus)<
@flux_file="./Fluxes_smeft/AllFluxesDUNE/FluxesPionNuBarModeFD.txt"
@time = 3.5 /* years */
@power = 11.0  /*e20 POT per year ****MATCH FLUX POWER TO TABLE ABOVE**** */
@eft_coeff_file = "./Fluxes_smeft/ProductionCoefficients_DUNE_separated/minus_ProdALLPionDecay.txt"
@norm =  8.956972e+15
@quark_flavors = 0 : 0
>

