
/* NUE APP */

//D

channel(#NuMode_Dplus_app_osc_nue)<
        @channel =      #DplusNuMode:     +:      m:      e:      #CC:    #app_nue_sig
       @post_smearing_efficiencies = copy(%post_app_FHC_nue_sig)
>

channel(#NuMode_Dplus_app_osc_nuebar)<
         @channel =      #DplusNuMode_minus:     -:      m:      e:      #CC:    #app_nuebar_sig
         @post_smearing_efficiencies = copy(%post_app_FHC_nuebar_sig)
>

channel(#NuMode_Dplus_app_bkg_nue)<
        @channel =      #DplusNuMode:     +:      e:      e:      #CC:    #app_nue_bkg
       @post_smearing_efficiencies = copy(%post_app_FHC_nue_bkg)
>

channel(#NuMode_Dplus_app_bkg_nuebar)<
         @channel =      #DplusNuMode_minus:     -:      e:      e:      #CC:    #app_nuebar_bkg
         @post_smearing_efficiencies = copy(%post_app_FHC_nuebar_bkg)
>

channel(#NuMode_Dplus_app_bkg_numu)<
        @channel =      #DplusNuMode:     +:      m:      m:      #CC:    #app_numu_bkg
       @post_smearing_efficiencies = copy(%post_app_FHC_numu_bkg)
>

channel(#NuMode_Dplus_app_bkg_numubar)<
         @channel =      #DplusNuMode_minus:     -:      m:      m:      #CC:    #app_numubar_bkg
         @post_smearing_efficiencies = copy(%post_app_FHC_numubar_bkg)
>


channel(#NuMode_Dplus_app_bkg_e_nutau)<
        @channel =      #DplusNuMode:     +:      e:      t:      #CCtau:    #app_nutau_bkg
       @post_smearing_efficiencies = copy(%post_app_FHC_nutau_bkg)
>

channel(#NuMode_Dplus_app_bkg_e_nutaubar)<
         @channel =      #DplusNuMode_minus:     -:      e:      t:      #CCtau:    #app_nutaubar_bkg
         @post_smearing_efficiencies = copy(%post_app_FHC_nutaubar_bkg)
>



channel(#NuMode_Dplus_app_bkg_nuNC)<
        @channel =      #DplusNuMode:     +:      m: NOSC_m:      #NC:    #app_NC_bkg
        @post_smearing_efficiencies = copy(%post_app_FHC_NC_bkg)
>

channel(#NuMode_Dplus_app_bkg_nubarNC)<
        @channel =      #DplusNuMode_minus:     -:      m: NOSC_m:      #NC:    #app_NC_bkg
        @post_smearing_efficiencies = copy(%post_app_FHC_aNC_bkg)
>

channel(#NuMode_Dplus_app_bkg_nuNC_e)<
        @channel =      #DplusNuMode:     +:      e: NOSC_e:      #NC:    #app_NC_bkg
        @post_smearing_efficiencies = copy(%post_app_FHC_NC_bkg)
>

channel(#NuMode_Dplus_app_bkg_nubarNC_e)<
        @channel =      #DplusNuMode_minus:     -:      e: NOSC_e:      #NC:    #app_NC_bkg
        @post_smearing_efficiencies = copy(%post_app_FHC_aNC_bkg)
>

channel(#NuMode_Dplus_app_bkg_nuNC_t)<
        @channel =      #DplusNuMode:     +:      t: NOSC_t:      #NC:    #app_NC_bkg
        @post_smearing_efficiencies = copy(%post_app_FHC_NC_bkg)
>

channel(#NuMode_Dplus_app_bkg_nubarNC_t)<
        @channel =      #DplusNuMode_minus:     -:      t: NOSC_t:      #NC:    #app_NC_bkg
        @post_smearing_efficiencies = copy(%post_app_FHC_aNC_bkg)
>

channel(#NuMode_Dplus_app_tau_e)<
        @channel =      #DplusNuMode:     +:      t:      e:      #CC:    #app_nue_bkg
       @post_smearing_efficiencies = copy(%post_app_FHC_nue_bkg)
>

channel(#NuMode_Dplus_app_tau_ebar)<
        @channel =      #DplusNuMode_minus:     -:      t:      e:      #CC:    #app_nue_bkg
       @post_smearing_efficiencies = copy(%post_app_FHC_nue_bkg)
>

channel(#NuMode_Dplus_app_tau_numu)<
        @channel =      #DplusNuMode:     +:      t:      m:      #CC:    #app_numu_bkg
       @post_smearing_efficiencies = copy(%post_app_FHC_nue_bkg)
>

channel(#NuMode_Dplus_app_tau_numubar)<
        @channel =      #DplusNuMode_minus:     -:      t:      m:      #CC:    #app_numu_bkg
       @post_smearing_efficiencies = copy(%post_app_FHC_nue_bkg)
>

channel(#NuMode_Dplus_app_bkg_nutau)<
        @channel =      #DplusNuMode:     +:      m:      t:      #CCtau:    #app_nutau_bkg
       @post_smearing_efficiencies = copy(%post_app_FHC_nutau_bkg)
>

channel(#NuMode_Dplus_app_bkg_nutaubar)<
         @channel =      #DplusNuMode_minus:     -:      m:      t:      #CCtau:    #app_nutaubar_bkg
         @post_smearing_efficiencies = copy(%post_app_FHC_nutaubar_bkg)
>


channel(#NuMode_Dplus_dis_nutau)<
        @channel =      #DplusNuMode:     +:      t:      t:      #CCtau:    #app_nutau_bkg
       @post_smearing_efficiencies = copy(%post_app_FHC_nutau_bkg)
>

channel(#NuMode_Dplus_dis_nutaubar)<
         @channel =      #DplusNuMode_minus:     -:      t:      t:      #CCtau:    #app_nutaubar_bkg
         @post_smearing_efficiencies = copy(%post_app_FHC_nutaubar_bkg)
>


channel(#NuMode_Dplus_app_bkg_e_numu)<
        @channel =      #DplusNuMode:     +:      e:      m:      #CC:    #app_numu_bkg
       @post_smearing_efficiencies = copy(%post_app_FHC_numu_bkg)
>

channel(#NuMode_Dplus_app_bkg_e_numubar)<
         @channel =      #DplusNuMode_minus:     -:      e:      m:      #CC:    #app_numubar_bkg
         @post_smearing_efficiencies = copy(%post_app_FHC_numubar_bkg)
>


channel(#NuMode_Dplus_app_mu_nutau)<
        @channel =      #DplusNuMode:     +:      m:      t:      #CCtau:    #tau_hadronicwithbias
>

channel(#NuMode_Dplus_app_mu_nutaubar)<
         @channel =      #DplusNuMode_minus:     -:      m:      t:      #CCtau:    #tau_hadronicwithbias
>

channel(#NuMode_Dplus_app_tau_nutau)<
        @channel =      #DplusNuMode:     +:      t:      t:      #CCtau:    #tau_hadronicwithbias
>

channel(#NuMode_Dplus_app_tau_nutaubar)<
         @channel =      #DplusNuMode_minus:     -:      t:      t:      #CCtau:    #tau_hadronicwithbias
>

channel(#NuMode_Dplus_app_e_nutau)<
        @channel =      #DplusNuMode:     +:      e:      t:      #CCtau:    #tau_hadronicwithbias
>

channel(#NuMode_Dplus_app_e_nutaubar)<
         @channel =      #DplusNuMode_minus:     -:      e:      t:      #CCtau:    #tau_hadronicwithbias
>

channel(#NuMode_Dplus_app_bkg_nuNC_noeff)<
        @channel =      #DplusNuMode:     +:      m: NOSC_m:      #NC:    #app_NC_bkg
>

channel(#NuMode_Dplus_app_bkg_nubarNC_noeff)<
        @channel =      #DplusNuMode_minus:     -:      m: NOSC_m:      #NC:    #app_NC_bkg
>

channel(#NuMode_Dplus_app_bkg_nuNC_e_noeff)<
        @channel =      #DplusNuMode:     +:      e: NOSC_e:      #NC:    #app_NC_bkg
>

channel(#NuMode_Dplus_app_bkg_nubarNC_e_noeff)<
        @channel =      #DplusNuMode_minus:     -:      e: NOSC_e:      #NC:    #app_NC_bkg
>

channel(#NuMode_Dplus_app_bkg_nuNC_t_noeff)<
        @channel =      #DplusNuMode:     +:      t: NOSC_t:      #NC:    #app_NC_bkg
>

channel(#NuMode_Dplus_app_bkg_nubarNC_t_noeff)<
        @channel =      #DplusNuMode_minus:     -:      t: NOSC_t:      #NC:    #app_NC_bkg
>




//Ds


channel(#NuMode_DS_app_osc_nue)<
        @channel =      #DSNuMode:     +:      m:      e:      #CC:    #app_nue_sig
       @post_smearing_efficiencies = copy(%post_app_FHC_nue_sig)
>

channel(#NuMode_DS_app_osc_nuebar)<
         @channel =      #DSNuMode_minus:     -:      m:      e:      #CC:    #app_nuebar_sig
         @post_smearing_efficiencies = copy(%post_app_FHC_nuebar_sig)
>

channel(#NuMode_DS_app_bkg_nue)<
        @channel =      #DSNuMode:     +:      e:      e:      #CC:    #app_nue_bkg
       @post_smearing_efficiencies = copy(%post_app_FHC_nue_bkg)
>

channel(#NuMode_DS_app_bkg_nuebar)<
         @channel =      #DSNuMode_minus:     -:      e:      e:      #CC:    #app_nuebar_bkg
         @post_smearing_efficiencies = copy(%post_app_FHC_nuebar_bkg)
>

channel(#NuMode_DS_app_bkg_numu)<
        @channel =      #DSNuMode:     +:      m:      m:      #CC:    #app_numu_bkg
       @post_smearing_efficiencies = copy(%post_app_FHC_numu_bkg)
>

channel(#NuMode_DS_app_bkg_numubar)<
         @channel =      #DSNuMode_minus:     -:      m:      m:      #CC:    #app_numubar_bkg
         @post_smearing_efficiencies = copy(%post_app_FHC_numubar_bkg)
>


channel(#NuMode_DS_app_bkg_e_nutau)<
        @channel =      #DSNuMode:     +:      e:      t:      #CCtau:    #app_nutau_bkg
       @post_smearing_efficiencies = copy(%post_app_FHC_nutau_bkg)
>

channel(#NuMode_DS_app_bkg_e_nutaubar)<
         @channel =      #DSNuMode_minus:     -:      e:      t:      #CCtau:    #app_nutaubar_bkg
         @post_smearing_efficiencies = copy(%post_app_FHC_nutaubar_bkg)
>



channel(#NuMode_DS_app_bkg_nuNC)<
        @channel =      #DSNuMode:     +:      m: NOSC_m:      #NC:    #app_NC_bkg
        @post_smearing_efficiencies = copy(%post_app_FHC_NC_bkg)
>

channel(#NuMode_DS_app_bkg_nubarNC)<
        @channel =      #DSNuMode_minus:     -:      m: NOSC_m:      #NC:    #app_NC_bkg
        @post_smearing_efficiencies = copy(%post_app_FHC_aNC_bkg)
>

channel(#NuMode_DS_app_bkg_nuNC_e)<
        @channel =      #DSNuMode:     +:      e: NOSC_e:      #NC:    #app_NC_bkg
        @post_smearing_efficiencies = copy(%post_app_FHC_NC_bkg)
>

channel(#NuMode_DS_app_bkg_nubarNC_e)<
        @channel =      #DSNuMode_minus:     -:      e: NOSC_e:      #NC:    #app_NC_bkg
        @post_smearing_efficiencies = copy(%post_app_FHC_aNC_bkg)
>

channel(#NuMode_DS_app_bkg_nuNC_t)<
        @channel =      #DSNuMode:     +:      t: NOSC_t:      #NC:    #app_NC_bkg
        @post_smearing_efficiencies = copy(%post_app_FHC_NC_bkg)
>

channel(#NuMode_DS_app_bkg_nubarNC_t)<
        @channel =      #DSNuMode_minus:     -:      t: NOSC_t:      #NC:    #app_NC_bkg
        @post_smearing_efficiencies = copy(%post_app_FHC_aNC_bkg)
>

channel(#NuMode_DS_app_tau_e)<
        @channel =      #DSNuMode:     +:      t:      e:      #CC:    #app_nue_bkg
       @post_smearing_efficiencies = copy(%post_app_FHC_nue_bkg)
>

channel(#NuMode_DS_app_tau_ebar)<
        @channel =      #DSNuMode_minus:     -:      t:      e:      #CC:    #app_nue_bkg
       @post_smearing_efficiencies = copy(%post_app_FHC_nue_bkg)
>

channel(#NuMode_DS_app_tau_numu)<
        @channel =      #DSNuMode:     +:      t:      m:      #CC:    #app_numu_bkg
       @post_smearing_efficiencies = copy(%post_app_FHC_nue_bkg)
>

channel(#NuMode_DS_app_tau_numubar)<
        @channel =      #DSNuMode_minus:     -:      t:      m:      #CC:    #app_numu_bkg
       @post_smearing_efficiencies = copy(%post_app_FHC_nue_bkg)
>

channel(#NuMode_DS_app_bkg_nutau)<
        @channel =      #DSNuMode:     +:      m:      t:      #CCtau:    #app_nutau_bkg
       @post_smearing_efficiencies = copy(%post_app_FHC_nutau_bkg)
>

channel(#NuMode_DS_app_bkg_nutaubar)<
         @channel =      #DSNuMode_minus:     -:      m:      t:      #CCtau:    #app_nutaubar_bkg
         @post_smearing_efficiencies = copy(%post_app_FHC_nutaubar_bkg)
>


channel(#NuMode_DS_dis_nutau)<
        @channel =      #DSNuMode:     +:      t:      t:      #CCtau:    #app_nutau_bkg
       @post_smearing_efficiencies = copy(%post_app_FHC_nutau_bkg)
>

channel(#NuMode_DS_dis_nutaubar)<
         @channel =      #DSNuMode_minus:     -:      t:      t:      #CCtau:    #app_nutaubar_bkg
         @post_smearing_efficiencies = copy(%post_app_FHC_nutaubar_bkg)
>


channel(#NuMode_DS_app_bkg_e_numu)<
        @channel =      #DSNuMode:     +:      e:      m:      #CC:    #app_numu_bkg
       @post_smearing_efficiencies = copy(%post_app_FHC_numu_bkg)
>

channel(#NuMode_DS_app_bkg_e_numubar)<
         @channel =      #DSNuMode_minus:     -:      e:      m:      #CC:    #app_numubar_bkg
         @post_smearing_efficiencies = copy(%post_app_FHC_numubar_bkg)
>

channel(#NuMode_DS_app_mu_nutau)<
        @channel =      #DSNuMode:     +:      m:      t:      #CCtau:    #tau_hadronicwithbias
>

channel(#NuMode_DS_app_mu_nutaubar)<
         @channel =      #DSNuMode_minus:     -:      m:      t:      #CCtau:    #tau_hadronicwithbias
>

channel(#NuMode_DS_app_tau_nutau)<
        @channel =      #DSNuMode:     +:      t:      t:      #CCtau:    #tau_hadronicwithbias
>

channel(#NuMode_DS_app_tau_nutaubar)<
         @channel =      #DSNuMode_minus:     -:      t:      t:      #CCtau:    #tau_hadronicwithbias
>

channel(#NuMode_DS_app_e_nutau)<
        @channel =      #DSNuMode:     +:      e:      t:      #CCtau:    #tau_hadronicwithbias
>

channel(#NuMode_DS_app_e_nutaubar)<
         @channel =      #DSNuMode_minus:     -:      e:      t:      #CCtau:    #tau_hadronicwithbias
>

channel(#NuMode_DS_app_bkg_nuNC_noeff)<
        @channel =      #DSNuMode:     +:      m: NOSC_m:      #NC:    #app_NC_bkg
>

channel(#NuMode_DS_app_bkg_nubarNC_noeff)<
        @channel =      #DSNuMode_minus:     -:      m: NOSC_m:      #NC:    #app_NC_bkg
>

channel(#NuMode_DS_app_bkg_nuNC_e_noeff)<
        @channel =      #DSNuMode:     +:      e: NOSC_e:      #NC:    #app_NC_bkg
>

channel(#NuMode_DS_app_bkg_nubarNC_e_noeff)<
        @channel =      #DSNuMode_minus:     -:      e: NOSC_e:      #NC:    #app_NC_bkg
>

channel(#NuMode_DS_app_bkg_nuNC_t_noeff)<
        @channel =      #DSNuMode:     +:      t: NOSC_t:      #NC:    #app_NC_bkg
>

channel(#NuMode_DS_app_bkg_nubarNC_t_noeff)<
        @channel =      #DSNuMode_minus:     -:      t: NOSC_t:      #NC:    #app_NC_bkg
>




//K

channel(#NuMode_Kaon_app_osc_nue)<
        @channel =      #KaonNuMode:     +:      m:      e:      #CC:    #app_nue_sig
       @post_smearing_efficiencies = copy(%post_app_FHC_nue_sig)
>

channel(#NuMode_Kaon_app_osc_nuebar)<
         @channel =      #KaonNuMode_minus:     -:      m:      e:      #CC:    #app_nuebar_sig
         @post_smearing_efficiencies = copy(%post_app_FHC_nuebar_sig)
>

channel(#NuMode_Kaon_app_bkg_nue)<
        @channel =      #KaonNuMode:     +:      e:      e:      #CC:    #app_nue_bkg
       @post_smearing_efficiencies = copy(%post_app_FHC_nue_bkg)
>

channel(#NuMode_Kaon_app_bkg_nuebar)<
         @channel =      #KaonNuMode_minus:     -:      e:      e:      #CC:    #app_nuebar_bkg
         @post_smearing_efficiencies = copy(%post_app_FHC_nuebar_bkg)
>

channel(#NuMode_Kaon_app_bkg_numu)<
        @channel =      #KaonNuMode:     +:      m:      m:      #CC:    #app_numu_bkg
       @post_smearing_efficiencies = copy(%post_app_FHC_numu_bkg)
>

channel(#NuMode_Kaon_app_bkg_numubar)<
         @channel =      #KaonNuMode_minus:     -:      m:      m:      #CC:    #app_numubar_bkg
         @post_smearing_efficiencies = copy(%post_app_FHC_numubar_bkg)
>


channel(#NuMode_Kaon_app_bkg_e_nutau)<
        @channel =      #KaonNuMode:     +:      e:      t:      #CCtau:    #app_nutau_bkg
       @post_smearing_efficiencies = copy(%post_app_FHC_nutau_bkg)
>

channel(#NuMode_Kaon_app_bkg_e_nutaubar)<
         @channel =      #KaonNuMode_minus:     -:      e:      t:      #CCtau:    #app_nutaubar_bkg
         @post_smearing_efficiencies = copy(%post_app_FHC_nutaubar_bkg)
>



channel(#NuMode_Kaon_app_bkg_nuNC)<
        @channel =      #KaonNuMode:     +:      m: NOSC_m:      #NC:    #app_NC_bkg
        @post_smearing_efficiencies = copy(%post_app_FHC_NC_bkg)
>

channel(#NuMode_Kaon_app_bkg_nubarNC)<
        @channel =      #KaonNuMode_minus:     -:      m: NOSC_m:      #NC:    #app_NC_bkg
        @post_smearing_efficiencies = copy(%post_app_FHC_aNC_bkg)
>

channel(#NuMode_Kaon_app_bkg_nuNC_e)<
        @channel =      #KaonNuMode:     +:      e: NOSC_e:      #NC:    #app_NC_bkg
        @post_smearing_efficiencies = copy(%post_app_FHC_NC_bkg)
>

channel(#NuMode_Kaon_app_bkg_nubarNC_e)<
        @channel =      #KaonNuMode_minus:     -:      e: NOSC_e:      #NC:    #app_NC_bkg
        @post_smearing_efficiencies = copy(%post_app_FHC_aNC_bkg)
>

channel(#NuMode_Kaon_app_bkg_nuNC_t)<
        @channel =      #KaonNuMode:     +:      t: NOSC_t:      #NC:    #app_NC_bkg
        @post_smearing_efficiencies = copy(%post_app_FHC_NC_bkg)
>

channel(#NuMode_Kaon_app_bkg_nubarNC_t)<
        @channel =      #KaonNuMode_minus:     -:      t: NOSC_t:      #NC:    #app_NC_bkg
        @post_smearing_efficiencies = copy(%post_app_FHC_aNC_bkg)
>

channel(#NuMode_Kaon_app_tau_e)<
        @channel =      #KaonNuMode:     +:      t:      e:      #CC:    #app_nue_bkg
       @post_smearing_efficiencies = copy(%post_app_FHC_nue_bkg)
>

channel(#NuMode_Kaon_app_tau_ebar)<
        @channel =      #KaonNuMode_minus:     -:      t:      e:      #CC:    #app_nue_bkg
       @post_smearing_efficiencies = copy(%post_app_FHC_nue_bkg)
>

channel(#NuMode_Kaon_app_tau_numu)<
        @channel =      #KaonNuMode:     +:      t:      m:      #CC:    #app_numu_bkg
       @post_smearing_efficiencies = copy(%post_app_FHC_nue_bkg)
>

channel(#NuMode_Kaon_app_tau_numubar)<
        @channel =      #KaonNuMode_minus:     -:      t:      m:      #CC:    #app_numu_bkg
       @post_smearing_efficiencies = copy(%post_app_FHC_nue_bkg)
>

channel(#NuMode_Kaon_app_bkg_nutau)<
        @channel =      #KaonNuMode:     +:      m:      t:      #CCtau:    #app_nutau_bkg
       @post_smearing_efficiencies = copy(%post_app_FHC_nutau_bkg)
>

channel(#NuMode_Kaon_app_bkg_nutaubar)<
         @channel =      #KaonNuMode_minus:     -:      m:      t:      #CCtau:    #app_nutaubar_bkg
         @post_smearing_efficiencies = copy(%post_app_FHC_nutaubar_bkg)
>


channel(#NuMode_Kaon_dis_nutau)<
        @channel =      #KaonNuMode:     +:      t:      t:      #CCtau:    #app_nutau_bkg
       @post_smearing_efficiencies = copy(%post_app_FHC_nutau_bkg)
>

channel(#NuMode_Kaon_dis_nutaubar)<
         @channel =      #KaonNuMode_minus:     -:      t:      t:      #CCtau:    #app_nutaubar_bkg
         @post_smearing_efficiencies = copy(%post_app_FHC_nutaubar_bkg)
>


channel(#NuMode_Kaon_app_bkg_e_numu)<
        @channel =      #KaonNuMode:     +:      e:      m:      #CC:    #app_numu_bkg
       @post_smearing_efficiencies = copy(%post_app_FHC_numu_bkg)
>

channel(#NuMode_Kaon_app_bkg_e_numubar)<
         @channel =      #KaonNuMode_minus:     -:      e:      m:      #CC:    #app_numubar_bkg
         @post_smearing_efficiencies = copy(%post_app_FHC_numubar_bkg)
>

channel(#NuMode_Kaon_app_mu_nutau)<
        @channel =      #KaonNuMode:     +:      m:      t:      #CCtau:    #tau_hadronicwithbias
>

channel(#NuMode_Kaon_app_mu_nutaubar)<
         @channel =      #KaonNuMode_minus:     -:      m:      t:      #CCtau:    #tau_hadronicwithbias
>

channel(#NuMode_Kaon_app_tau_nutau)<
        @channel =      #KaonNuMode:     +:      t:      t:      #CCtau:    #tau_hadronicwithbias
>

channel(#NuMode_Kaon_app_tau_nutaubar)<
         @channel =      #KaonNuMode_minus:     -:      t:      t:      #CCtau:    #tau_hadronicwithbias
>

channel(#NuMode_Kaon_app_e_nutau)<
        @channel =      #KaonNuMode:     +:      e:      t:      #CCtau:    #tau_hadronicwithbias
>

channel(#NuMode_Kaon_app_e_nutaubar)<
         @channel =      #KaonNuMode_minus:     -:      e:      t:      #CCtau:    #tau_hadronicwithbias
>

channel(#NuMode_Kaon_app_bkg_nuNC_noeff)<
        @channel =      #KaonNuMode:     +:      m: NOSC_m:      #NC:    #app_NC_bkg
>

channel(#NuMode_Kaon_app_bkg_nubarNC_noeff)<
        @channel =      #KaonNuMode_minus:     -:      m: NOSC_m:      #NC:    #app_NC_bkg
>

channel(#NuMode_Kaon_app_bkg_nuNC_e_noeff)<
        @channel =      #KaonNuMode:     +:      e: NOSC_e:      #NC:    #app_NC_bkg
>

channel(#NuMode_Kaon_app_bkg_nubarNC_e_noeff)<
        @channel =      #KaonNuMode_minus:     -:      e: NOSC_e:      #NC:    #app_NC_bkg
>

channel(#NuMode_Kaon_app_bkg_nuNC_t_noeff)<
        @channel =      #KaonNuMode:     +:      t: NOSC_t:      #NC:    #app_NC_bkg
>

channel(#NuMode_Kaon_app_bkg_nubarNC_t_noeff)<
        @channel =      #KaonNuMode_minus:     -:      t: NOSC_t:      #NC:    #app_NC_bkg
>





//Pion

channel(#NuMode_Pion_app_osc_nue)<
        @channel =      #PionNuMode:     +:      m:      e:      #CC:    #app_nue_sig
       @post_smearing_efficiencies = copy(%post_app_FHC_nue_sig)
>

channel(#NuMode_Pion_app_osc_nuebar)<
         @channel =      #PionNuMode_minus:     -:      m:      e:      #CC:    #app_nuebar_sig
         @post_smearing_efficiencies = copy(%post_app_FHC_nuebar_sig)
>

channel(#NuMode_Pion_app_bkg_nue)<
        @channel =      #PionNuMode:     +:      e:      e:      #CC:    #app_nue_bkg
       @post_smearing_efficiencies = copy(%post_app_FHC_nue_bkg)
>

channel(#NuMode_Pion_app_bkg_nuebar)<
         @channel =      #PionNuMode_minus:     -:      e:      e:      #CC:    #app_nuebar_bkg
         @post_smearing_efficiencies = copy(%post_app_FHC_nuebar_bkg)
>

channel(#NuMode_Pion_app_bkg_numu)<
        @channel =      #PionNuMode:     +:      m:      m:      #CC:    #app_numu_bkg
       @post_smearing_efficiencies = copy(%post_app_FHC_numu_bkg)
>

channel(#NuMode_Pion_app_bkg_numubar)<
         @channel =      #PionNuMode_minus:     -:      m:      m:      #CC:    #app_numubar_bkg
         @post_smearing_efficiencies = copy(%post_app_FHC_numubar_bkg)
>


channel(#NuMode_Pion_app_bkg_e_nutau)<
        @channel =      #PionNuMode:     +:      e:      t:      #CCtau:    #app_nutau_bkg
       @post_smearing_efficiencies = copy(%post_app_FHC_nutau_bkg)
>

channel(#NuMode_Pion_app_bkg_e_nutaubar)<
         @channel =      #PionNuMode_minus:     -:      e:      t:      #CCtau:    #app_nutaubar_bkg
         @post_smearing_efficiencies = copy(%post_app_FHC_nutaubar_bkg)
>



channel(#NuMode_Pion_app_bkg_nuNC)<
        @channel =      #PionNuMode:     +:      m: NOSC_m:      #NC:    #app_NC_bkg
        @post_smearing_efficiencies = copy(%post_app_FHC_NC_bkg)
>

channel(#NuMode_Pion_app_bkg_nubarNC)<
        @channel =      #PionNuMode_minus:     -:      m: NOSC_m:      #NC:    #app_NC_bkg
        @post_smearing_efficiencies = copy(%post_app_FHC_aNC_bkg)
>

channel(#NuMode_Pion_app_bkg_nuNC_e)<
        @channel =      #PionNuMode:     +:      e: NOSC_e:      #NC:    #app_NC_bkg
        @post_smearing_efficiencies = copy(%post_app_FHC_NC_bkg)
>

channel(#NuMode_Pion_app_bkg_nubarNC_e)<
        @channel =      #PionNuMode_minus:     -:      e: NOSC_e:      #NC:    #app_NC_bkg
        @post_smearing_efficiencies = copy(%post_app_FHC_aNC_bkg)
>

channel(#NuMode_Pion_app_bkg_nuNC_t)<
        @channel =      #PionNuMode:     +:      t: NOSC_t:      #NC:    #app_NC_bkg
        @post_smearing_efficiencies = copy(%post_app_FHC_NC_bkg)
>

channel(#NuMode_Pion_app_bkg_nubarNC_t)<
        @channel =      #PionNuMode_minus:     -:      t: NOSC_t:      #NC:    #app_NC_bkg
        @post_smearing_efficiencies = copy(%post_app_FHC_aNC_bkg)
>

channel(#NuMode_Pion_app_tau_e)<
        @channel =      #PionNuMode:     +:      t:      e:      #CC:    #app_nue_bkg
       @post_smearing_efficiencies = copy(%post_app_FHC_nue_bkg)
>

channel(#NuMode_Pion_app_tau_ebar)<
        @channel =      #PionNuMode_minus:     -:      t:      e:      #CC:    #app_nue_bkg
       @post_smearing_efficiencies = copy(%post_app_FHC_nue_bkg)
>

channel(#NuMode_Pion_app_tau_numu)<
        @channel =      #PionNuMode:     +:      t:      m:      #CC:    #app_numu_bkg
       @post_smearing_efficiencies = copy(%post_app_FHC_nue_bkg)
>

channel(#NuMode_Pion_app_tau_numubar)<
        @channel =      #PionNuMode_minus:     -:      t:      m:      #CC:    #app_numu_bkg
       @post_smearing_efficiencies = copy(%post_app_FHC_nue_bkg)
>

channel(#NuMode_Pion_app_bkg_nutau)<
        @channel =      #PionNuMode:     +:      m:      t:      #CCtau:    #app_nutau_bkg
       @post_smearing_efficiencies = copy(%post_app_FHC_nutau_bkg)
>

channel(#NuMode_Pion_app_bkg_nutaubar)<
         @channel =      #PionNuMode_minus:     -:      m:      t:      #CCtau:    #app_nutaubar_bkg
         @post_smearing_efficiencies = copy(%post_app_FHC_nutaubar_bkg)
>


channel(#NuMode_Pion_dis_nutau)<
        @channel =      #PionNuMode:     +:      t:      t:      #CCtau:    #app_nutau_bkg
       @post_smearing_efficiencies = copy(%post_app_FHC_nutau_bkg)
>

channel(#NuMode_Pion_dis_nutaubar)<
         @channel =      #PionNuMode_minus:     -:      t:      t:      #CCtau:    #app_nutaubar_bkg
         @post_smearing_efficiencies = copy(%post_app_FHC_nutaubar_bkg)
>


channel(#NuMode_Pion_app_bkg_e_numu)<
        @channel =      #PionNuMode:     +:      e:      m:      #CC:    #app_numu_bkg
       @post_smearing_efficiencies = copy(%post_app_FHC_numu_bkg)
>

channel(#NuMode_Pion_app_bkg_e_numubar)<
         @channel =      #PionNuMode_minus:     -:      e:      m:      #CC:    #app_numubar_bkg
         @post_smearing_efficiencies = copy(%post_app_FHC_numubar_bkg)
>

channel(#NuMode_Pion_app_mu_nutau)<
        @channel =      #PionNuMode:     +:      m:      t:      #CCtau:    #tau_hadronicwithbias
>

channel(#NuMode_Pion_app_mu_nutaubar)<
         @channel =      #PionNuMode_minus:     -:      m:      t:      #CCtau:    #tau_hadronicwithbias
>

channel(#NuMode_Pion_app_tau_nutau)<
        @channel =      #PionNuMode:     +:      t:      t:      #CCtau:    #tau_hadronicwithbias
>

channel(#NuMode_Pion_app_tau_nutaubar)<
         @channel =      #PionNuMode_minus:     -:      t:      t:      #CCtau:    #tau_hadronicwithbias
>

channel(#NuMode_Pion_app_e_nutau)<
        @channel =      #PionNuMode:     +:      e:      t:      #CCtau:    #tau_hadronicwithbias
>

channel(#NuMode_Pion_app_e_nutaubar)<
         @channel =      #PionNuMode_minus:     -:      e:      t:      #CCtau:    #tau_hadronicwithbias
>

channel(#NuMode_Pion_app_bkg_nuNC_noeff)<
        @channel =      #PionNuMode:     +:      m: NOSC_m:      #NC:    #app_NC_bkg
>

channel(#NuMode_Pion_app_bkg_nubarNC_noeff)<
        @channel =      #PionNuMode_minus:     -:      m: NOSC_m:      #NC:    #app_NC_bkg
>

channel(#NuMode_Pion_app_bkg_nuNC_e_noeff)<
        @channel =      #PionNuMode:     +:      e: NOSC_e:      #NC:    #app_NC_bkg
>

channel(#NuMode_Pion_app_bkg_nubarNC_e_noeff)<
        @channel =      #PionNuMode_minus:     -:      e: NOSC_e:      #NC:    #app_NC_bkg
>

channel(#NuMode_Pion_app_bkg_nuNC_t_noeff)<
        @channel =      #PionNuMode:     +:      t: NOSC_t:      #NC:    #app_NC_bkg
>

channel(#NuMode_Pion_app_bkg_nubarNC_t_noeff)<
        @channel =      #PionNuMode_minus:     -:      t: NOSC_t:      #NC:    #app_NC_bkg
>




//Muon


channel(#NuMode_Muon_app_osc_nue)<
        @channel =      #MuonNuMode:     +:      m:      e:      #CC:    #app_nue_sig
       @post_smearing_efficiencies = copy(%post_app_FHC_nue_sig)
>

channel(#NuMode_Muon_app_osc_nuebar)<
         @channel =      #MuonNuMode_minus:     -:      m:      e:      #CC:    #app_nuebar_sig
         @post_smearing_efficiencies = copy(%post_app_FHC_nuebar_sig)
>

channel(#NuMode_Muon_app_bkg_nue)<
        @channel =      #MuonNuMode:     +:      e:      e:      #CC:    #app_nue_bkg
       @post_smearing_efficiencies = copy(%post_app_FHC_nue_bkg)
>

channel(#NuMode_Muon_app_bkg_nuebar)<
         @channel =      #MuonNuMode_minus:     -:      e:      e:      #CC:    #app_nuebar_bkg
         @post_smearing_efficiencies = copy(%post_app_FHC_nuebar_bkg)
>

channel(#NuMode_Muon_app_bkg_numu)<
        @channel =      #MuonNuMode:     +:      m:      m:      #CC:    #app_numu_bkg
       @post_smearing_efficiencies = copy(%post_app_FHC_numu_bkg)
>

channel(#NuMode_Muon_app_bkg_numubar)<
         @channel =      #MuonNuMode_minus:     -:      m:      m:      #CC:    #app_numubar_bkg
         @post_smearing_efficiencies = copy(%post_app_FHC_numubar_bkg)
>


channel(#NuMode_Muon_app_bkg_e_nutau)<
        @channel =      #MuonNuMode:     +:      e:      t:      #CCtau:    #app_nutau_bkg
       @post_smearing_efficiencies = copy(%post_app_FHC_nutau_bkg)
>

channel(#NuMode_Muon_app_bkg_e_nutaubar)<
         @channel =      #MuonNuMode_minus:     -:      e:      t:      #CCtau:    #app_nutaubar_bkg
         @post_smearing_efficiencies = copy(%post_app_FHC_nutaubar_bkg)
>



channel(#NuMode_Muon_app_bkg_nuNC)<
        @channel =      #MuonNuMode:     +:      m: NOSC_m:      #NC:    #app_NC_bkg
        @post_smearing_efficiencies = copy(%post_app_FHC_NC_bkg)
>

channel(#NuMode_Muon_app_bkg_nubarNC)<
        @channel =      #MuonNuMode_minus:     -:      m: NOSC_m:      #NC:    #app_NC_bkg
        @post_smearing_efficiencies = copy(%post_app_FHC_aNC_bkg)
>

channel(#NuMode_Muon_app_bkg_nuNC_e)<
        @channel =      #MuonNuMode:     +:      e: NOSC_e:      #NC:    #app_NC_bkg
        @post_smearing_efficiencies = copy(%post_app_FHC_NC_bkg)
>

channel(#NuMode_Muon_app_bkg_nubarNC_e)<
        @channel =      #MuonNuMode_minus:     -:      e: NOSC_e:      #NC:    #app_NC_bkg
        @post_smearing_efficiencies = copy(%post_app_FHC_aNC_bkg)
>

channel(#NuMode_Muon_app_bkg_nuNC_t)<
        @channel =      #MuonNuMode:     +:      t: NOSC_t:      #NC:    #app_NC_bkg
        @post_smearing_efficiencies = copy(%post_app_FHC_NC_bkg)
>

channel(#NuMode_Muon_app_bkg_nubarNC_t)<
        @channel =      #MuonNuMode_minus:     -:      t: NOSC_t:      #NC:    #app_NC_bkg
        @post_smearing_efficiencies = copy(%post_app_FHC_aNC_bkg)
>

channel(#NuMode_Muon_app_tau_e)<
        @channel =      #MuonNuMode:     +:      t:      e:      #CC:    #app_nue_bkg
       @post_smearing_efficiencies = copy(%post_app_FHC_nue_bkg)
>

channel(#NuMode_Muon_app_tau_ebar)<
        @channel =      #MuonNuMode_minus:     -:      t:      e:      #CC:    #app_nue_bkg
       @post_smearing_efficiencies = copy(%post_app_FHC_nue_bkg)
>

channel(#NuMode_Muon_app_tau_numu)<
        @channel =      #MuonNuMode:     +:      t:      m:      #CC:    #app_numu_bkg
       @post_smearing_efficiencies = copy(%post_app_FHC_nue_bkg)
>

channel(#NuMode_Muon_app_tau_numubar)<
        @channel =      #MuonNuMode_minus:     -:      t:      m:      #CC:    #app_numu_bkg
       @post_smearing_efficiencies = copy(%post_app_FHC_nue_bkg)
>

channel(#NuMode_Muon_app_bkg_nutau)<
        @channel =      #MuonNuMode:     +:      m:      t:      #CCtau:    #app_nutau_bkg
       @post_smearing_efficiencies = copy(%post_app_FHC_nutau_bkg)
>

channel(#NuMode_Muon_app_bkg_nutaubar)<
         @channel =      #MuonNuMode_minus:     -:      m:      t:      #CCtau:    #app_nutaubar_bkg
         @post_smearing_efficiencies = copy(%post_app_FHC_nutaubar_bkg)
>


channel(#NuMode_Muon_dis_nutau)<
        @channel =      #MuonNuMode:     +:      t:      t:      #CCtau:    #app_nutau_bkg
       @post_smearing_efficiencies = copy(%post_app_FHC_nutau_bkg)
>

channel(#NuMode_Muon_dis_nutaubar)<
         @channel =      #MuonNuMode_minus:     -:      t:      t:      #CCtau:    #app_nutaubar_bkg
         @post_smearing_efficiencies = copy(%post_app_FHC_nutaubar_bkg)
>


channel(#NuMode_Muon_app_bkg_e_numu)<
        @channel =      #MuonNuMode:     +:      e:      m:      #CC:    #app_numu_bkg
       @post_smearing_efficiencies = copy(%post_app_FHC_numu_bkg)
>

channel(#NuMode_Muon_app_bkg_e_numubar)<
         @channel =      #MuonNuMode_minus:     -:      e:      m:      #CC:    #app_numubar_bkg
         @post_smearing_efficiencies = copy(%post_app_FHC_numubar_bkg)
>

channel(#NuMode_Muon_app_mu_nutau)<
        @channel =      #MuonNuMode:     +:      m:      t:      #CCtau:    #tau_hadronicwithbias
>

channel(#NuMode_Muon_app_mu_nutaubar)<
         @channel =      #MuonNuMode_minus:     -:      m:      t:      #CCtau:    #tau_hadronicwithbias
>

channel(#NuMode_Muon_app_tau_nutau)<
        @channel =      #MuonNuMode:     +:      t:      t:      #CCtau:    #tau_hadronicwithbias
>

channel(#NuMode_Muon_app_tau_nutaubar)<
         @channel =      #MuonNuMode_minus:     -:      t:      t:      #CCtau:    #tau_hadronicwithbias
>

channel(#NuMode_Muon_app_e_nutau)<
        @channel =      #MuonNuMode:     +:      e:      t:      #CCtau:    #tau_hadronicwithbias
>

channel(#NuMode_Muon_app_e_nutaubar)<
         @channel =      #MuonNuMode_minus:     -:      e:      t:      #CCtau:    #tau_hadronicwithbias
>

channel(#NuMode_Muon_app_bkg_nuNC_noeff)<
        @channel =      #MuonNuMode:     +:      m: NOSC_m:      #NC:    #app_NC_bkg
>

channel(#NuMode_Muon_app_bkg_nubarNC_noeff)<
        @channel =      #MuonNuMode_minus:     -:      m: NOSC_m:      #NC:    #app_NC_bkg
>

channel(#NuMode_Muon_app_bkg_nuNC_e_noeff)<
        @channel =      #MuonNuMode:     +:      e: NOSC_e:      #NC:    #app_NC_bkg
>

channel(#NuMode_Muon_app_bkg_nubarNC_e_noeff)<
        @channel =      #MuonNuMode_minus:     -:      e: NOSC_e:      #NC:    #app_NC_bkg
>

channel(#NuMode_Muon_app_bkg_nuNC_t_noeff)<
        @channel =      #MuonNuMode:     +:      t: NOSC_t:      #NC:    #app_NC_bkg
>

channel(#NuMode_Muon_app_bkg_nubarNC_t_noeff)<
        @channel =      #MuonNuMode_minus:     -:      t: NOSC_t:      #NC:    #app_NC_bkg
>





/* NUEBAR APP */


//D

channel(#NuBarMode_Dplus_app_osc_nue)<
        @channel =      #DplusNuBarMode:     +:      m:      e:      #CC:    #app_nue_sig
       @post_smearing_efficiencies = copy(%post_app_RHC_nue_sig)
>

channel(#NuBarMode_Dplus_app_osc_nuebar)<
         @channel =      #DplusNuBarMode_minus:     -:      m:      e:      #CC:    #app_nuebar_sig
         @post_smearing_efficiencies = copy(%post_app_RHC_nuebar_sig)
>

channel(#NuBarMode_Dplus_app_bkg_nue)<
        @channel =      #DplusNuBarMode:     +:      e:      e:      #CC:    #app_nue_bkg
       @post_smearing_efficiencies = copy(%post_app_RHC_nue_bkg)
>

channel(#NuBarMode_Dplus_app_bkg_nuebar)<
         @channel =      #DplusNuBarMode_minus:     -:      e:      e:      #CC:    #app_nuebar_bkg
         @post_smearing_efficiencies = copy(%post_app_RHC_nuebar_bkg)
>

channel(#NuBarMode_Dplus_app_bkg_numu)<
        @channel =      #DplusNuBarMode:     +:      m:      m:      #CC:    #app_numu_bkg
       @post_smearing_efficiencies = copy(%post_app_RHC_numu_bkg)
>

channel(#NuBarMode_Dplus_app_bkg_numubar)<
         @channel =      #DplusNuBarMode_minus:     -:      m:      m:      #CC:    #app_numubar_bkg
         @post_smearing_efficiencies = copy(%post_app_RHC_numubar_bkg)
>


channel(#NuBarMode_Dplus_app_bkg_e_nutau)<
        @channel =      #DplusNuBarMode:     +:      e:      t:      #CCtau:    #app_nutau_bkg
       @post_smearing_efficiencies = copy(%post_app_RHC_nutau_bkg)
>

channel(#NuBarMode_Dplus_app_bkg_e_nutaubar)<
         @channel =      #DplusNuBarMode_minus:     -:      e:      t:      #CCtau:    #app_nutaubar_bkg
         @post_smearing_efficiencies = copy(%post_app_RHC_nutaubar_bkg)
>



channel(#NuBarMode_Dplus_app_bkg_nuNC)<
        @channel =      #DplusNuBarMode:     +:      m: NOSC_m:      #NC:    #app_NC_bkg
        @post_smearing_efficiencies = copy(%post_app_RHC_NC_bkg)
>

channel(#NuBarMode_Dplus_app_bkg_nubarNC)<
        @channel =      #DplusNuBarMode_minus:     -:      m: NOSC_m:      #NC:    #app_NC_bkg
        @post_smearing_efficiencies = copy(%post_app_RHC_aNC_bkg)
>

channel(#NuBarMode_Dplus_app_bkg_nuNC_e)<
        @channel =      #DplusNuBarMode:     +:      e: NOSC_e:      #NC:    #app_NC_bkg
        @post_smearing_efficiencies = copy(%post_app_RHC_NC_bkg)
>

channel(#NuBarMode_Dplus_app_bkg_nubarNC_e)<
        @channel =      #DplusNuBarMode_minus:     -:      e: NOSC_e:      #NC:    #app_NC_bkg
        @post_smearing_efficiencies = copy(%post_app_RHC_aNC_bkg)
>

channel(#NuBarMode_Dplus_app_bkg_nuNC_t)<
        @channel =      #DplusNuBarMode:     +:      t: NOSC_t:      #NC:    #app_NC_bkg
        @post_smearing_efficiencies = copy(%post_app_RHC_NC_bkg)
>

channel(#NuBarMode_Dplus_app_bkg_nubarNC_t)<
        @channel =      #DplusNuBarMode_minus:     -:      t: NOSC_t:      #NC:    #app_NC_bkg
        @post_smearing_efficiencies = copy(%post_app_RHC_aNC_bkg)
>

channel(#NuBarMode_Dplus_app_tau_e)<
        @channel =      #DplusNuBarMode:     +:      t:      e:      #CC:    #app_nue_bkg
       @post_smearing_efficiencies = copy(%post_app_RHC_nue_bkg)
>

channel(#NuBarMode_Dplus_app_tau_ebar)<
        @channel =      #DplusNuBarMode_minus:     -:      t:      e:      #CC:    #app_nue_bkg
       @post_smearing_efficiencies = copy(%post_app_RHC_nue_bkg)
>

channel(#NuBarMode_Dplus_app_tau_numu)<
        @channel =      #DplusNuBarMode:     +:      t:      m:      #CC:    #app_numu_bkg
       @post_smearing_efficiencies = copy(%post_app_RHC_nue_bkg)
>

channel(#NuBarMode_Dplus_app_tau_numubar)<
        @channel =      #DplusNuBarMode_minus:     -:      t:      m:      #CC:    #app_numu_bkg
       @post_smearing_efficiencies = copy(%post_app_RHC_nue_bkg)
>

channel(#NuBarMode_Dplus_app_bkg_nutau)<
        @channel =      #DplusNuBarMode:     +:      m:      t:      #CCtau:    #app_nutau_bkg
       @post_smearing_efficiencies = copy(%post_app_RHC_nutau_bkg)
>

channel(#NuBarMode_Dplus_app_bkg_nutaubar)<
         @channel =      #DplusNuBarMode_minus:     -:      m:      t:      #CCtau:    #app_nutaubar_bkg
         @post_smearing_efficiencies = copy(%post_app_RHC_nutaubar_bkg)
>


channel(#NuBarMode_Dplus_dis_nutau)<
        @channel =      #DplusNuBarMode:     +:      t:      t:      #CCtau:    #app_nutau_bkg
       @post_smearing_efficiencies = copy(%post_app_RHC_nutau_bkg)
>

channel(#NuBarMode_Dplus_dis_nutaubar)<
         @channel =      #DplusNuBarMode_minus:     -:      t:      t:      #CCtau:    #app_nutaubar_bkg
         @post_smearing_efficiencies = copy(%post_app_RHC_nutaubar_bkg)
>


channel(#NuBarMode_Dplus_app_bkg_e_numu)<
        @channel =      #DplusNuBarMode:     +:      e:      m:      #CC:    #app_numu_bkg
       @post_smearing_efficiencies = copy(%post_app_RHC_numu_bkg)
>

channel(#NuBarMode_Dplus_app_bkg_e_numubar)<
         @channel =      #DplusNuBarMode_minus:     -:      e:      m:      #CC:    #app_numubar_bkg
         @post_smearing_efficiencies = copy(%post_app_RHC_numubar_bkg)
>


channel(#NuBarMode_Dplus_app_mu_nutau)<
        @channel =      #DplusNuBarMode:     +:      m:      t:      #CCtau:    #tau_hadronicwithbias
>

channel(#NuBarMode_Dplus_app_mu_nutaubar)<
         @channel =      #DplusNuBarMode_minus:     -:      m:      t:      #CCtau:    #tau_hadronicwithbias
>

channel(#NuBarMode_Dplus_app_tau_nutau)<
        @channel =      #DplusNuBarMode:     +:      t:      t:      #CCtau:    #tau_hadronicwithbias
>

channel(#NuBarMode_Dplus_app_tau_nutaubar)<
         @channel =      #DplusNuBarMode_minus:     -:      t:      t:      #CCtau:    #tau_hadronicwithbias
>

channel(#NuBarMode_Dplus_app_e_nutau)<
        @channel =      #DplusNuBarMode:     +:      e:      t:      #CCtau:    #tau_hadronicwithbias
>

channel(#NuBarMode_Dplus_app_e_nutaubar)<
         @channel =      #DplusNuBarMode_minus:     -:      e:      t:      #CCtau:    #tau_hadronicwithbias
>

channel(#NuBarMode_Dplus_app_bkg_nuNC_noeff)<
        @channel =      #DplusNuBarMode:     +:      m: NOSC_m:      #NC:    #app_NC_bkg
>

channel(#NuBarMode_Dplus_app_bkg_nubarNC_noeff)<
        @channel =      #DplusNuBarMode_minus:     -:      m: NOSC_m:      #NC:    #app_NC_bkg
>

channel(#NuBarMode_Dplus_app_bkg_nuNC_e_noeff)<
        @channel =      #DplusNuBarMode:     +:      e: NOSC_e:      #NC:    #app_NC_bkg
>

channel(#NuBarMode_Dplus_app_bkg_nubarNC_e_noeff)<
        @channel =      #DplusNuBarMode_minus:     -:      e: NOSC_e:      #NC:    #app_NC_bkg
>

channel(#NuBarMode_Dplus_app_bkg_nuNC_t_noeff)<
        @channel =      #DplusNuBarMode:     +:      t: NOSC_t:      #NC:    #app_NC_bkg
>

channel(#NuBarMode_Dplus_app_bkg_nubarNC_t_noeff)<
        @channel =      #DplusNuBarMode_minus:     -:      t: NOSC_t:      #NC:    #app_NC_bkg
>




//Ds


channel(#NuBarMode_DS_app_osc_nue)<
        @channel =      #DSNuBarMode:     +:      m:      e:      #CC:    #app_nue_sig
       @post_smearing_efficiencies = copy(%post_app_RHC_nue_sig)
>

channel(#NuBarMode_DS_app_osc_nuebar)<
         @channel =      #DSNuBarMode_minus:     -:      m:      e:      #CC:    #app_nuebar_sig
         @post_smearing_efficiencies = copy(%post_app_RHC_nuebar_sig)
>

channel(#NuBarMode_DS_app_bkg_nue)<
        @channel =      #DSNuBarMode:     +:      e:      e:      #CC:    #app_nue_bkg
       @post_smearing_efficiencies = copy(%post_app_RHC_nue_bkg)
>

channel(#NuBarMode_DS_app_bkg_nuebar)<
         @channel =      #DSNuBarMode_minus:     -:      e:      e:      #CC:    #app_nuebar_bkg
         @post_smearing_efficiencies = copy(%post_app_RHC_nuebar_bkg)
>

channel(#NuBarMode_DS_app_bkg_numu)<
        @channel =      #DSNuBarMode:     +:      m:      m:      #CC:    #app_numu_bkg
       @post_smearing_efficiencies = copy(%post_app_RHC_numu_bkg)
>

channel(#NuBarMode_DS_app_bkg_numubar)<
         @channel =      #DSNuBarMode_minus:     -:      m:      m:      #CC:    #app_numubar_bkg
         @post_smearing_efficiencies = copy(%post_app_RHC_numubar_bkg)
>


channel(#NuBarMode_DS_app_bkg_e_nutau)<
        @channel =      #DSNuBarMode:     +:      e:      t:      #CCtau:    #app_nutau_bkg
       @post_smearing_efficiencies = copy(%post_app_RHC_nutau_bkg)
>

channel(#NuBarMode_DS_app_bkg_e_nutaubar)<
         @channel =      #DSNuBarMode_minus:     -:      e:      t:      #CCtau:    #app_nutaubar_bkg
         @post_smearing_efficiencies = copy(%post_app_RHC_nutaubar_bkg)
>



channel(#NuBarMode_DS_app_bkg_nuNC)<
        @channel =      #DSNuBarMode:     +:      m: NOSC_m:      #NC:    #app_NC_bkg
        @post_smearing_efficiencies = copy(%post_app_RHC_NC_bkg)
>

channel(#NuBarMode_DS_app_bkg_nubarNC)<
        @channel =      #DSNuBarMode_minus:     -:      m: NOSC_m:      #NC:    #app_NC_bkg
        @post_smearing_efficiencies = copy(%post_app_RHC_aNC_bkg)
>

channel(#NuBarMode_DS_app_bkg_nuNC_e)<
        @channel =      #DSNuBarMode:     +:      e: NOSC_e:      #NC:    #app_NC_bkg
        @post_smearing_efficiencies = copy(%post_app_RHC_NC_bkg)
>

channel(#NuBarMode_DS_app_bkg_nubarNC_e)<
        @channel =      #DSNuBarMode_minus:     -:      e: NOSC_e:      #NC:    #app_NC_bkg
        @post_smearing_efficiencies = copy(%post_app_RHC_aNC_bkg)
>

channel(#NuBarMode_DS_app_bkg_nuNC_t)<
        @channel =      #DSNuBarMode:     +:      t: NOSC_t:      #NC:    #app_NC_bkg
        @post_smearing_efficiencies = copy(%post_app_RHC_NC_bkg)
>

channel(#NuBarMode_DS_app_bkg_nubarNC_t)<
        @channel =      #DSNuBarMode_minus:     -:      t: NOSC_t:      #NC:    #app_NC_bkg
        @post_smearing_efficiencies = copy(%post_app_RHC_aNC_bkg)
>

channel(#NuBarMode_DS_app_tau_e)<
        @channel =      #DSNuBarMode:     +:      t:      e:      #CC:    #app_nue_bkg
       @post_smearing_efficiencies = copy(%post_app_RHC_nue_bkg)
>

channel(#NuBarMode_DS_app_tau_ebar)<
        @channel =      #DSNuBarMode_minus:     -:      t:      e:      #CC:    #app_nue_bkg
       @post_smearing_efficiencies = copy(%post_app_RHC_nue_bkg)
>

channel(#NuBarMode_DS_app_tau_numu)<
        @channel =      #DSNuBarMode:     +:      t:      m:      #CC:    #app_numu_bkg
       @post_smearing_efficiencies = copy(%post_app_RHC_nue_bkg)
>

channel(#NuBarMode_DS_app_tau_numubar)<
        @channel =      #DSNuBarMode_minus:     -:      t:      m:      #CC:    #app_numu_bkg
       @post_smearing_efficiencies = copy(%post_app_RHC_nue_bkg)
>

channel(#NuBarMode_DS_app_bkg_nutau)<
        @channel =      #DSNuBarMode:     +:      m:      t:      #CCtau:    #app_nutau_bkg
       @post_smearing_efficiencies = copy(%post_app_RHC_nutau_bkg)
>

channel(#NuBarMode_DS_app_bkg_nutaubar)<
         @channel =      #DSNuBarMode_minus:     -:      m:      t:      #CCtau:    #app_nutaubar_bkg
         @post_smearing_efficiencies = copy(%post_app_RHC_nutaubar_bkg)
>


channel(#NuBarMode_DS_dis_nutau)<
        @channel =      #DSNuBarMode:     +:      t:      t:      #CCtau:    #app_nutau_bkg
       @post_smearing_efficiencies = copy(%post_app_RHC_nutau_bkg)
>

channel(#NuBarMode_DS_dis_nutaubar)<
         @channel =      #DSNuBarMode_minus:     -:      t:      t:      #CCtau:    #app_nutaubar_bkg
         @post_smearing_efficiencies = copy(%post_app_RHC_nutaubar_bkg)
>


channel(#NuBarMode_DS_app_bkg_e_numu)<
        @channel =      #DSNuBarMode:     +:      e:      m:      #CC:    #app_numu_bkg
       @post_smearing_efficiencies = copy(%post_app_RHC_numu_bkg)
>

channel(#NuBarMode_DS_app_bkg_e_numubar)<
         @channel =      #DSNuBarMode_minus:     -:      e:      m:      #CC:    #app_numubar_bkg
         @post_smearing_efficiencies = copy(%post_app_RHC_numubar_bkg)
>

channel(#NuBarMode_DS_app_mu_nutau)<
        @channel =      #DSNuBarMode:     +:      m:      t:      #CCtau:    #tau_hadronicwithbias
>

channel(#NuBarMode_DS_app_mu_nutaubar)<
         @channel =      #DSNuBarMode_minus:     -:      m:      t:      #CCtau:    #tau_hadronicwithbias
>

channel(#NuBarMode_DS_app_tau_nutau)<
        @channel =      #DSNuBarMode:     +:      t:      t:      #CCtau:    #tau_hadronicwithbias
>

channel(#NuBarMode_DS_app_tau_nutaubar)<
         @channel =      #DSNuBarMode_minus:     -:      t:      t:      #CCtau:    #tau_hadronicwithbias
>

channel(#NuBarMode_DS_app_e_nutau)<
        @channel =      #DSNuBarMode:     +:      e:      t:      #CCtau:    #tau_hadronicwithbias
>

channel(#NuBarMode_DS_app_e_nutaubar)<
         @channel =      #DSNuBarMode_minus:     -:      e:      t:      #CCtau:    #tau_hadronicwithbias
>

channel(#NuBarMode_DS_app_bkg_nuNC_noeff)<
        @channel =      #DSNuBarMode:     +:      m: NOSC_m:      #NC:    #app_NC_bkg
>

channel(#NuBarMode_DS_app_bkg_nubarNC_noeff)<
        @channel =      #DSNuBarMode_minus:     -:      m: NOSC_m:      #NC:    #app_NC_bkg
>

channel(#NuBarMode_DS_app_bkg_nuNC_e_noeff)<
        @channel =      #DSNuBarMode:     +:      e: NOSC_e:      #NC:    #app_NC_bkg
>

channel(#NuBarMode_DS_app_bkg_nubarNC_e_noeff)<
        @channel =      #DSNuBarMode_minus:     -:      e: NOSC_e:      #NC:    #app_NC_bkg
>

channel(#NuBarMode_DS_app_bkg_nuNC_t_noeff)<
        @channel =      #DSNuBarMode:     +:      t: NOSC_t:      #NC:    #app_NC_bkg
>

channel(#NuBarMode_DS_app_bkg_nubarNC_t_noeff)<
        @channel =      #DSNuBarMode_minus:     -:      t: NOSC_t:      #NC:    #app_NC_bkg
>




//K

channel(#NuBarMode_Kaon_app_osc_nue)<
        @channel =      #KaonNuBarMode:     +:      m:      e:      #CC:    #app_nue_sig
       @post_smearing_efficiencies = copy(%post_app_RHC_nue_sig)
>

channel(#NuBarMode_Kaon_app_osc_nuebar)<
         @channel =      #KaonNuBarMode_minus:     -:      m:      e:      #CC:    #app_nuebar_sig
         @post_smearing_efficiencies = copy(%post_app_RHC_nuebar_sig)
>

channel(#NuBarMode_Kaon_app_bkg_nue)<
        @channel =      #KaonNuBarMode:     +:      e:      e:      #CC:    #app_nue_bkg
       @post_smearing_efficiencies = copy(%post_app_RHC_nue_bkg)
>

channel(#NuBarMode_Kaon_app_bkg_nuebar)<
         @channel =      #KaonNuBarMode_minus:     -:      e:      e:      #CC:    #app_nuebar_bkg
         @post_smearing_efficiencies = copy(%post_app_RHC_nuebar_bkg)
>

channel(#NuBarMode_Kaon_app_bkg_numu)<
        @channel =      #KaonNuBarMode:     +:      m:      m:      #CC:    #app_numu_bkg
       @post_smearing_efficiencies = copy(%post_app_RHC_numu_bkg)
>

channel(#NuBarMode_Kaon_app_bkg_numubar)<
         @channel =      #KaonNuBarMode_minus:     -:      m:      m:      #CC:    #app_numubar_bkg
         @post_smearing_efficiencies = copy(%post_app_RHC_numubar_bkg)
>


channel(#NuBarMode_Kaon_app_bkg_e_nutau)<
        @channel =      #KaonNuBarMode:     +:      e:      t:      #CCtau:    #app_nutau_bkg
       @post_smearing_efficiencies = copy(%post_app_RHC_nutau_bkg)
>

channel(#NuBarMode_Kaon_app_bkg_e_nutaubar)<
         @channel =      #KaonNuBarMode_minus:     -:      e:      t:      #CCtau:    #app_nutaubar_bkg
         @post_smearing_efficiencies = copy(%post_app_RHC_nutaubar_bkg)
>



channel(#NuBarMode_Kaon_app_bkg_nuNC)<
        @channel =      #KaonNuBarMode:     +:      m: NOSC_m:      #NC:    #app_NC_bkg
        @post_smearing_efficiencies = copy(%post_app_RHC_NC_bkg)
>

channel(#NuBarMode_Kaon_app_bkg_nubarNC)<
        @channel =      #KaonNuBarMode_minus:     -:      m: NOSC_m:      #NC:    #app_NC_bkg
        @post_smearing_efficiencies = copy(%post_app_RHC_aNC_bkg)
>

channel(#NuBarMode_Kaon_app_bkg_nuNC_e)<
        @channel =      #KaonNuBarMode:     +:      e: NOSC_e:      #NC:    #app_NC_bkg
        @post_smearing_efficiencies = copy(%post_app_RHC_NC_bkg)
>

channel(#NuBarMode_Kaon_app_bkg_nubarNC_e)<
        @channel =      #KaonNuBarMode_minus:     -:      e: NOSC_e:      #NC:    #app_NC_bkg
        @post_smearing_efficiencies = copy(%post_app_RHC_aNC_bkg)
>

channel(#NuBarMode_Kaon_app_bkg_nuNC_t)<
        @channel =      #KaonNuBarMode:     +:      t: NOSC_t:      #NC:    #app_NC_bkg
        @post_smearing_efficiencies = copy(%post_app_RHC_NC_bkg)
>

channel(#NuBarMode_Kaon_app_bkg_nubarNC_t)<
        @channel =      #KaonNuBarMode_minus:     -:      t: NOSC_t:      #NC:    #app_NC_bkg
        @post_smearing_efficiencies = copy(%post_app_RHC_aNC_bkg)
>

channel(#NuBarMode_Kaon_app_tau_e)<
        @channel =      #KaonNuBarMode:     +:      t:      e:      #CC:    #app_nue_bkg
       @post_smearing_efficiencies = copy(%post_app_RHC_nue_bkg)
>

channel(#NuBarMode_Kaon_app_tau_ebar)<
        @channel =      #KaonNuBarMode_minus:     -:      t:      e:      #CC:    #app_nue_bkg
       @post_smearing_efficiencies = copy(%post_app_RHC_nue_bkg)
>

channel(#NuBarMode_Kaon_app_tau_numu)<
        @channel =      #KaonNuBarMode:     +:      t:      m:      #CC:    #app_numu_bkg
       @post_smearing_efficiencies = copy(%post_app_RHC_nue_bkg)
>

channel(#NuBarMode_Kaon_app_tau_numubar)<
        @channel =      #KaonNuBarMode_minus:     -:      t:      m:      #CC:    #app_numu_bkg
       @post_smearing_efficiencies = copy(%post_app_RHC_nue_bkg)
>

channel(#NuBarMode_Kaon_app_bkg_nutau)<
        @channel =      #KaonNuBarMode:     +:      m:      t:      #CCtau:    #app_nutau_bkg
       @post_smearing_efficiencies = copy(%post_app_RHC_nutau_bkg)
>

channel(#NuBarMode_Kaon_app_bkg_nutaubar)<
         @channel =      #KaonNuBarMode_minus:     -:      m:      t:      #CCtau:    #app_nutaubar_bkg
         @post_smearing_efficiencies = copy(%post_app_RHC_nutaubar_bkg)
>


channel(#NuBarMode_Kaon_dis_nutau)<
        @channel =      #KaonNuBarMode:     +:      t:      t:      #CCtau:    #app_nutau_bkg
       @post_smearing_efficiencies = copy(%post_app_RHC_nutau_bkg)
>

channel(#NuBarMode_Kaon_dis_nutaubar)<
         @channel =      #KaonNuBarMode_minus:     -:      t:      t:      #CCtau:    #app_nutaubar_bkg
         @post_smearing_efficiencies = copy(%post_app_RHC_nutaubar_bkg)
>


channel(#NuBarMode_Kaon_app_bkg_e_numu)<
        @channel =      #KaonNuBarMode:     +:      e:      m:      #CC:    #app_numu_bkg
       @post_smearing_efficiencies = copy(%post_app_RHC_numu_bkg)
>

channel(#NuBarMode_Kaon_app_bkg_e_numubar)<
         @channel =      #KaonNuBarMode_minus:     -:      e:      m:      #CC:    #app_numubar_bkg
         @post_smearing_efficiencies = copy(%post_app_RHC_numubar_bkg)
>

channel(#NuBarMode_Kaon_app_mu_nutau)<
        @channel =      #KaonNuBarMode:     +:      m:      t:      #CCtau:    #tau_hadronicwithbias
>

channel(#NuBarMode_Kaon_app_mu_nutaubar)<
         @channel =      #KaonNuBarMode_minus:     -:      m:      t:      #CCtau:    #tau_hadronicwithbias
>

channel(#NuBarMode_Kaon_app_tau_nutau)<
        @channel =      #KaonNuBarMode:     +:      t:      t:      #CCtau:    #tau_hadronicwithbias
>

channel(#NuBarMode_Kaon_app_tau_nutaubar)<
         @channel =      #KaonNuBarMode_minus:     -:      t:      t:      #CCtau:    #tau_hadronicwithbias
>

channel(#NuBarMode_Kaon_app_e_nutau)<
        @channel =      #KaonNuBarMode:     +:      e:      t:      #CCtau:    #tau_hadronicwithbias
>

channel(#NuBarMode_Kaon_app_e_nutaubar)<
         @channel =      #KaonNuBarMode_minus:     -:      e:      t:      #CCtau:    #tau_hadronicwithbias
>

channel(#NuBarMode_Kaon_app_bkg_nuNC_noeff)<
        @channel =      #KaonNuBarMode:     +:      m: NOSC_m:      #NC:    #app_NC_bkg
>

channel(#NuBarMode_Kaon_app_bkg_nubarNC_noeff)<
        @channel =      #KaonNuBarMode_minus:     -:      m: NOSC_m:      #NC:    #app_NC_bkg
>

channel(#NuBarMode_Kaon_app_bkg_nuNC_e_noeff)<
        @channel =      #KaonNuBarMode:     +:      e: NOSC_e:      #NC:    #app_NC_bkg
>

channel(#NuBarMode_Kaon_app_bkg_nubarNC_e_noeff)<
        @channel =      #KaonNuBarMode_minus:     -:      e: NOSC_e:      #NC:    #app_NC_bkg
>

channel(#NuBarMode_Kaon_app_bkg_nuNC_t_noeff)<
        @channel =      #KaonNuBarMode:     +:      t: NOSC_t:      #NC:    #app_NC_bkg
>

channel(#NuBarMode_Kaon_app_bkg_nubarNC_t_noeff)<
        @channel =      #KaonNuBarMode_minus:     -:      t: NOSC_t:      #NC:    #app_NC_bkg
>





//Pion

channel(#NuBarMode_Pion_app_osc_nue)<
        @channel =      #PionNuBarMode:     +:      m:      e:      #CC:    #app_nue_sig
       @post_smearing_efficiencies = copy(%post_app_RHC_nue_sig)
>

channel(#NuBarMode_Pion_app_osc_nuebar)<
         @channel =      #PionNuBarMode_minus:     -:      m:      e:      #CC:    #app_nuebar_sig
         @post_smearing_efficiencies = copy(%post_app_RHC_nuebar_sig)
>

channel(#NuBarMode_Pion_app_bkg_nue)<
        @channel =      #PionNuBarMode:     +:      e:      e:      #CC:    #app_nue_bkg
       @post_smearing_efficiencies = copy(%post_app_RHC_nue_bkg)
>

channel(#NuBarMode_Pion_app_bkg_nuebar)<
         @channel =      #PionNuBarMode_minus:     -:      e:      e:      #CC:    #app_nuebar_bkg
         @post_smearing_efficiencies = copy(%post_app_RHC_nuebar_bkg)
>

channel(#NuBarMode_Pion_app_bkg_numu)<
        @channel =      #PionNuBarMode:     +:      m:      m:      #CC:    #app_numu_bkg
       @post_smearing_efficiencies = copy(%post_app_RHC_numu_bkg)
>

channel(#NuBarMode_Pion_app_bkg_numubar)<
         @channel =      #PionNuBarMode_minus:     -:      m:      m:      #CC:    #app_numubar_bkg
         @post_smearing_efficiencies = copy(%post_app_RHC_numubar_bkg)
>


channel(#NuBarMode_Pion_app_bkg_e_nutau)<
        @channel =      #PionNuBarMode:     +:      e:      t:      #CCtau:    #app_nutau_bkg
       @post_smearing_efficiencies = copy(%post_app_RHC_nutau_bkg)
>

channel(#NuBarMode_Pion_app_bkg_e_nutaubar)<
         @channel =      #PionNuBarMode_minus:     -:      e:      t:      #CCtau:    #app_nutaubar_bkg
         @post_smearing_efficiencies = copy(%post_app_RHC_nutaubar_bkg)
>



channel(#NuBarMode_Pion_app_bkg_nuNC)<
        @channel =      #PionNuBarMode:     +:      m: NOSC_m:      #NC:    #app_NC_bkg
        @post_smearing_efficiencies = copy(%post_app_RHC_NC_bkg)
>

channel(#NuBarMode_Pion_app_bkg_nubarNC)<
        @channel =      #PionNuBarMode_minus:     -:      m: NOSC_m:      #NC:    #app_NC_bkg
        @post_smearing_efficiencies = copy(%post_app_RHC_aNC_bkg)
>

channel(#NuBarMode_Pion_app_bkg_nuNC_e)<
        @channel =      #PionNuBarMode:     +:      e: NOSC_e:      #NC:    #app_NC_bkg
        @post_smearing_efficiencies = copy(%post_app_RHC_NC_bkg)
>

channel(#NuBarMode_Pion_app_bkg_nubarNC_e)<
        @channel =      #PionNuBarMode_minus:     -:      e: NOSC_e:      #NC:    #app_NC_bkg
        @post_smearing_efficiencies = copy(%post_app_RHC_aNC_bkg)
>

channel(#NuBarMode_Pion_app_bkg_nuNC_t)<
        @channel =      #PionNuBarMode:     +:      t: NOSC_t:      #NC:    #app_NC_bkg
        @post_smearing_efficiencies = copy(%post_app_RHC_NC_bkg)
>

channel(#NuBarMode_Pion_app_bkg_nubarNC_t)<
        @channel =      #PionNuBarMode_minus:     -:      t: NOSC_t:      #NC:    #app_NC_bkg
        @post_smearing_efficiencies = copy(%post_app_RHC_aNC_bkg)
>

channel(#NuBarMode_Pion_app_tau_e)<
        @channel =      #PionNuBarMode:     +:      t:      e:      #CC:    #app_nue_bkg
       @post_smearing_efficiencies = copy(%post_app_RHC_nue_bkg)
>

channel(#NuBarMode_Pion_app_tau_ebar)<
        @channel =      #PionNuBarMode_minus:     -:      t:      e:      #CC:    #app_nue_bkg
       @post_smearing_efficiencies = copy(%post_app_RHC_nue_bkg)
>

channel(#NuBarMode_Pion_app_tau_numu)<
        @channel =      #PionNuBarMode:     +:      t:      m:      #CC:    #app_numu_bkg
       @post_smearing_efficiencies = copy(%post_app_RHC_nue_bkg)
>

channel(#NuBarMode_Pion_app_tau_numubar)<
        @channel =      #PionNuBarMode_minus:     -:      t:      m:      #CC:    #app_numu_bkg
       @post_smearing_efficiencies = copy(%post_app_RHC_nue_bkg)
>

channel(#NuBarMode_Pion_app_bkg_nutau)<
        @channel =      #PionNuBarMode:     +:      m:      t:      #CCtau:    #app_nutau_bkg
       @post_smearing_efficiencies = copy(%post_app_RHC_nutau_bkg)
>

channel(#NuBarMode_Pion_app_bkg_nutaubar)<
         @channel =      #PionNuBarMode_minus:     -:      m:      t:      #CCtau:    #app_nutaubar_bkg
         @post_smearing_efficiencies = copy(%post_app_RHC_nutaubar_bkg)
>


channel(#NuBarMode_Pion_dis_nutau)<
        @channel =      #PionNuBarMode:     +:      t:      t:      #CCtau:    #app_nutau_bkg
       @post_smearing_efficiencies = copy(%post_app_RHC_nutau_bkg)
>

channel(#NuBarMode_Pion_dis_nutaubar)<
         @channel =      #PionNuBarMode_minus:     -:      t:      t:      #CCtau:    #app_nutaubar_bkg
         @post_smearing_efficiencies = copy(%post_app_RHC_nutaubar_bkg)
>


channel(#NuBarMode_Pion_app_bkg_e_numu)<
        @channel =      #PionNuBarMode:     +:      e:      m:      #CC:    #app_numu_bkg
       @post_smearing_efficiencies = copy(%post_app_RHC_numu_bkg)
>

channel(#NuBarMode_Pion_app_bkg_e_numubar)<
         @channel =      #PionNuBarMode_minus:     -:      e:      m:      #CC:    #app_numubar_bkg
         @post_smearing_efficiencies = copy(%post_app_RHC_numubar_bkg)
>

channel(#NuBarMode_Pion_app_mu_nutau)<
        @channel =      #PionNuBarMode:     +:      m:      t:      #CCtau:    #tau_hadronicwithbias
>

channel(#NuBarMode_Pion_app_mu_nutaubar)<
         @channel =      #PionNuBarMode_minus:     -:      m:      t:      #CCtau:    #tau_hadronicwithbias
>

channel(#NuBarMode_Pion_app_tau_nutau)<
        @channel =      #PionNuBarMode:     +:      t:      t:      #CCtau:    #tau_hadronicwithbias
>

channel(#NuBarMode_Pion_app_tau_nutaubar)<
         @channel =      #PionNuBarMode_minus:     -:      t:      t:      #CCtau:    #tau_hadronicwithbias
>

channel(#NuBarMode_Pion_app_e_nutau)<
        @channel =      #PionNuBarMode:     +:      e:      t:      #CCtau:    #tau_hadronicwithbias
>

channel(#NuBarMode_Pion_app_e_nutaubar)<
         @channel =      #PionNuBarMode_minus:     -:      e:      t:      #CCtau:    #tau_hadronicwithbias
>

channel(#NuBarMode_Pion_app_bkg_nuNC_noeff)<
        @channel =      #PionNuBarMode:     +:      m: NOSC_m:      #NC:    #app_NC_bkg
>

channel(#NuBarMode_Pion_app_bkg_nubarNC_noeff)<
        @channel =      #PionNuBarMode_minus:     -:      m: NOSC_m:      #NC:    #app_NC_bkg
>

channel(#NuBarMode_Pion_app_bkg_nuNC_e_noeff)<
        @channel =      #PionNuBarMode:     +:      e: NOSC_e:      #NC:    #app_NC_bkg
>

channel(#NuBarMode_Pion_app_bkg_nubarNC_e_noeff)<
        @channel =      #PionNuBarMode_minus:     -:      e: NOSC_e:      #NC:    #app_NC_bkg
>

channel(#NuBarMode_Pion_app_bkg_nuNC_t_noeff)<
        @channel =      #PionNuBarMode:     +:      t: NOSC_t:      #NC:    #app_NC_bkg
>

channel(#NuBarMode_Pion_app_bkg_nubarNC_t_noeff)<
        @channel =      #PionNuBarMode_minus:     -:      t: NOSC_t:      #NC:    #app_NC_bkg
>




//Muon


channel(#NuBarMode_Muon_app_osc_nue)<
        @channel =      #MuonNuBarMode:     +:      m:      e:      #CC:    #app_nue_sig
       @post_smearing_efficiencies = copy(%post_app_RHC_nue_sig)
>

channel(#NuBarMode_Muon_app_osc_nuebar)<
         @channel =      #MuonNuBarMode_minus:     -:      m:      e:      #CC:    #app_nuebar_sig
         @post_smearing_efficiencies = copy(%post_app_RHC_nuebar_sig)
>

channel(#NuBarMode_Muon_app_bkg_nue)<
        @channel =      #MuonNuBarMode:     +:      e:      e:      #CC:    #app_nue_bkg
       @post_smearing_efficiencies = copy(%post_app_RHC_nue_bkg)
>

channel(#NuBarMode_Muon_app_bkg_nuebar)<
         @channel =      #MuonNuBarMode_minus:     -:      e:      e:      #CC:    #app_nuebar_bkg
         @post_smearing_efficiencies = copy(%post_app_RHC_nuebar_bkg)
>

channel(#NuBarMode_Muon_app_bkg_numu)<
        @channel =      #MuonNuBarMode:     +:      m:      m:      #CC:    #app_numu_bkg
       @post_smearing_efficiencies = copy(%post_app_RHC_numu_bkg)
>

channel(#NuBarMode_Muon_app_bkg_numubar)<
         @channel =      #MuonNuBarMode_minus:     -:      m:      m:      #CC:    #app_numubar_bkg
         @post_smearing_efficiencies = copy(%post_app_RHC_numubar_bkg)
>


channel(#NuBarMode_Muon_app_bkg_e_nutau)<
        @channel =      #MuonNuBarMode:     +:      e:      t:      #CCtau:    #app_nutau_bkg
       @post_smearing_efficiencies = copy(%post_app_RHC_nutau_bkg)
>

channel(#NuBarMode_Muon_app_bkg_e_nutaubar)<
         @channel =      #MuonNuBarMode_minus:     -:      e:      t:      #CCtau:    #app_nutaubar_bkg
         @post_smearing_efficiencies = copy(%post_app_RHC_nutaubar_bkg)
>



channel(#NuBarMode_Muon_app_bkg_nuNC)<
        @channel =      #MuonNuBarMode:     +:      m: NOSC_m:      #NC:    #app_NC_bkg
        @post_smearing_efficiencies = copy(%post_app_RHC_NC_bkg)
>

channel(#NuBarMode_Muon_app_bkg_nubarNC)<
        @channel =      #MuonNuBarMode_minus:     -:      m: NOSC_m:      #NC:    #app_NC_bkg
        @post_smearing_efficiencies = copy(%post_app_RHC_aNC_bkg)
>

channel(#NuBarMode_Muon_app_bkg_nuNC_e)<
        @channel =      #MuonNuBarMode:     +:      e: NOSC_e:      #NC:    #app_NC_bkg
        @post_smearing_efficiencies = copy(%post_app_RHC_NC_bkg)
>

channel(#NuBarMode_Muon_app_bkg_nubarNC_e)<
        @channel =      #MuonNuBarMode_minus:     -:      e: NOSC_e:      #NC:    #app_NC_bkg
        @post_smearing_efficiencies = copy(%post_app_RHC_aNC_bkg)
>

channel(#NuBarMode_Muon_app_bkg_nuNC_t)<
        @channel =      #MuonNuBarMode:     +:      t: NOSC_t:      #NC:    #app_NC_bkg
        @post_smearing_efficiencies = copy(%post_app_RHC_NC_bkg)
>

channel(#NuBarMode_Muon_app_bkg_nubarNC_t)<
        @channel =      #MuonNuBarMode_minus:     -:      t: NOSC_t:      #NC:    #app_NC_bkg
        @post_smearing_efficiencies = copy(%post_app_RHC_aNC_bkg)
>

channel(#NuBarMode_Muon_app_tau_e)<
        @channel =      #MuonNuBarMode:     +:      t:      e:      #CC:    #app_nue_bkg
       @post_smearing_efficiencies = copy(%post_app_RHC_nue_bkg)
>

channel(#NuBarMode_Muon_app_tau_ebar)<
        @channel =      #MuonNuBarMode_minus:     -:      t:      e:      #CC:    #app_nue_bkg
       @post_smearing_efficiencies = copy(%post_app_RHC_nue_bkg)
>

channel(#NuBarMode_Muon_app_tau_numu)<
        @channel =      #MuonNuBarMode:     +:      t:      m:      #CC:    #app_numu_bkg
       @post_smearing_efficiencies = copy(%post_app_RHC_nue_bkg)
>

channel(#NuBarMode_Muon_app_tau_numubar)<
        @channel =      #MuonNuBarMode_minus:     -:      t:      m:      #CC:    #app_numu_bkg
       @post_smearing_efficiencies = copy(%post_app_RHC_nue_bkg)
>

channel(#NuBarMode_Muon_app_bkg_nutau)<
        @channel =      #MuonNuBarMode:     +:      m:      t:      #CCtau:    #app_nutau_bkg
       @post_smearing_efficiencies = copy(%post_app_RHC_nutau_bkg)
>

channel(#NuBarMode_Muon_app_bkg_nutaubar)<
         @channel =      #MuonNuBarMode_minus:     -:      m:      t:      #CCtau:    #app_nutaubar_bkg
         @post_smearing_efficiencies = copy(%post_app_RHC_nutaubar_bkg)
>


channel(#NuBarMode_Muon_dis_nutau)<
        @channel =      #MuonNuBarMode:     +:      t:      t:      #CCtau:    #app_nutau_bkg
       @post_smearing_efficiencies = copy(%post_app_RHC_nutau_bkg)
>

channel(#NuBarMode_Muon_dis_nutaubar)<
         @channel =      #MuonNuBarMode_minus:     -:      t:      t:      #CCtau:    #app_nutaubar_bkg
         @post_smearing_efficiencies = copy(%post_app_RHC_nutaubar_bkg)
>


channel(#NuBarMode_Muon_app_bkg_e_numu)<
        @channel =      #MuonNuBarMode:     +:      e:      m:      #CC:    #app_numu_bkg
       @post_smearing_efficiencies = copy(%post_app_RHC_numu_bkg)
>

channel(#NuBarMode_Muon_app_bkg_e_numubar)<
         @channel =      #MuonNuBarMode_minus:     -:      e:      m:      #CC:    #app_numubar_bkg
         @post_smearing_efficiencies = copy(%post_app_RHC_numubar_bkg)
>

channel(#NuBarMode_Muon_app_mu_nutau)<
        @channel =      #MuonNuBarMode:     +:      m:      t:      #CCtau:    #tau_hadronicwithbias
>

channel(#NuBarMode_Muon_app_mu_nutaubar)<
         @channel =      #MuonNuBarMode_minus:     -:      m:      t:      #CCtau:    #tau_hadronicwithbias
>

channel(#NuBarMode_Muon_app_tau_nutau)<
        @channel =      #MuonNuBarMode:     +:      t:      t:      #CCtau:    #tau_hadronicwithbias
>

channel(#NuBarMode_Muon_app_tau_nutaubar)<
         @channel =      #MuonNuBarMode_minus:     -:      t:      t:      #CCtau:    #tau_hadronicwithbias
>

channel(#NuBarMode_Muon_app_e_nutau)<
        @channel =      #MuonNuBarMode:     +:      e:      t:      #CCtau:    #tau_hadronicwithbias
>

channel(#NuBarMode_Muon_app_e_nutaubar)<
         @channel =      #MuonNuBarMode_minus:     -:      e:      t:      #CCtau:    #tau_hadronicwithbias
>

channel(#NuBarMode_Muon_app_bkg_nuNC_noeff)<
        @channel =      #MuonNuBarMode:     +:      m: NOSC_m:      #NC:    #app_NC_bkg
>

channel(#NuBarMode_Muon_app_bkg_nubarNC_noeff)<
        @channel =      #MuonNuBarMode_minus:     -:      m: NOSC_m:      #NC:    #app_NC_bkg
>

channel(#NuBarMode_Muon_app_bkg_nuNC_e_noeff)<
        @channel =      #MuonNuBarMode:     +:      e: NOSC_e:      #NC:    #app_NC_bkg
>

channel(#NuBarMode_Muon_app_bkg_nubarNC_e_noeff)<
        @channel =      #MuonNuBarMode_minus:     -:      e: NOSC_e:      #NC:    #app_NC_bkg
>

channel(#NuBarMode_Muon_app_bkg_nuNC_t_noeff)<
        @channel =      #MuonNuBarMode:     +:      t: NOSC_t:      #NC:    #app_NC_bkg
>

channel(#NuBarMode_Muon_app_bkg_nubarNC_t_noeff)<
        @channel =      #MuonNuBarMode_minus:     -:      t: NOSC_t:      #NC:    #app_NC_bkg
>









/* NUMU DIS */

//D

channel(#NuMode_Dplus_dis_sig_numu)<
        @channel = #DplusNuMode:     +:      m:      m:      #CC:    #dis_numu_sig
        @post_smearing_efficiencies = copy(%post_dis_FHC_numu_sig)
>

channel(#NuMode_Dplus_dis_sig_numubar)<
        @channel = #DplusNuMode_minus:     -:      m:      m:      #CC:    #dis_numubar_sig
        @post_smearing_efficiencies = copy(%post_dis_FHC_numubar_sig)
>

channel(#NuMode_Dplus_dis_bkg_nutau)<
        @channel = #DplusNuMode:     +:      m:      t:      #CCtau:    #dis_nutau_bkg
        @post_smearing_efficiencies = copy(%post_dis_FHC_nutau_bkg)
>

channel(#NuMode_Dplus_dis_bkg_nutaubar)<
        @channel = #DplusNuMode_minus:     -:      m:      t:      #CCtau:    #dis_nutaubar_bkg
        @post_smearing_efficiencies = copy(%post_dis_FHC_nutaubar_bkg)
>

channel(#NuMode_Dplus_dis_bkg_nuNC)<
        @channel =      #DplusNuMode:     +:      m: NOSC_m:      #NC:    #dis_NC_bkg
        @post_smearing_efficiencies = copy(%post_dis_FHC_NC_bkg)
>

channel(#NuMode_Dplus_dis_bkg_nubarNC)<
        @channel =      #DplusNuMode_minus:     -:      m: NOSC_m:      #NC:    #dis_aNC_bkg
        @post_smearing_efficiencies = copy(%post_dis_FHC_NC_bkg)
>

channel(#NuMode_Dplus_dis_bkg_nuNC_e)<
        @channel =      #DplusNuMode:     +:      e: NOSC_e:      #NC:    #dis_NC_bkg
        @post_smearing_efficiencies = copy(%post_dis_FHC_NC_bkg)
>

channel(#NuMode_Dplus_dis_bkg_nubarNC_e)<
        @channel =      #DplusNuMode_minus:     -:      e: NOSC_e:      #NC:    #dis_aNC_bkg
        @post_smearing_efficiencies = copy(%post_dis_FHC_NC_bkg)
>

channel(#NuMode_Dplus_dis_bkg_nuNC_tau)<
        @channel =      #DplusNuMode:     +:      t: NOSC_t:      #NC:    #dis_NC_bkg
        @post_smearing_efficiencies = copy(%post_dis_FHC_NC_bkg)
>

channel(#NuMode_Dplus_dis_bkg_nubarNC_tau)<
        @channel =      #DplusNuMode_minus:     -:      t: NOSC_t:      #NC:    #dis_aNC_bkg
        @post_smearing_efficiencies = copy(%post_dis_FHC_NC_bkg)
>



channel(#NuMode_Dplus_dis_bkg_e_nutau)<
        @channel = #DplusNuMode:     +:      e:      t:      #CCtau:    #dis_nutau_bkg
        @post_smearing_efficiencies = copy(%post_dis_FHC_nutau_bkg)
>

channel(#NuMode_Dplus_dis_bkg_e_nutaubar)<
        @channel = #DplusNuMode_minus:     -:      e:      t:      #CCtau:    #dis_nutaubar_bkg
        @post_smearing_efficiencies = copy(%post_dis_FHC_nutaubar_bkg)
>

channel(#NuMode_Dplus_dis_sig_e_numu)<
        @channel = #DplusNuMode:     +:      e:      m:      #CC:    #dis_numu_sig
        @post_smearing_efficiencies = copy(%post_dis_FHC_numu_sig)
>

channel(#NuMode_Dplus_dis_sig_e_numubar)<
        @channel = #DplusNuMode_minus:     -:      e:      m:      #CC:    #dis_numubar_sig
        @post_smearing_efficiencies = copy(%post_dis_FHC_numubar_sig)
>

channel(#NuMode_Dplus_dis_sig_t_numu)<
        @channel = #DplusNuMode:     +:      t:      m:      #CC:    #dis_numu_sig
        @post_smearing_efficiencies = copy(%post_dis_FHC_numu_sig)
>

channel(#NuMode_Dplus_dis_sig_t_numubar)<
        @channel = #DplusNuMode_minus:     -:      t:      m:      #CC:    #dis_numubar_sig
        @post_smearing_efficiencies = copy(%post_dis_FHC_numubar_sig)
>


//Ds

channel(#NuMode_DS_dis_sig_numu)<
        @channel =      #DSNuMode:     +:      m:      m:      #CC:    #dis_numu_sig
        @post_smearing_efficiencies = copy(%post_dis_FHC_numu_sig)
>

channel(#NuMode_DS_dis_sig_numubar)<
        @channel =      #DSNuMode_minus:     -:      m:      m:      #CC:    #dis_numubar_sig
        @post_smearing_efficiencies = copy(%post_dis_FHC_numubar_sig)
>

channel(#NuMode_DS_dis_bkg_nutau)<
        @channel =      #DSNuMode:     +:      m:      t:      #CCtau:    #dis_nutau_bkg
        @post_smearing_efficiencies = copy(%post_dis_FHC_nutau_bkg)
>

channel(#NuMode_DS_dis_bkg_nutaubar)<
        @channel =      #DSNuMode_minus:     -:      m:      t:      #CCtau:    #dis_nutaubar_bkg
        @post_smearing_efficiencies = copy(%post_dis_FHC_nutaubar_bkg)
>

channel(#NuMode_DS_dis_bkg_nuNC)<
        @channel =      #DSNuMode:     +:      m: NOSC_m:      #NC:    #dis_NC_bkg
        @post_smearing_efficiencies = copy(%post_dis_FHC_NC_bkg)
>

channel(#NuMode_DS_dis_bkg_nubarNC)<
        @channel =      #DSNuMode_minus:     -:      m: NOSC_m:      #NC:    #dis_aNC_bkg
        @post_smearing_efficiencies = copy(%post_dis_FHC_NC_bkg)
>

channel(#NuMode_DS_dis_bkg_nuNC_e)<
        @channel =      #DSNuMode:     +:      e: NOSC_e:      #NC:    #dis_NC_bkg
        @post_smearing_efficiencies = copy(%post_dis_FHC_NC_bkg)
>

channel(#NuMode_DS_dis_bkg_nubarNC_e)<
        @channel =      #DSNuMode_minus:     -:      e: NOSC_e:      #NC:    #dis_aNC_bkg
        @post_smearing_efficiencies = copy(%post_dis_FHC_NC_bkg)
>

channel(#NuMode_DS_dis_bkg_nuNC_tau)<
        @channel =      #DSNuMode:     +:      t: NOSC_t:      #NC:    #dis_NC_bkg
        @post_smearing_efficiencies = copy(%post_dis_FHC_NC_bkg)
>

channel(#NuMode_DS_dis_bkg_nubarNC_tau)<
        @channel =      #DSNuMode_minus:     -:      t: NOSC_t:      #NC:    #dis_aNC_bkg
        @post_smearing_efficiencies = copy(%post_dis_FHC_NC_bkg)
>

channel(#NuMode_DS_dis_bkg_e_nutau)<
        @channel = #DSNuMode:     +:      e:      t:      #CCtau:    #dis_nutau_bkg
        @post_smearing_efficiencies = copy(%post_dis_FHC_nutau_bkg)
>

channel(#NuMode_DS_dis_bkg_e_nutaubar)<
        @channel = #DSNuMode_minus:     -:      e:      t:      #CCtau:    #dis_nutaubar_bkg
        @post_smearing_efficiencies = copy(%post_dis_FHC_nutaubar_bkg)
>

channel(#NuMode_DS_dis_sig_e_numu)<
        @channel = #DSNuMode:     +:      e:      m:      #CC:    #dis_numu_sig
        @post_smearing_efficiencies = copy(%post_dis_FHC_numu_sig)
>

channel(#NuMode_DS_dis_sig_e_numubar)<
        @channel = #DSNuMode_minus:     -:      e:      m:      #CC:    #dis_numubar_sig
        @post_smearing_efficiencies = copy(%post_dis_FHC_numubar_sig)
>

channel(#NuMode_DS_dis_sig_t_numu)<
        @channel = #DSNuMode:     +:      t:      m:      #CC:    #dis_numu_sig
        @post_smearing_efficiencies = copy(%post_dis_FHC_numu_sig)
>

channel(#NuMode_DS_dis_sig_t_numubar)<
        @channel = #DSNuMode_minus:     -:      t:      m:      #CC:    #dis_numubar_sig
        @post_smearing_efficiencies = copy(%post_dis_FHC_numubar_sig)
>


//Kaon

channel(#NuMode_Kaon_dis_sig_numu)<
        @channel =      #KaonNuMode:     +:      m:      m:      #CC:    #dis_numu_sig
        @post_smearing_efficiencies = copy(%post_dis_FHC_numu_sig)
>

channel(#NuMode_Kaon_dis_sig_numubar)<
        @channel =      #KaonNuMode_minus:     -:      m:      m:      #CC:    #dis_numubar_sig
        @post_smearing_efficiencies = copy(%post_dis_FHC_numubar_sig)
>

channel(#NuMode_Kaon_dis_bkg_nutau)<
        @channel =      #KaonNuMode:     +:      m:      t:      #CCtau:    #dis_nutau_bkg
        @post_smearing_efficiencies = copy(%post_dis_FHC_nutau_bkg)
>

channel(#NuMode_Kaon_dis_bkg_nutaubar)<
        @channel =      #KaonNuMode:      -:      m:      t:      #CCtau:    #dis_nutaubar_bkg
        @post_smearing_efficiencies = copy(%post_dis_FHC_nutaubar_bkg)
>

channel(#NuMode_Kaon_dis_bkg_nuNC)<
        @channel =      #KaonNuMode:     +:      m: NOSC_m:      #NC:    #dis_NC_bkg
        @post_smearing_efficiencies = copy(%post_dis_FHC_NC_bkg)
>

channel(#NuMode_Kaon_dis_bkg_nubarNC)<
        @channel =      #KaonNuMode_minus:     -:      m: NOSC_m:      #NC:    #dis_aNC_bkg
        @post_smearing_efficiencies = copy(%post_dis_FHC_NC_bkg)
>

channel(#NuMode_Kaon_dis_bkg_nuNC_e)<
        @channel =      #KaonNuMode:     +:      e: NOSC_e:      #NC:    #dis_NC_bkg
        @post_smearing_efficiencies = copy(%post_dis_FHC_NC_bkg)
>

channel(#NuMode_Kaon_dis_bkg_nubarNC_e)<
        @channel =      #KaonNuMode_minus:     -:      e: NOSC_e:      #NC:    #dis_aNC_bkg
        @post_smearing_efficiencies = copy(%post_dis_FHC_NC_bkg)
>

channel(#NuMode_Kaon_dis_bkg_nuNC_tau)<
        @channel =      #KaonNuMode:     +:      t: NOSC_t:      #NC:    #dis_NC_bkg
        @post_smearing_efficiencies = copy(%post_dis_FHC_NC_bkg)
>

channel(#NuMode_Kaon_dis_bkg_nubarNC_tau)<
        @channel =      #KaonNuMode_minus:     -:      t: NOSC_t:      #NC:    #dis_aNC_bkg
        @post_smearing_efficiencies = copy(%post_dis_FHC_NC_bkg)
>


channel(#NuMode_Kaon_dis_bkg_e_nutau)<
        @channel = #KaonNuMode:     +:      e:      t:      #CCtau:    #dis_nutau_bkg
        @post_smearing_efficiencies = copy(%post_dis_FHC_nutau_bkg)
>

channel(#NuMode_Kaon_dis_bkg_e_nutaubar)<
        @channel = #KaonNuMode_minus:     -:      e:      t:      #CCtau:    #dis_nutaubar_bkg
        @post_smearing_efficiencies = copy(%post_dis_FHC_nutaubar_bkg)
>

channel(#NuMode_Kaon_dis_sig_e_numu)<
        @channel = #KaonNuMode:     +:      e:      m:      #CC:    #dis_numu_sig
        @post_smearing_efficiencies = copy(%post_dis_FHC_numu_sig)
>

channel(#NuMode_Kaon_dis_sig_e_numubar)<
        @channel = #KaonNuMode_minus:     -:      e:      m:      #CC:    #dis_numubar_sig
        @post_smearing_efficiencies = copy(%post_dis_FHC_numubar_sig)
>

channel(#NuMode_Kaon_dis_sig_t_numu)<
        @channel = #KaonNuMode:     +:      t:      m:      #CC:    #dis_numu_sig
        @post_smearing_efficiencies = copy(%post_dis_FHC_numu_sig)
>

channel(#NuMode_Kaon_dis_sig_t_numubar)<
        @channel = #KaonNuMode_minus:     -:      t:      m:      #CC:    #dis_numubar_sig
        @post_smearing_efficiencies = copy(%post_dis_FHC_numubar_sig)
>


//Pion

channel(#NuMode_Pion_dis_sig_numu)<
        @channel =      #PionNuMode:     +:      m:      m:      #CC:    #dis_numu_sig
        @post_smearing_efficiencies = copy(%post_dis_FHC_numu_sig)
>

channel(#NuMode_Pion_dis_sig_numubar)<
        @channel =      #PionNuMode_minus:     -:      m:      m:      #CC:    #dis_numubar_sig
        @post_smearing_efficiencies = copy(%post_dis_FHC_numubar_sig)
>

channel(#NuMode_Pion_dis_bkg_nutau)<
        @channel =      #PionNuMode:     +:      m:      t:      #CCtau:    #dis_nutau_bkg
        @post_smearing_efficiencies = copy(%post_dis_FHC_nutau_bkg)
>

channel(#NuMode_Pion_dis_bkg_nutaubar)<
        @channel =      #PionNuMode_minus:     -:      m:      t:      #CCtau:   #dis_nutaubar_bkg
        @post_smearing_efficiencies = copy(%post_dis_FHC_nutaubar_bkg)
>

channel(#NuMode_Pion_dis_bkg_nuNC)<
        @channel =      #PionNuMode:     +:      m: NOSC_m:      #NC:    #dis_NC_bkg
        @post_smearing_efficiencies = copy(%post_dis_FHC_NC_bkg)
>

channel(#NuMode_Pion_dis_bkg_nubarNC)<
        @channel =      #PionNuMode_minus:     -:      m: NOSC_m:      #NC:    #dis_aNC_bkg
        @post_smearing_efficiencies = copy(%post_dis_FHC_NC_bkg)
>

channel(#NuMode_Pion_dis_bkg_nuNC_e)<
        @channel =      #PionNuMode:     +:      e: NOSC_e:      #NC:    #dis_NC_bkg
        @post_smearing_efficiencies = copy(%post_dis_FHC_NC_bkg)
>

channel(#NuMode_Pion_dis_bkg_nubarNC_e)<
        @channel =      #PionNuMode_minus:     -:      e: NOSC_e:      #NC:    #dis_aNC_bkg
        @post_smearing_efficiencies = copy(%post_dis_FHC_NC_bkg)
>

channel(#NuMode_Pion_dis_bkg_nuNC_tau)<
        @channel =      #PionNuMode:     +:      t: NOSC_t:      #NC:    #dis_NC_bkg
        @post_smearing_efficiencies = copy(%post_dis_FHC_NC_bkg)
>

channel(#NuMode_Pion_dis_bkg_nubarNC_tau)<
        @channel =      #PionNuMode_minus:     -:      t: NOSC_t:      #NC:    #dis_aNC_bkg
        @post_smearing_efficiencies = copy(%post_dis_FHC_NC_bkg)
>

channel(#NuMode_Pion_dis_bkg_e_nutau)<
        @channel = #PionNuMode:     +:      e:      t:      #CCtau:    #dis_nutau_bkg
        @post_smearing_efficiencies = copy(%post_dis_FHC_nutau_bkg)
>

channel(#NuMode_Pion_dis_bkg_e_nutaubar)<
        @channel = #PionNuMode_minus:     -:      e:      t:      #CCtau:    #dis_nutaubar_bkg
        @post_smearing_efficiencies = copy(%post_dis_FHC_nutaubar_bkg)
>

channel(#NuMode_Pion_dis_sig_e_numu)<
        @channel = #PionNuMode:     +:      e:      m:      #CC:    #dis_numu_sig
        @post_smearing_efficiencies = copy(%post_dis_FHC_numu_sig)
>

channel(#NuMode_Pion_dis_sig_e_numubar)<
        @channel = #PionNuMode_minus:     -:      e:      m:      #CC:    #dis_numubar_sig
        @post_smearing_efficiencies = copy(%post_dis_FHC_numubar_sig)
>

channel(#NuMode_Pion_dis_sig_t_numu)<
        @channel = #PionNuMode:     +:      t:      m:      #CC:    #dis_numu_sig
        @post_smearing_efficiencies = copy(%post_dis_FHC_numu_sig)
>

channel(#NuMode_Pion_dis_sig_t_numubar)<
        @channel = #PionNuMode_minus:     -:      t:      m:      #CC:    #dis_numubar_sig
        @post_smearing_efficiencies = copy(%post_dis_FHC_numubar_sig)
>

//Muon 

channel(#NuMode_Muon_dis_sig_numu)<
        @channel =      #MuonNuMode:     +:      m:      m:      #CC:    #dis_numu_sig
        @post_smearing_efficiencies = copy(%post_dis_FHC_numu_sig)
>

channel(#NuMode_Muon_dis_sig_numubar)<
        @channel =      #MuonNuMode_minus:     -:      m:      m:      #CC:    #dis_numubar_sig
        @post_smearing_efficiencies = copy(%post_dis_FHC_numubar_sig)
>

channel(#NuMode_Muon_dis_bkg_nutau)<
        @channel =      #MuonNuMode:     +:      m:      t:      #CCtau:    #dis_nutau_bkg
        @post_smearing_efficiencies = copy(%post_dis_FHC_nutau_bkg)
>

channel(#NuMode_Muon_dis_bkg_nutaubar)<
        @channel =      #MuonNuMode_minus:     -:      m:      t:      #CCtau:   #dis_nutaubar_bkg
        @post_smearing_efficiencies = copy(%post_dis_FHC_nutaubar_bkg)
>

channel(#NuMode_Muon_dis_bkg_nuNC)<
        @channel =      #MuonNuMode:     +:      m: NOSC_m:      #NC:    #dis_NC_bkg
        @post_smearing_efficiencies = copy(%post_dis_FHC_NC_bkg)
>

channel(#NuMode_Muon_dis_bkg_nubarNC)<
        @channel =      #MuonNuMode_minus:     -:      m: NOSC_m:      #NC:    #dis_aNC_bkg
        @post_smearing_efficiencies = copy(%post_dis_FHC_NC_bkg)
>

channel(#NuMode_Muon_dis_bkg_nuNC_e)<
        @channel =      #MuonNuMode:     +:      e: NOSC_e:      #NC:    #dis_NC_bkg
        @post_smearing_efficiencies = copy(%post_dis_FHC_NC_bkg)
>

channel(#NuMode_Muon_dis_bkg_nubarNC_e)<
        @channel =      #MuonNuMode_minus:     -:      e: NOSC_e:      #NC:    #dis_aNC_bkg
        @post_smearing_efficiencies = copy(%post_dis_FHC_NC_bkg)
>

channel(#NuMode_Muon_dis_bkg_nuNC_tau)<
        @channel =      #MuonNuMode:     +:      t: NOSC_t:      #NC:    #dis_NC_bkg
        @post_smearing_efficiencies = copy(%post_dis_FHC_NC_bkg)
>

channel(#NuMode_Muon_dis_bkg_nubarNC_tau)<
        @channel =      #MuonNuMode_minus:     -:      t: NOSC_t:      #NC:    #dis_aNC_bkg
        @post_smearing_efficiencies = copy(%post_dis_FHC_NC_bkg)
>



channel(#NuMode_Muon_dis_bkg_e_nutau)<
        @channel = #MuonNuMode:     +:      e:      t:      #CCtau:    #dis_nutau_bkg
        @post_smearing_efficiencies = copy(%post_dis_FHC_nutau_bkg)
>

channel(#NuMode_Muon_dis_bkg_e_nutaubar)<
        @channel = #MuonNuMode_minus:     -:      e:      t:      #CCtau:    #dis_nutaubar_bkg
        @post_smearing_efficiencies = copy(%post_dis_FHC_nutaubar_bkg)
>

channel(#NuMode_Muon_dis_sig_e_numu)<
        @channel = #MuonNuMode:     +:      e:      m:      #CC:    #dis_numu_sig
        @post_smearing_efficiencies = copy(%post_dis_FHC_numu_sig)
>

channel(#NuMode_Muon_dis_sig_e_numubar)<
        @channel = #MuonNuMode_minus:     -:      e:      m:      #CC:    #dis_numubar_sig
        @post_smearing_efficiencies = copy(%post_dis_FHC_numubar_sig)
>

channel(#NuMode_Muon_dis_sig_t_numu)<
        @channel = #MuonNuMode:     +:      t:      m:      #CC:    #dis_numu_sig
        @post_smearing_efficiencies = copy(%post_dis_FHC_numu_sig)
>

channel(#NuMode_Muon_dis_sig_t_numubar)<
        @channel = #MuonNuMode_minus:     -:      t:      m:      #CC:    #dis_numubar_sig
        @post_smearing_efficiencies = copy(%post_dis_FHC_numubar_sig)
>




/* NUMUBAR DIS */

//D

channel(#NuBarMode_Dplus_dis_sig_numu)<
        @channel = #DplusNuBarMode:     +:      m:      m:      #CC:    #dis_numu_sig
        @post_smearing_efficiencies = copy(%post_dis_RHC_numu_sig)
>

channel(#NuBarMode_Dplus_dis_sig_numubar)<
        @channel = #DplusNuBarMode_minus:     -:      m:      m:      #CC:    #dis_numubar_sig
        @post_smearing_efficiencies = copy(%post_dis_RHC_numubar_sig)
>

channel(#NuBarMode_Dplus_dis_bkg_nutau)<
        @channel = #DplusNuBarMode:     +:      m:      t:      #CCtau:    #dis_nutau_bkg
        @post_smearing_efficiencies = copy(%post_dis_RHC_nutau_bkg)
>

channel(#NuBarMode_Dplus_dis_bkg_nutaubar)<
        @channel = #DplusNuBarMode_minus:     -:      m:      t:      #CCtau:    #dis_nutaubar_bkg
        @post_smearing_efficiencies = copy(%post_dis_RHC_nutaubar_bkg)
>

channel(#NuBarMode_Dplus_dis_bkg_nuNC)<
        @channel =      #DplusNuBarMode:     +:      m: NOSC_m:      #NC:    #dis_NC_bkg
        @post_smearing_efficiencies = copy(%post_dis_RHC_NC_bkg)
>

channel(#NuBarMode_Dplus_dis_bkg_nubarNC)<
        @channel =      #DplusNuBarMode_minus:     -:      m: NOSC_m:      #NC:    #dis_aNC_bkg
        @post_smearing_efficiencies = copy(%post_dis_RHC_NC_bkg)
>

channel(#NuBarMode_Dplus_dis_bkg_nuNC_e)<
        @channel =      #DplusNuBarMode:     +:      e: NOSC_e:      #NC:    #dis_NC_bkg
        @post_smearing_efficiencies = copy(%post_dis_RHC_NC_bkg)
>

channel(#NuBarMode_Dplus_dis_bkg_nubarNC_e)<
        @channel =      #DplusNuBarMode_minus:     -:      e: NOSC_e:      #NC:    #dis_aNC_bkg
        @post_smearing_efficiencies = copy(%post_dis_RHC_NC_bkg)
>

channel(#NuBarMode_Dplus_dis_bkg_nuNC_tau)<
        @channel =      #DplusNuBarMode:     +:      t: NOSC_t:      #NC:    #dis_NC_bkg
        @post_smearing_efficiencies = copy(%post_dis_RHC_NC_bkg)
>

channel(#NuBarMode_Dplus_dis_bkg_nubarNC_tau)<
        @channel =      #DplusNuBarMode_minus:     -:      t: NOSC_t:      #NC:    #dis_aNC_bkg
        @post_smearing_efficiencies = copy(%post_dis_RHC_NC_bkg)
>



channel(#NuBarMode_Dplus_dis_bkg_e_nutau)<
        @channel = #DplusNuBarMode:     +:      e:      t:      #CCtau:    #dis_nutau_bkg
        @post_smearing_efficiencies = copy(%post_dis_RHC_nutau_bkg)
>

channel(#NuBarMode_Dplus_dis_bkg_e_nutaubar)<
        @channel = #DplusNuBarMode_minus:     -:      e:      t:      #CCtau:    #dis_nutaubar_bkg
        @post_smearing_efficiencies = copy(%post_dis_RHC_nutaubar_bkg)
>

channel(#NuBarMode_Dplus_dis_sig_e_numu)<
        @channel = #DplusNuBarMode:     +:      e:      m:      #CC:    #dis_numu_sig
        @post_smearing_efficiencies = copy(%post_dis_RHC_numu_sig)
>

channel(#NuBarMode_Dplus_dis_sig_e_numubar)<
        @channel = #DplusNuBarMode_minus:     -:      e:      m:      #CC:    #dis_numubar_sig
        @post_smearing_efficiencies = copy(%post_dis_RHC_numubar_sig)
>

channel(#NuBarMode_Dplus_dis_sig_t_numu)<
        @channel = #DplusNuBarMode:     +:      t:      m:      #CC:    #dis_numu_sig
        @post_smearing_efficiencies = copy(%post_dis_RHC_numu_sig)
>

channel(#NuBarMode_Dplus_dis_sig_t_numubar)<
        @channel = #DplusNuBarMode_minus:     -:      t:      m:      #CC:    #dis_numubar_sig
        @post_smearing_efficiencies = copy(%post_dis_RHC_numubar_sig)
>


//Ds

channel(#NuBarMode_DS_dis_sig_numu)<
        @channel =      #DSNuBarMode:     +:      m:      m:      #CC:    #dis_numu_sig
        @post_smearing_efficiencies = copy(%post_dis_RHC_numu_sig)
>

channel(#NuBarMode_DS_dis_sig_numubar)<
        @channel =      #DSNuBarMode_minus:     -:      m:      m:      #CC:    #dis_numubar_sig
        @post_smearing_efficiencies = copy(%post_dis_RHC_numubar_sig)
>

channel(#NuBarMode_DS_dis_bkg_nutau)<
        @channel =      #DSNuBarMode:     +:      m:      t:      #CCtau:    #dis_nutau_bkg
        @post_smearing_efficiencies = copy(%post_dis_RHC_nutau_bkg)
>

channel(#NuBarMode_DS_dis_bkg_nutaubar)<
        @channel =      #DSNuBarMode_minus:     -:      m:      t:      #CCtau:    #dis_nutaubar_bkg
        @post_smearing_efficiencies = copy(%post_dis_RHC_nutaubar_bkg)
>

channel(#NuBarMode_DS_dis_bkg_nuNC)<
        @channel =      #DSNuBarMode:     +:      m: NOSC_m:      #NC:    #dis_NC_bkg
        @post_smearing_efficiencies = copy(%post_dis_RHC_NC_bkg)
>

channel(#NuBarMode_DS_dis_bkg_nubarNC)<
        @channel =      #DSNuBarMode_minus:     -:      m: NOSC_m:      #NC:    #dis_aNC_bkg
        @post_smearing_efficiencies = copy(%post_dis_RHC_NC_bkg)
>

channel(#NuBarMode_DS_dis_bkg_nuNC_e)<
        @channel =      #DSNuBarMode:     +:      e: NOSC_e:      #NC:    #dis_NC_bkg
        @post_smearing_efficiencies = copy(%post_dis_RHC_NC_bkg)
>

channel(#NuBarMode_DS_dis_bkg_nubarNC_e)<
        @channel =      #DSNuBarMode_minus:     -:      e: NOSC_e:      #NC:    #dis_aNC_bkg
        @post_smearing_efficiencies = copy(%post_dis_RHC_NC_bkg)
>

channel(#NuBarMode_DS_dis_bkg_nuNC_tau)<
        @channel =      #DSNuBarMode:     +:      t: NOSC_t:      #NC:    #dis_NC_bkg
        @post_smearing_efficiencies = copy(%post_dis_RHC_NC_bkg)
>

channel(#NuBarMode_DS_dis_bkg_nubarNC_tau)<
        @channel =      #DSNuBarMode_minus:     -:      t: NOSC_t:      #NC:    #dis_aNC_bkg
        @post_smearing_efficiencies = copy(%post_dis_RHC_NC_bkg)
>

channel(#NuBarMode_DS_dis_bkg_e_nutau)<
        @channel = #DSNuBarMode:     +:      e:      t:      #CCtau:    #dis_nutau_bkg
        @post_smearing_efficiencies = copy(%post_dis_RHC_nutau_bkg)
>

channel(#NuBarMode_DS_dis_bkg_e_nutaubar)<
        @channel = #DSNuBarMode_minus:     -:      e:      t:      #CCtau:    #dis_nutaubar_bkg
        @post_smearing_efficiencies = copy(%post_dis_RHC_nutaubar_bkg)
>

channel(#NuBarMode_DS_dis_sig_e_numu)<
        @channel = #DSNuBarMode:     +:      e:      m:      #CC:    #dis_numu_sig
        @post_smearing_efficiencies = copy(%post_dis_RHC_numu_sig)
>

channel(#NuBarMode_DS_dis_sig_e_numubar)<
        @channel = #DSNuBarMode_minus:     -:      e:      m:      #CC:    #dis_numubar_sig
        @post_smearing_efficiencies = copy(%post_dis_RHC_numubar_sig)
>

channel(#NuBarMode_DS_dis_sig_t_numu)<
        @channel = #DSNuBarMode:     +:      t:      m:      #CC:    #dis_numu_sig
        @post_smearing_efficiencies = copy(%post_dis_RHC_numu_sig)
>

channel(#NuBarMode_DS_dis_sig_t_numubar)<
        @channel = #DSNuBarMode_minus:     -:      t:      m:      #CC:    #dis_numubar_sig
        @post_smearing_efficiencies = copy(%post_dis_RHC_numubar_sig)
>


//Kaon

channel(#NuBarMode_Kaon_dis_sig_numu)<
        @channel =      #KaonNuBarMode:     +:      m:      m:      #CC:    #dis_numu_sig
        @post_smearing_efficiencies = copy(%post_dis_RHC_numu_sig)
>

channel(#NuBarMode_Kaon_dis_sig_numubar)<
        @channel =      #KaonNuBarMode_minus:     -:      m:      m:      #CC:    #dis_numubar_sig
        @post_smearing_efficiencies = copy(%post_dis_RHC_numubar_sig)
>

channel(#NuBarMode_Kaon_dis_bkg_nutau)<
        @channel =      #KaonNuBarMode:     +:      m:      t:      #CCtau:    #dis_nutau_bkg
        @post_smearing_efficiencies = copy(%post_dis_RHC_nutau_bkg)
>

channel(#NuBarMode_Kaon_dis_bkg_nutaubar)<
        @channel =      #KaonNuBarMode:      -:      m:      t:      #CCtau:    #dis_nutaubar_bkg
        @post_smearing_efficiencies = copy(%post_dis_RHC_nutaubar_bkg)
>

channel(#NuBarMode_Kaon_dis_bkg_nuNC)<
        @channel =      #KaonNuBarMode:     +:      m: NOSC_m:      #NC:    #dis_NC_bkg
        @post_smearing_efficiencies = copy(%post_dis_RHC_NC_bkg)
>

channel(#NuBarMode_Kaon_dis_bkg_nubarNC)<
        @channel =      #KaonNuBarMode_minus:     -:      m: NOSC_m:      #NC:    #dis_aNC_bkg
        @post_smearing_efficiencies = copy(%post_dis_RHC_NC_bkg)
>

channel(#NuBarMode_Kaon_dis_bkg_nuNC_e)<
        @channel =      #KaonNuBarMode:     +:      e: NOSC_e:      #NC:    #dis_NC_bkg
        @post_smearing_efficiencies = copy(%post_dis_RHC_NC_bkg)
>

channel(#NuBarMode_Kaon_dis_bkg_nubarNC_e)<
        @channel =      #KaonNuBarMode_minus:     -:      e: NOSC_e:      #NC:    #dis_aNC_bkg
        @post_smearing_efficiencies = copy(%post_dis_RHC_NC_bkg)
>

channel(#NuBarMode_Kaon_dis_bkg_nuNC_tau)<
        @channel =      #KaonNuBarMode:     +:      t: NOSC_t:      #NC:    #dis_NC_bkg
        @post_smearing_efficiencies = copy(%post_dis_RHC_NC_bkg)
>

channel(#NuBarMode_Kaon_dis_bkg_nubarNC_tau)<
        @channel =      #KaonNuBarMode_minus:     -:      t: NOSC_t:      #NC:    #dis_aNC_bkg
        @post_smearing_efficiencies = copy(%post_dis_RHC_NC_bkg)
>


channel(#NuBarMode_Kaon_dis_bkg_e_nutau)<
        @channel = #KaonNuBarMode:     +:      e:      t:      #CCtau:    #dis_nutau_bkg
        @post_smearing_efficiencies = copy(%post_dis_RHC_nutau_bkg)
>

channel(#NuBarMode_Kaon_dis_bkg_e_nutaubar)<
        @channel = #KaonNuBarMode_minus:     -:      e:      t:      #CCtau:    #dis_nutaubar_bkg
        @post_smearing_efficiencies = copy(%post_dis_RHC_nutaubar_bkg)
>

channel(#NuBarMode_Kaon_dis_sig_e_numu)<
        @channel = #KaonNuBarMode:     +:      e:      m:      #CC:    #dis_numu_sig
        @post_smearing_efficiencies = copy(%post_dis_RHC_numu_sig)
>

channel(#NuBarMode_Kaon_dis_sig_e_numubar)<
        @channel = #KaonNuBarMode_minus:     -:      e:      m:      #CC:    #dis_numubar_sig
        @post_smearing_efficiencies = copy(%post_dis_RHC_numubar_sig)
>

channel(#NuBarMode_Kaon_dis_sig_t_numu)<
        @channel = #KaonNuBarMode:     +:      t:      m:      #CC:    #dis_numu_sig
        @post_smearing_efficiencies = copy(%post_dis_RHC_numu_sig)
>

channel(#NuBarMode_Kaon_dis_sig_t_numubar)<
        @channel = #KaonNuBarMode_minus:     -:      t:      m:      #CC:    #dis_numubar_sig
        @post_smearing_efficiencies = copy(%post_dis_RHC_numubar_sig)
>


//Pion

channel(#NuBarMode_Pion_dis_sig_numu)<
        @channel =      #PionNuBarMode:     +:      m:      m:      #CC:    #dis_numu_sig
        @post_smearing_efficiencies = copy(%post_dis_RHC_numu_sig)
>

channel(#NuBarMode_Pion_dis_sig_numubar)<
        @channel =      #PionNuBarMode_minus:     -:      m:      m:      #CC:    #dis_numubar_sig
        @post_smearing_efficiencies = copy(%post_dis_RHC_numubar_sig)
>

channel(#NuBarMode_Pion_dis_bkg_nutau)<
        @channel =      #PionNuBarMode:     +:      m:      t:      #CCtau:    #dis_nutau_bkg
        @post_smearing_efficiencies = copy(%post_dis_RHC_nutau_bkg)
>

channel(#NuBarMode_Pion_dis_bkg_nutaubar)<
        @channel =      #PionNuBarMode_minus:     -:      m:      t:      #CCtau:   #dis_nutaubar_bkg
        @post_smearing_efficiencies = copy(%post_dis_RHC_nutaubar_bkg)
>

channel(#NuBarMode_Pion_dis_bkg_nuNC)<
        @channel =      #PionNuBarMode:     +:      m: NOSC_m:      #NC:    #dis_NC_bkg
        @post_smearing_efficiencies = copy(%post_dis_RHC_NC_bkg)
>

channel(#NuBarMode_Pion_dis_bkg_nubarNC)<
        @channel =      #PionNuBarMode_minus:     -:      m: NOSC_m:      #NC:    #dis_aNC_bkg
        @post_smearing_efficiencies = copy(%post_dis_RHC_NC_bkg)
>

channel(#NuBarMode_Pion_dis_bkg_nuNC_e)<
        @channel =      #PionNuBarMode:     +:      e: NOSC_e:      #NC:    #dis_NC_bkg
        @post_smearing_efficiencies = copy(%post_dis_RHC_NC_bkg)
>

channel(#NuBarMode_Pion_dis_bkg_nubarNC_e)<
        @channel =      #PionNuBarMode_minus:     -:      e: NOSC_e:      #NC:    #dis_aNC_bkg
        @post_smearing_efficiencies = copy(%post_dis_RHC_NC_bkg)
>

channel(#NuBarMode_Pion_dis_bkg_nuNC_tau)<
        @channel =      #PionNuBarMode:     +:      t: NOSC_t:      #NC:    #dis_NC_bkg
        @post_smearing_efficiencies = copy(%post_dis_RHC_NC_bkg)
>

channel(#NuBarMode_Pion_dis_bkg_nubarNC_tau)<
        @channel =      #PionNuBarMode_minus:     -:      t: NOSC_t:      #NC:    #dis_aNC_bkg
        @post_smearing_efficiencies = copy(%post_dis_RHC_NC_bkg)
>

channel(#NuBarMode_Pion_dis_bkg_e_nutau)<
        @channel = #PionNuBarMode:     +:      e:      t:      #CCtau:    #dis_nutau_bkg
        @post_smearing_efficiencies = copy(%post_dis_RHC_nutau_bkg)
>

channel(#NuBarMode_Pion_dis_bkg_e_nutaubar)<
        @channel = #PionNuBarMode_minus:     -:      e:      t:      #CCtau:    #dis_nutaubar_bkg
        @post_smearing_efficiencies = copy(%post_dis_RHC_nutaubar_bkg)
>

channel(#NuBarMode_Pion_dis_sig_e_numu)<
        @channel = #PionNuBarMode:     +:      e:      m:      #CC:    #dis_numu_sig
        @post_smearing_efficiencies = copy(%post_dis_RHC_numu_sig)
>

channel(#NuBarMode_Pion_dis_sig_e_numubar)<
        @channel = #PionNuBarMode_minus:     -:      e:      m:      #CC:    #dis_numubar_sig
        @post_smearing_efficiencies = copy(%post_dis_RHC_numubar_sig)
>

channel(#NuBarMode_Pion_dis_sig_t_numu)<
        @channel = #PionNuBarMode:     +:      t:      m:      #CC:    #dis_numu_sig
        @post_smearing_efficiencies = copy(%post_dis_RHC_numu_sig)
>

channel(#NuBarMode_Pion_dis_sig_t_numubar)<
        @channel = #PionNuBarMode_minus:     -:      t:      m:      #CC:    #dis_numubar_sig
        @post_smearing_efficiencies = copy(%post_dis_RHC_numubar_sig)
>

//Muon 

channel(#NuBarMode_Muon_dis_sig_numu)<
        @channel =      #MuonNuBarMode:     +:      m:      m:      #CC:    #dis_numu_sig
        @post_smearing_efficiencies = copy(%post_dis_RHC_numu_sig)
>

channel(#NuBarMode_Muon_dis_sig_numubar)<
        @channel =      #MuonNuBarMode_minus:     -:      m:      m:      #CC:    #dis_numubar_sig
        @post_smearing_efficiencies = copy(%post_dis_RHC_numubar_sig)
>

channel(#NuBarMode_Muon_dis_bkg_nutau)<
        @channel =      #MuonNuBarMode:     +:      m:      t:      #CCtau:    #dis_nutau_bkg
        @post_smearing_efficiencies = copy(%post_dis_RHC_nutau_bkg)
>

channel(#NuBarMode_Muon_dis_bkg_nutaubar)<
        @channel =      #MuonNuBarMode_minus:     -:      m:      t:      #CCtau:   #dis_nutaubar_bkg
        @post_smearing_efficiencies = copy(%post_dis_RHC_nutaubar_bkg)
>

channel(#NuBarMode_Muon_dis_bkg_nuNC)<
        @channel =      #MuonNuBarMode:     +:      m: NOSC_m:      #NC:    #dis_NC_bkg
        @post_smearing_efficiencies = copy(%post_dis_RHC_NC_bkg)
>

channel(#NuBarMode_Muon_dis_bkg_nubarNC)<
        @channel =      #MuonNuBarMode_minus:     -:      m: NOSC_m:      #NC:    #dis_aNC_bkg
        @post_smearing_efficiencies = copy(%post_dis_RHC_NC_bkg)
>

channel(#NuBarMode_Muon_dis_bkg_nuNC_e)<
        @channel =      #MuonNuBarMode:     +:      e: NOSC_e:      #NC:    #dis_NC_bkg
        @post_smearing_efficiencies = copy(%post_dis_RHC_NC_bkg)
>

channel(#NuBarMode_Muon_dis_bkg_nubarNC_e)<
        @channel =      #MuonNuBarMode_minus:     -:      e: NOSC_e:      #NC:    #dis_aNC_bkg
        @post_smearing_efficiencies = copy(%post_dis_RHC_NC_bkg)
>

channel(#NuBarMode_Muon_dis_bkg_nuNC_tau)<
        @channel =      #MuonNuBarMode:     +:      t: NOSC_t:      #NC:    #dis_NC_bkg
        @post_smearing_efficiencies = copy(%post_dis_RHC_NC_bkg)
>

channel(#NuBarMode_Muon_dis_bkg_nubarNC_tau)<
        @channel =      #MuonNuBarMode_minus:     -:      t: NOSC_t:      #NC:    #dis_aNC_bkg
        @post_smearing_efficiencies = copy(%post_dis_RHC_NC_bkg)
>



channel(#NuBarMode_Muon_dis_bkg_e_nutau)<
        @channel = #MuonNuBarMode:     +:      e:      t:      #CCtau:    #dis_nutau_bkg
        @post_smearing_efficiencies = copy(%post_dis_RHC_nutau_bkg)
>

channel(#NuBarMode_Muon_dis_bkg_e_nutaubar)<
        @channel = #MuonNuBarMode_minus:     -:      e:      t:      #CCtau:    #dis_nutaubar_bkg
        @post_smearing_efficiencies = copy(%post_dis_RHC_nutaubar_bkg)
>

channel(#NuBarMode_Muon_dis_sig_e_numu)<
        @channel = #MuonNuBarMode:     +:      e:      m:      #CC:    #dis_numu_sig
        @post_smearing_efficiencies = copy(%post_dis_RHC_numu_sig)
>

channel(#NuBarMode_Muon_dis_sig_e_numubar)<
        @channel = #MuonNuBarMode_minus:     -:      e:      m:      #CC:    #dis_numubar_sig
        @post_smearing_efficiencies = copy(%post_dis_RHC_numubar_sig)
>

channel(#NuBarMode_Muon_dis_sig_t_numu)<
        @channel = #MuonNuBarMode:     +:      t:      m:      #CC:    #dis_numu_sig
        @post_smearing_efficiencies = copy(%post_dis_RHC_numu_sig)
>

channel(#NuBarMode_Muon_dis_sig_t_numubar)<
        @channel = #MuonNuBarMode_minus:     -:      t:      m:      #CC:    #dis_numubar_sig
        @post_smearing_efficiencies = copy(%post_dis_RHC_numubar_sig)
>


