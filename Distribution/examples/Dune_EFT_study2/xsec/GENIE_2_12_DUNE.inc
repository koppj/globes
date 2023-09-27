
$version="3.2.13"

/* Cross sections */
cross(#CC)<
        @cross_file = "./xsec/xsec_cc.dat"
        @eft_coeff_file = "./xsec/detection_coefficients/NeutrinosCCQEOxygen.txt"
        @quark_flavors = 0 : 0
>

cross(#CCtau)<
        @cross_file = "./xsec/xsec_tau.dat" //They have set the tau cross section to 0 in the other file
        @eft_coeff_file = "./xsec/detection_coefficients/NeutrinosCCQEOxygen.txt"
        @quark_flavors = 0 : 0
>

cross(#NC)<
        @cross_file = "./xsec/xsec_nc.dat"
>


