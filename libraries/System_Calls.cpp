#include "System_Calls.hpp"

#include <utility>


void SU2_DTR(prob_settings settings,std::string filename, const std::string& su2_conf, int nt ){


    /*------------Defining all necessary strings----------------*/

    //New cfg file SU2 for residual evaluation and modification for SU2_DTR
    std::string root_conf;
    root_conf.assign ( su2_conf, 0, su2_conf.size() - 4);
    std::string su2_conf_new = root_conf + "-reseval.cfg";


    Modify_su2_cfg( su2_conf, su2_conf_new, std::move(filename), settings, nt);

    //String to launch SU2_DTR from terminal
    std::string su2dtr_string = "./SU2_CFD " + su2_conf_new + " > SU2_" + settings.flag_method[0] + "_" + std::to_string (nt+1) +".log";
    std::string su2sol_string = "./SU2_SOL " + su2_conf_new + "> SU2_SOL" + settings.flag_method[0] + "_" +std::to_string(nt+1)+ ".log";
    int len_s = su2dtr_string.length();
    int len_ss= su2sol_string.length();
    char su2_sys_call[len_s + 1];
    char su2sol_sys_call[len_ss +1];
    strcpy(su2_sys_call, su2dtr_string.c_str());
    strcpy(su2sol_sys_call,su2sol_string.c_str());



    if (settings.flag_rec== "YES"){
        std::cout << "Calling SU2_CFD for residual evaluation and writing file to history " << std::endl;
        std::system(su2_sys_call);
    }else if (settings.flag_rec=="NO"){
        std::cout <<" Calling SU2_SOL for visualization" <<std::endl;
        std::system(su2sol_sys_call);
    }

    std::cout << std::endl;
}

