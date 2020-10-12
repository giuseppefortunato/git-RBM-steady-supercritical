//
// Created by giuseppe on 01/10/2020.
//
#include "Post-Process.hpp"



void Write_History_ResError_global(prob_settings settings){

    //Defining Common Variables
    std::cout<<"Write global file of residual errors"<<std::endl;
    std::string header_history, value_history; //to write headers in global files
    std::ifstream history_su2;
    std::ofstream history_global;
    //create the name for the history global file

    std::string filename_history_global={};
    if(settings.flag_method[0]== "ISOMAP")
        filename_history_global= "history_global_rbm_"+settings.flag_method[0]+"_"+std::to_string(settings.r_isomap)+".csv";
    else if (settings.flag_method[0]=="POD")
        filename_history_global= "history_global_rbm_"+settings.flag_method[0]+"_"+std::to_string(settings.r)+".csv";


    //Getting filename
    std::string filename_history_su2 ={};
    history_global.open(filename_history_global);
    for (int nt=0; nt < settings.param_rec_1.size(); nt++) {
       if (settings.flag_method[0] == "ISOMAP") {
           filename_history_su2 = "history_" + settings.flag_method[0] + "_" + std::to_string(settings.r_isomap) + "_" +
                                  std::to_string(nt + 1) + ".csv";
       }else if (settings.flag_method[0] == "POD"){
           filename_history_su2 = "history_" + settings.flag_method[0] + "_" + std::to_string(settings.r) + "_" +
                                  std::to_string(nt + 1) + ".csv";
       }


        //se nt=0, apri il file globale e scrivi l'header
        history_su2.open(filename_history_su2);
        if(history_su2.is_open()){
            std::string linedata1, linedata2;
            //Getting row of headers
            getline(history_su2,linedata1);
            header_history = linedata1;
            //Getting values
            getline(history_su2,linedata2);
            value_history = linedata2;
        } else {
            std::cout << "Unable to open SU2 single history file. Exiting ..." << std::endl;
            exit(EXIT_FAILURE);
        }
        history_su2.close();
        //history_global.open(filename_history_global);
        if(history_global.is_open()){
            if( nt==0 ) {
                history_global << header_history << std::endl;
                history_global << std::to_string(nt+1) << "_" << std::to_string(settings.param_rec_1[nt]) << "_" << std::to_string(settings.param_rec_2[nt])
                               << "," << value_history << std::endl;
            }else{
                history_global << std::to_string(nt + 1) << "_" << std::to_string(settings.param_rec_1[nt]) << "_"
                               << std::to_string(settings.param_rec_2[nt])
                               << "," << value_history << std::endl;
            }
        }else{
            std::cout << "Unable to open global file for writing. Exiting... " << std::endl;
            exit(EXIT_FAILURE);
        }

        //String for removing useless solutions (single history files that now are included in the global one)
        std::string rmf_string = "rm -f " + filename_history_su2;
        int len_s = rmf_string.length();
        char rmf_sys_call[len_s + 20];
        strcpy(rmf_sys_call, rmf_string.c_str());
        std::system(rmf_sys_call);

        //string for removing the output.log from SU2_CFD / SU2_SOL---->if I want to see these files, comment these lines.
        //std::string rmf_string_su2log= "rm -f SU2_"+settings.flag_method[0]+"_"+std::to_string(nt+1)+".log";
        //len_s= rmf_string_su2log.length();
        //rmf_sys_call[len_s + 20];
        //strcpy(rmf_sys_call, rmf_string_su2log.c_str());
        //std::system(rmf_sys_call);

    }


    history_global.close();
    std::cout<<"DONE"<<std::endl;
}



