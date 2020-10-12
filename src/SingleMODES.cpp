#include "Extract_Basis.hpp"
#include "read_Inputs.hpp"
#include "Generate_snset.hpp"
#include "System_Calls.hpp"
#include "Post-Process.hpp"
#include "Reconstruction.hpp"
#include "write_Outputs.hpp"
#include "manifold_learning.hpp"
#include <time.h>

int main(int argc, char *argv[]) {

    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << "-----------Single-MODES start-------------" << std::endl << std::endl;

    prob_settings settings;
    std::string filecfg = argv[1];
    //std::string su2_conf = argv[2];     //for gems i don't need this parameter because i don't want call SU2


    //Reading configuration file
    Read_cfg( filecfg, settings );
    Config_stream ( settings );

    // Calculate number of grid points
    std::string root_inputfile;
    root_inputfile.assign ( settings.in_file, 0, settings.in_file.size() - 4);
    std::string input_format;
    input_format.assign ( settings.in_file, settings.in_file.size() - 3, 3);

    std::stringstream buffer;
    buffer << std::setfill('0') << std::setw(5) << std::to_string(settings.nstart);
    std::string file_1 = root_inputfile + "_" + buffer.str() + "." + input_format;
    //provvisorio for supercritical
    //std::stringstream buffer_2;
    //buffer_2 << std::setfill('0') << std::setw(5) << std::to_string(1);
    //std::string file_cord = root_inputfile+ "_"+buffer_2.str()+"."+input_format;

    //DEFINE PARAMETERS

    Eigen::MatrixXd param(settings.Ns, settings.Np);
    Eigen::RowVectorXd norme(settings.Np);
    std::string file_parameters = argv[2];                          //if gems this the second argument passed, otherwise it is the third
    define_parameters(settings, param, file_parameters);

    Eigen::RowVectorXd supercritical(settings.Np);
    for (int i=0; i<settings.Np; i++)
        supercritical(i)= settings.supercritical_rec[i];



    //COMPUTING PARAMETERS NORM
    for (int i=0; i<settings.Np; i++)
        norme(i)= param.col(i).norm();


    int Nr = N_gridpoints ( file_1 );
    std::cout << "Number of grid points : " << Nr << std::endl;
    int nC = settings.Cols.size();
    int s_Nf = 1;   //Number of values for the SPOD filter (POD included)
    std::vector<int> Nf(s_Nf);
    Nf[0] = 0;
    int nfj = 0;
    
    // Reading coordinates
    std::cout << "Reading Coordinates ... \t ";
    Eigen::MatrixXd Coords = read_col( file_1, Nr, settings.Cols_coords );          //file_1 in origine, file_cord se supercritical
    std::cout << "Done " << std::endl;

    // Create matrix of snapshots
    clock_t start,end;
    double tempo;
    start = clock();
    std::cout << "Storing snapshot Matrix ... \n ";
    Eigen::MatrixXd sn_set = generate_snap_matrix( Nr, settings.Ns, settings.Ds, settings.nstart,
                                        settings.Cols,
                                        settings.in_file,
                                        settings.flag_prob);

    Eigen::VectorXd mean = sn_set.rowwise().mean();
//    Eigen::VectorXd Ic = IC(settings, ,Nr)
    std::cout<<"done"<<std::endl;
    end = clock();
    tempo = ((double)(end-start))/CLOCKS_PER_SEC;
    std::cout<<"tempo per leggere la matrice degli snnapshot:"<< std::setprecision(10)<<tempo<<std::endl;

    if ( settings.flag_mean == "YES" ) {
        std::cout << "Subtracting mean from snapshots ... " << std::endl << std::endl;
        for ( int i = 0; i < settings.Ns; i++ )
            sn_set.col(i) -= mean;
    } else if ( settings.flag_mean == "IC") {
        std::cout << "Mean is the only reference solution implemented in SingleMODES\n Exiting ... " << std::endl << std::endl;
        exit(EXIT_FAILURE);
    }

    Eigen::VectorXd K_pc(settings.Ns);

   /* if (settings.flag_method[0] == "POD"){
        if(settings.analysis == "STEADY") {
            std::cout << "-------------------" << std::endl;
            std::cout << "POD-STEADY ANALYSIS" << std::endl;
            std::cout << "-------------------" << std::endl;
        }else if( settings.analysis == "UNSTEADY"){
            std::cout << "---------------------" << std::endl;
            std::cout << "POD-UNSTEADY ANALYSIS" << std::endl;
            std::cout << "---------------------" << std::endl;
        }

        Eigen::MatrixXd Sn_Cons_time = Eigen::MatrixXd::Zero(nC*Nr, 3);
        std::vector<Eigen::MatrixXd> Phi(nC);
        std::vector<Eigen::VectorXd> lambda(nC);
        std::vector<Eigen::MatrixXd> coefs = {};

        //std::vector< std::vector<rbf> > surr_coefs(nC);
        Eigen::MatrixXd eig_vec(settings.Ns, settings.Ns);
        std::vector<int> Nm(nC);
        int N_notZero;
        //Check only for POD for now
        for (int i = 0; i < nC; i++) {
            Phi[i] = Eigen::MatrixXd::Zero(Nr,settings.Ns);
            lambda[i] = Eigen::VectorXd::Zero(settings.Ns);
        }

        std::cout << "Extracting basis ... " << std::endl << std::endl;
        for (int ncons=0; ncons<nC; ncons++){
           std::cout << "Processing conservative variable " << ncons << std::endl;
           Phi[ncons] = SPOD_basis( sn_set.middleRows(ncons*Nr,Nr),
                                          lambda[ncons], K_pc, eig_vec,
                                          Nf[nfj],
                                          settings.flag_bc,
                                          settings.flag_filter,
                                          settings.sigma);

           std::cout << " Done! " << std::endl << std::endl;
           N_notZero = Phi[ncons].cols();
           if (settings.r == 0) Nm[ncons] = Nmod(settings.En, K_pc);
           else Nm[ncons] = std::min(settings.r, N_notZero);
           std::cout << "Number of modes used in reconstruction " << Nm[ncons] << std::endl;
           coefs.push_back(eig_vec.transpose());
           //surr_coefs[ncons] = getSurrCoefs(t_vec, eig_vec, settings.flag_interp);

        }

        Eigen::MatrixXd Field(Nr, nC);


        for (int nt = 0; nt < settings.param_rec_1.size(); nt++) {                  // nt < settings.param_rec_1.size() original
            std::cout << "Reconstructing field for:" << settings.param_rec_1[nt]
                          << "\t" << settings.param_rec_2[nt]
                          << "\t" << settings.param_rec_3[nt]
                          << "\t" << settings.param_rec_4[nt]
                          << "\t" << settings.param_rec_5[nt]
                          << "\t" << settings.param_rec_6[nt] << std::endl;
            for (int ncons = 0; ncons < nC; ncons++) {
                    Eigen::MatrixXd Rec = Reconstruction_POD(param, norme, K_pc, lambda[ncons], coefs[ncons],
                                                             Phi[ncons],
                                                             settings.param_rec_1[nt],
                                                             settings.param_rec_2[nt],
                                                             settings.param_rec_3[nt],
                                                             settings.param_rec_4[nt],
                                                             settings.param_rec_5[nt],
                                                             settings.param_rec_6[nt],
                                                             Nm[ncons],
                                                             settings.Np,
                                                             settings.flag_interp,
                                                             supercritical);
                    Field.col(ncons) = Rec;
                    std::cout << "Done" << std::endl;
            }

            std::cout<<"creating reconstructed fields file"<<std::endl;
            std::string filename;
            filename = write_reconstruction_file(Field, Coords, nt, nC, Nm[0], settings);
            SU2_DTR(settings, filename, su2_conf, nt);
            std::cout <<"DONE"<<std::endl;

        }

        for(int nt =0; nt<settings.param_rec_1.size(); nt++){


        }

        Write_History_ResError_global(settings);

        std::cout<<"------------------END OF POD--------------------" <<std::endl;
    }*/

   if (settings.flag_method[0] == "ISOMAP"){
        if( settings.analysis== "STEADY"){
        std::cout<<"----------------------"<<std::endl;
        std::cout<<"--ISOMAP STEADY + I --"<<std::endl;
        std::cout<<"----------------------"<<std::endl;
        }else if (settings.analysis== "UNSTEADY"){
        std::cout<<"------------------------"<<std::endl;
        std::cout<<"--ISOMAP UNSTEADY + I --"<<std::endl;
        std::cout<<"------------------------"<<std::endl;
        }


        Eigen::MatrixXd snapshot_star_total = Eigen::MatrixXd::Zero(Nr, nC);
        Eigen::MatrixXd sn_set_partial(Nr, sn_set.cols());
        Eigen::VectorXd norma(sn_set.cols());
        Eigen::MatrixXd y(settings.r_isomap*nC, sn_set.cols());
        Eigen::VectorXi k(nC);
        Eigen::VectorXd kruskal_min = Eigen::VectorXd::Ones(nC);


        //Loop for each conservative variable---> ISOMAP
        for(int q=0; q< nC; q++){
            std::cout << "ISOMAP FOR CONSERVATIVE VARIABLE NUMBER " << "\t" << q << std::endl;
            ISOMAP(sn_set, settings, Nr, q, K_pc, k, kruskal_min, y);
        }

        std::cout<<"DONE "<<std::endl;
        std::cout << "Kruskal stress minimo \t " << kruskal_min << "\t per k=" << k << std::endl;
        std::cout << y << std::endl<<std::endl<<std::endl;
        std::cout<<"END OF ISOMAP..."<<std::endl<<std::endl;

        std::cout<<"------------START BACK-MAPPING SECTION--------------"<<std::endl<<std::endl;

        for(int nt=0; nt < 1; nt++){                  //nt < settings.param_rec_1.size()...1 nel caso di supercritical
            std::cout << "RECONSTRUCTING FIELD FOR:\t" << settings.param_rec_1[nt]
                      << "\t" << settings.param_rec_2[nt]
                      << "\t" << settings.param_rec_3[nt]
                      << "\t" << settings.param_rec_4[nt]
                      << "\t" << settings.param_rec_5[nt]
                      << "\t" << settings.param_rec_6[nt] << std::endl;

            for (int q=0; q<nC; q++) {

                sn_set_partial = sn_set.block(Nr * q, 0, Nr, sn_set.cols());
                for (int i = 0; i < sn_set.cols(); i++) {
                    norma(i) = sn_set_partial.col(i).norm();
                    sn_set_partial.col(i) = sn_set_partial.col(i) / norma(i);
                }
                Eigen::MatrixXd y_part(settings.r_isomap, sn_set.cols());
                for (int j = 0; j < settings.r_isomap; j++) {
                    y_part.row(j) = y.row((q * settings.r_isomap) + j);
                }
                //std::cout<<y_part<<std::endl;

                std::vector<double> y_star = ISOMAP_surr_coeff(param, norme, y_part,
                                                               settings.param_rec_1[nt],
                                                               settings.param_rec_2[nt],
                                                               settings.param_rec_3[nt],
                                                               settings.param_rec_4[nt],
                                                               settings.param_rec_5[nt],
                                                               settings.param_rec_6[nt],
                                                               settings.Np,
                                                               settings.flag_interp,
                                                               supercritical);

                std::cout << "values RBF" << std::endl;
                for(int i = 0; i < y_star.size(); i++)
                    std::cout << y_star[i] << std::endl;

                Eigen::MatrixXd snapshot_star = create_rec_fields(settings,y_star,y_part,sn_set_partial,
                                                                 norma);

                snapshot_star_total.col(q) = snapshot_star;
            }

            std::cout<<"snap*"<<std::endl<<snapshot_star_total.topRows(5)<<std::endl;
            std::cout<<"creating reconstructed fields file"<<std::endl;
            std::string filename;
            filename = write_reconstruction_file(snapshot_star_total, Coords, nt, nC, settings.r_isomap, settings);
            //SU2_DTR(settings, filename, su2_conf, nt);
            //std::cout << "DONE" << std::endl;

        }

        //Write_History_ResError_global(settings);


        std::cout <<"-----END OF ISOMAP+BACK-MAPPING------"<<std::endl<< std::endl;
    }

    std::cout << std::endl;
    std::cout<<"------------------------------------------" <<std::endl;
    std::cout<<"------------Single-MODES end--------------" <<std::endl;
    std::cout<<"------------------------------------------" <<std::endl;

   return 0;

}
