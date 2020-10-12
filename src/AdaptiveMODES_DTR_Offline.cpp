/*
CODE FOR RESIDUAL EVALUATION OF DIFFERENT RBM TECHNIQUES USING SU2 CODE
INPUT ARGUMENTS
Config File RBM + Config File SU2
*/

#include "Extract_Basis.hpp"
#include "read_Inputs.hpp"
#include "Pre-Process.hpp"
#include "System_Calls.hpp"
#include "Post-Process.hpp"
#include "Generate_snset.hpp"
#include "Reconstruction.hpp"
#include "write_Outputs.hpp"
#include "manifold_learning.hpp"


int main( int argc, char *argv[] )
{

    std::cout << "-----------MODES Adaptive Residual Evaluation starts-------------" << std::endl << std::endl;

    std::cout << "Initializing common variables ... " << std::endl << std::endl;

    prob_settings settings;
    std::string filecfg = argv[1];
    std::string su2_conf = argv[2];
    Read_cfg( filecfg, settings );
    std::string root_outputfile;
    root_outputfile.assign ( settings.out_file, 0, settings.out_file.size() - 4);
    std::string root_inputfile;
    root_inputfile.assign ( settings.in_file, 0, settings.in_file.size() - 4);
    std::string input_format;
    input_format.assign ( settings.in_file, settings.in_file.size() - 3, 3);
    bool direct_error = true;


    Eigen::MatrixXd param(settings.Ns, settings.Np);
    Eigen::RowVectorXd norme(settings.Np);
    std::string file_parameters =argv[3];
    define_parameters(settings, param, file_parameters);

    //COMPUTING PARAMETERS NORM
    for (int i=0; i<settings.Np; i++)
        norme(i)= param.col(i).norm();

    //CONSTRAIN: DESIGN SPACE MUST HAVE THE SAME POINT OF EVALUATION IN EACH AXIS
    Eigen::MatrixXd axis = Eigen::MatrixXd::Zero(50,6); //6 perchè il numero di parametri max che posso leggere da modes
    for (int i=0; i < settings.Np; i++) {
        std::cout <<" define min e max of your design space for the parameter #"<< i+1<< std::endl;
        double low;
        double high;
        std::cin >> low;
        std::cin >> high;
        axis.col(i) = Eigen::VectorXd::LinSpaced(50, low, high);
    }

    std::cout<< axis.topRows(10)<<std::endl;


    //hard coded for now.
    Eigen::MatrixXd design_space_test= Eigen::MatrixXd::Zero(axis.rows()^settings.Np, 6 );
    int interval = 0;
        for (int i = 0; i < axis.rows(); i++) {
            for (int j = 0; j < axis.rows(); j++) {
                design_space_test(interval+j, 0) = axis(i,0);
                design_space_test(interval+j,1) = axis(j,1);
            }
            interval = interval+axis.rows();
        }
    std::cout << design_space_test.topRows(10)<<std::endl;


    /*

    int s_Nf = 1;   //Number of values for the SPOD filter (POD included)
    std::vector<int> Nf(s_Nf);
    Nf[0] = 0;
    int nfj = 0;

    int nC = settings.Cols.size();
    std::vector<double> Dt_res = settings.Dt_res;

    std::stringstream buffer;
    buffer << std::setfill('0') << std::setw(5) << std::to_string(settings.nstart);
    std::string file_1 = root_inputfile + "_" + buffer.str() + "." + input_format;

    // Calculate number of grid points
    int Nr = N_gridpoints ( file_1 );
    std::cout << "Number of grid points : " << Nr << std::endl;

    std::cout << "Reading Coordinates ... \t ";
    Eigen::MatrixXd Coords = read_col( file_1, Nr, settings.Cols_coords );
    std::cout << "Done " << std::endl;

// How we upload the snapshot matrix in the most efficient way?
// By now one big igen Matrix where each column has all the vectors of the conservative variables
    std::cout << "Storing snapshot Matrix ... \n ";
    Eigen::MatrixXd sn_set = generate_snap_matrix( Nr, settings.Ns, settings.Ds, settings.nstart,
                                                   settings.Cols,
                                                   settings.in_file,
                                                   settings.flag_prob);

    std::cout << std::endl;

    std::cout << "Computing mean/Initial Condition of CFD solution ... " << std::endl;
    //Defining Initial condition
    Eigen::VectorXd mean = sn_set.rowwise().mean();
    Eigen::VectorXd Ic = Eigen::VectorXd::Zero(nC*Nr);

    std::string binary = "NO";
    Eigen::VectorXd svd_cum_sum(settings.Ns);

    if ( settings.flag_mean == "IC" ) {
        std::cout << "Subtracting Initial condition" << std::endl << std::endl;
        Ic = IC(sn_set,settings,nC,Nr);
    } else {
        std::cout << "Using data without subtracting any reference state" << std::endl << std::endl;
    }


    auto methods = settings.flag_method;
    //Defining common scope for POD-SPOD
    std::vector<std::string>::iterator itPOD;
    itPOD = std::find (methods.begin(), methods.end(), "POD");
    if (itPOD != methods.end()) {
        std::cout << "--------------------------------------" << std::endl ;
        std::cout << "-------Performin POD ResEval----------" << std::endl ;
        std::cout << "--------------------------------------" << std::endl ;
        //Vector of MatrixXd where to store the evolution in time of conservative variables
        Eigen::MatrixXd Sn_Cons_time = Eigen::MatrixXd::Zero(nC*Nr, 3);
        std::vector<Eigen::MatrixXd> Phi(nC);
        std::vector<Eigen::VectorXd> lambda(nC);

        std::vector< std::vector<rbf> > surr_coefs(nC);
        Eigen::VectorXd K_pc(settings.Ns);
        Eigen::MatrixXd eig_vec(settings.Ns, settings.Ns);
        std::vector<int> Nm(nC);
        int N_notZero;
        //Check only for POD for now
        for (int i = 0; i < nC; i++) {
            Phi[i] = Eigen::MatrixXd::Zero(Nr,settings.Ns);
            lambda[i] = Eigen::VectorXd::Zero(settings.Ns);
        }
        std::cout << std::endl << "Extraction of the basis" << std::endl << std::endl;
        for ( int ncons = 0; ncons < nC; ncons ++ ) {
            std::cout << "Processing conservative variable " << ncons << std::endl;
            Phi[ncons] = SPOD_basis(sn_set.middleRows(ncons * Nr, Nr),
                                    lambda[ncons], K_pc, eig_vec,
                                    Nf[nfj],
                                    settings.flag_bc,
                                    settings.flag_filter,
                                    settings.sigma);
            N_notZero = Phi[ncons].cols();
            if (settings.r == 0) Nm[ncons] = Nmod(settings.En, K_pc);
            else Nm[ncons] = std::min(settings.r, N_notZero);
            std::cout << "Number of modes used in reconstruction " << Nm[ncons] << std::endl;
            surr_coefs[ncons] = getSurrCoefs(t_vec, eig_vec, settings.flag_interp);
        }
        std::cout << std::endl;

//        std::cout << "Computing SPOD " << Nf[nfj] << " reconstruction for each conservative variable ... " << "\n";

        for ( int idtr = 0; idtr < settings.Dt_res.size(); idtr++ ) {
            std::cout << " --------------DT_RES = " << settings.Dt_res[idtr] << "--------------"<< std::endl;

            for (int itr = 0; itr < settings.t_res.size(); itr++) {
                std::cout << "Computing residuals at time t = " << settings.t_res[itr] << std::endl;

                for (int ncons = 0; ncons < nC; ncons++) {

                    Eigen::MatrixXd coef_t(3, Nm[ncons]);

                    std::vector<double> tr(1);
                    std::vector<double> t_evaluate = {settings.t_res[itr] - 2.0 * settings.Dt_res[idtr],
                                                      settings.t_res[itr] - settings.Dt_res[idtr],
                                                      settings.t_res[itr]};

                    for (int j = 0; j < 3; j++) {
                        tr[0] = t_evaluate[j];
                        for (int i = 0; i < Nm[ncons]; i++)
                            surr_coefs[ncons][i].evaluate(tr, coef_t(j, i));
                    }

                //    }
                    Eigen::MatrixXd Sig = Eigen::MatrixXd::Zero(Nm[ncons], Nm[ncons]);
                    for (int i = 0; i < Nm[ncons]; i++)
                        Sig(i, i) = std::sqrt(lambda[ncons](i));
                    Sn_Cons_time.middleRows(ncons * Nr, Nr) = Phi[ncons].leftCols(Nm[ncons]) * Sig * coef_t.transpose();
                }

                if (settings.flag_mean == "IC") {
                    for (int it = 0; it < 3; it++)
                        Sn_Cons_time.col(it) += Ic;
                }

                //Launching SU2_DTR and saving errors and Residuals to file
                int iter = std::round(settings.t_res[itr]/settings.Dt_cfd);
                Write_Restart_Cons_Time(Sn_Cons_time, Coords, settings.out_file, iter, nC, settings.alpha, settings.beta, binary, "POD");
                SU2_DTR(settings, su2_conf, "POD", idtr, itr);
                Write_History_ResError(settings, "POD", idtr, itr);
                std::cout << std::endl;
            }
        }
    }

    std::vector<std::string>::iterator itPOD_STEADY;
    itPOD_STEADY = std::find(methods.begin(),methods.end(),"POD_STEADY");
    if (itPOD_STEADY != methods.end()){
        std::cout << "--------------------------------------" << std::endl ;
        std::cout << "-----Performin POD-steady ResEval-----" << std::endl ;
        std::cout << "--------------------------------------" << std::endl ;

        Eigen::MatrixXd Sn_Cons_time = Eigen::MatrixXd::Zero(nC*Nr, 3);
        std::vector<Eigen::MatrixXd> Phi(nC);
        std::vector<Eigen::VectorXd> lambda(nC);
        std::vector<Eigen::MatrixXd> coefs = {};

        std::vector< std::vector<rbf> > surr_coefs(nC);
        Eigen::VectorXd K_pc(settings.Ns);
        Eigen::MatrixXd eig_vec(settings.Ns, settings.Ns);
        std::vector<int> Nm(nC);
        int N_notZero;
        //Check only for POD for now
        for (int i = 0; i < nC; i++) {
            Phi[i] = Eigen::MatrixXd::Zero(Nr,settings.Ns);
            lambda[i] = Eigen::VectorXd::Zero(settings.Ns);
        }
        std::cout << std::endl << "Extraction of the basis" << std::endl << std::endl;
        for ( int ncons = 0; ncons < nC; ncons ++ ) {
            std::cout << "Processing conservative variable " << ncons << std::endl;
            Phi[ncons] = SPOD_basis(sn_set.middleRows(ncons * Nr, Nr),
                                    lambda[ncons], K_pc, eig_vec,
                                    Nf[nfj],
                                    settings.flag_bc,
                                    settings.flag_filter,
                                    settings.sigma);
            N_notZero = Phi[ncons].cols();
            if (settings.r == 0) Nm[ncons] = Nmod(settings.En, K_pc);
            else Nm[ncons] = std::min(settings.r, N_notZero);
            std::cout << "Number of modes used in reconstruction " << Nm[ncons] << std::endl;
            coefs.push_back(eig_vec.transpose());
            //surr_coefs[ncons] = getSurrCoefs(t_vec, eig_vec, settings.flag_interp);
        }
        std::cout << std::endl;
        Eigen::MatrixXd Field(Nr, nC);
        for (int nt=0; nt< design_space_test.rows(); nt++){
            std::cout << "RECONSTRUCTION FOR MACH=" << design_space_test(nt,0) << "\t AND ALPHA=" //nt instead of 1
                      << design_space_test(nt,1) <<std::endl;
            for (int ncons=0; ncons<nC; ncons++) {
                Eigen::MatrixXd Rec = Reconstruction_steady_POD(param,
                                                                K_pc, lambda[ncons], coefs[ncons],
                                                                Phi[ncons], design_space_test(nt,1),
                                                                alpha_norm, mach_norm,
                                                                design_space_test(nt,0),
                                                                Nm[ncons],
                                                                settings.flag_interp);
                Field.col(ncons)= Rec;
                std::cout << "Done for ncons= "<<ncons<< std::endl;

            }
        std::cout<<Field.topRows(5)<<std::endl;

        //Launching SU2_DTR and saving errors and Residuals to file
        Write_Restart_Cons_Time_steady(Field, Coords, settings.out_file, nC, design_space_test(nt,1), settings.beta,
                                           design_space_test(nt,0), binary, "POD_STEADY");
        SU2_DTR_steady(settings, su2_conf, "POD_STEADY", design_space_test(nt,0), design_space_test(nt,1), nt);
        Write_History_ResError_steady(settings, "POD_STEADY", design_space_test(nt,0), design_space_test(nt,1), nt);
        std::cout <<"DONE"<<std::endl;
        }

        std::cout<<"END OF POD_STEADY"<<std::endl;

    }


    std::vector<std::string>::iterator itISOMAP;
    itISOMAP= std::find(methods.begin(),methods.end(),"ISOMAP");
    if (itISOMAP != methods.end()){
        std::cout << "--------------------------------------" << std::endl ;
        std::cout << "-----Performin ISOMAP+I ResEval-------" << std::endl ;
        std::cout << "--------------------------------------" << std::endl ;

        Eigen::MatrixXd snapshot_star_total = Eigen::MatrixXd::Zero(Nr, nC);
        Eigen::MatrixXd sn_set_partial(Nr, sn_set.cols());

        Eigen::VectorXd norma(sn_set.cols());
        Eigen::MatrixXd y(settings.r_isomap*nC, sn_set.cols());
        Eigen::VectorXi k(nC);
        Eigen::VectorXd kruskal_min = Eigen::VectorXd::Ones(nC);

        // LOOP FOR EACH CONSERVATIVE VARIABLE!
        std::cout << "ISOMAP FOR THE SELECTED CONSERVATIVE VARIABLES...";
        for (int q = 0; q < nC; q++) {
            //std::cout << "ISOMAP FOR CONSERVATIVE VARIABLE NUMBER " << "\t" << q << std::endl;
            sn_set_partial = sn_set.block(Nr * q, 0, Nr, sn_set.cols());
            for (int i = 0; i < sn_set.cols(); i++) {
                norma(i) = sn_set_partial.col(i).norm();
                sn_set_partial.col(i) = sn_set_partial.col(i) / norma(i);
            }

            //std::cout << "number of features:\t\t\t\t\t" << sn_set_partial.rows() << std::endl;
            //std::cout << "number of snapshots:\t\t\t\t\t" << sn_set_partial.cols() << std::endl;

            //std::cout << "Preliminary Step --------------->euclidean distances" << std::endl << std::endl;
            Eigen::MatrixXd euclidean = euclidean_distance(sn_set_partial, sn_set_partial.rows(),
                                                           sn_set_partial.cols());
            //std::cout<<"distanze euclidee"<<std::endl<<euclidean<<std::endl;


            Eigen::VectorXd kruskal_stress(settings.neighbors.size());
            Eigen::MatrixXd euclidean_y(euclidean.rows(), euclidean.cols());
            Eigen::VectorXd K_pc(settings.Ns);

            for (int i = 0; i < settings.neighbors.size(); i++) {

                //std::cout << "1° Step-------->grafo delle distanze per k=" << settings.neighbors[i] << std::endl;
                Eigen::MatrixXd graph = KNN(euclidean, sn_set_partial.cols(), settings.neighbors[i]);
                //std::cout<<"Grafo delle distanze"<<std::endl<<graph<<std::endl;


                //std::cout << "2° Step-------->shortest path dijkastra's algorithm per k=" << settings.neighbors[i]
                //         << std::endl;
                Eigen::MatrixXd distance_input = dijkstraalgorithm(graph, sn_set_partial.cols());
                //std::cout<<"Matrix of geodesic distances "<<endl<<distance_input<<std::endl;


                //std::cout << "3° Step-------->Multidimensional scaling for k=" << settings.neighbors[i] << std::endl;
                Eigen::MatrixXd reduced = multidim_scaling(distance_input, sn_set_partial.cols(), settings.r_isomap,
                                                           settings.En, K_pc);
                // std::cout<<"autovettori in uscita dal MDS"<<endl<<reduced<<std::endl;

                //std::cout << "4° Step-------->Euclidean distance for kruskal stress analysis with k="
                //         << settings.neighbors[i] << std::endl;
                Eigen::MatrixXd distance_output = euclidean_distance(reduced, reduced.rows(), reduced.cols());
                //std::cout<<"La matrice delle distanze con dimensiona ridotta è: "<<endl<<distance_output<<std::endl;

                kruskal_stress(i) = kruskal_evaluation(distance_input, distance_output, sn_set_partial.cols());

                if (kruskal_stress(i) < kruskal_min(q)) {
                    for (int j = 0; j < settings.r_isomap; j++) {
                        y.row((q * settings.r_isomap) + j) = reduced.row(j);
                    }
                    //euclidean_y = distance_output;
                    k(q) = settings.neighbors[i];
                    kruskal_min(q) = kruskal_stress(i);
                }
                //std::cout << "kruskal_stress"<<'\t'<< kruskal_stress(i) << std::endl;
            }

        }
        std::cout<<"DONE "<<std::endl;
        //std::cout << "Kruskal stress minimo \t " << kruskal_min << "\t per k=" << k << std::endl;
        //std::cout << y << std::endl<<std::endl<<std::endl;
        std::cout<<"END OF ISOMAP..."<<std::endl<<std::endl;



        std::cout<<"------------START BACK-MAPPING SECTION--------------"<<std::endl<<std::endl;




        for(int nt=0; nt< design_space_test.rows(); nt++) {             //nt<design_space_test.rows()
            std::cout << "RECONSTRUCTION FOR MACH=" << design_space_test(nt,0) << "\t AND ALPHA="
                      << design_space_test(nt,1) <<std::endl;
            for (int q = 0; q < nC; q++) {
                //std::cout << "Conservative variable number" << "\t" << q << std::endl;
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
                std::vector<double> y_star = ISOMAP_surr_coeff(param,
                                                               y_part,
                                                               design_space_test(nt,1),
                                                               design_space_test(nt,0),
                                                               settings.flag_interp,
                                                               alpha_norm, mach_norm);

                //std::cout << "values RBF" << std::endl;
                //for (int i = 0; i < y_star.size(); i++)
                //   std::cout << y_star[i] << "\t";

                // nearest points for reconstruction y*
                Eigen::VectorXd embedding_sorted = Eigen::VectorXd::Zero(y.cols());
                Eigen::VectorXi y_star_index_nearest = embedding_nearest(y_star,
                                                                         y_part,
                                                                         embedding_sorted);

                //std::cout << "embedding distance" << std::endl << embedding_sorted.transpose() << std::endl;
                //std::cout << "index of embedding neighbors" << std::endl << y_star_index_nearest.transpose()
                          //<< std::endl;


                // EVALUATION OF PENALTY TERM
                Eigen::MatrixXd C = Eigen::MatrixXd::Zero(settings.k_rec, settings.k_rec);
                penalty_term(C, settings.k_rec, embedding_sorted);
                //std::cout<<C<<std::endl;

                //EVALUATION OF GRAM MATRIX
                Eigen::MatrixXd G(settings.k_rec, settings.k_rec);
                gram_matrix(G, settings.k_rec, y_star_index_nearest, y_part, y_star);
                //std::cout<<"Gram Matrix"<<std::endl<<G<<std::endl;

                Eigen::MatrixXd G_tilde = G + C;
                Eigen::VectorXd weigths;
                weigths = lagrange(G_tilde);
                //std::cout << "vector of weights" << std::endl << weigths.transpose() << std::endl;


                Eigen::MatrixXd snapshot_star = Eigen::MatrixXd::Zero(sn_set_partial.rows(), 1);
                for (int index = 0; index < settings.k_rec; index++) {
                    int snap = y_star_index_nearest(index);
                    snapshot_star = snapshot_star + (weigths(index) * sn_set_partial.col(snap)* norma(snap));
                }

                snapshot_star_total.col(q) = snapshot_star;
            }
            std::cout<<snapshot_star_total.topRows(5)<<std::endl;

            //Launching SU2_DTR and saving errors and Residuals to file
            Write_Restart_Cons_Time_steady(snapshot_star_total, Coords, settings.out_file, nC, design_space_test(nt,1), settings.beta,
                                           design_space_test(nt,0),binary, "ISOMAP");
            SU2_DTR_steady(settings, su2_conf, "ISOMAP", design_space_test(nt,0), design_space_test(nt,1), nt);
            Write_History_ResError_steady(settings, "ISOMAP", design_space_test(nt,0), design_space_test(nt,1), nt);
            std::cout << std::endl;
            std::cout <<"DONE"<<std::endl;

        }


        std::cout <<"END OF BACK-MAPPING  "<<std::endl;

    }*/



    std::cout << "MODES Adaptive Residual Evaluation ends" << std::endl << std::endl;
    return 0;

}
