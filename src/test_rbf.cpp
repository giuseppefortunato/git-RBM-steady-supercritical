//
// test for RADIAL BASIS FUNCTIONS-july 2020
//

#include "Extract_Basis.hpp"
#include "read_Inputs.hpp"
#include "Generate_snset.hpp"
#include "Reconstruction.hpp"
#include "write_Outputs.hpp"
#include "manifold_learning.hpp"

int main(int argc, char *argv[]) {

    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << "-----------Single-MODES-TEST_RBF start-------------" << std::endl << std::endl;

    std::string filecfg = argv[1];
    prob_settings settings;

    //Reading configuration file
    Read_cfg(filecfg, settings);
    Config_stream(settings);

    // Calculate number of grid points
    std::string root_inputfile;
    root_inputfile.assign(settings.in_file, 0, settings.in_file.size() - 4);
    std::string input_format;
    input_format.assign(settings.in_file, settings.in_file.size() - 3, 3);

    std::stringstream buffer;
    buffer << std::setfill('0') << std::setw(5) << std::to_string(settings.nstart);
    std::string file_1 = root_inputfile + "_" + buffer.str() + "." + input_format;

    int Nr = N_gridpoints(file_1);
    std::cout << "Number of grid points : " << Nr << std::endl;

    // Reading coordinates
    std::cout << "Reading Coordinates ... \t ";
    Eigen::MatrixXd Coords = read_col(file_1, Nr, settings.Cols_coords);
    std::cout << "Done " << std::endl;

    // Create matrix of snapshots
    std::cout << "Storing snapshot Matrix ... \n ";
    Eigen::MatrixXd sn_set = generate_snap_matrix(Nr, settings.Ns, settings.Ds, settings.nstart,
                                                  settings.Cols,
                                                  settings.in_file,
                                                  settings.flag_prob);

    Eigen::VectorXd mean = sn_set.rowwise().mean();
//    Eigen::VectorXd Ic = IC(settings, ,Nr)

    if (settings.flag_mean == "YES") {
        std::cout << "Subtracting mean from snapshots ... " << std::endl << std::endl;
        for (int i = 0; i < settings.Ns; i++)
            sn_set.col(i) -= mean;
    } else if (settings.flag_mean == "IC") {
        std::cout << "Mean is the only reference solution implemented in SingleMODES\n Exiting ... " << std::endl
                  << std::endl;
        exit(EXIT_FAILURE);
    }

    Eigen::VectorXd K_pc(settings.Ns);
    if (settings.flag_method[0] == "ISOMAP") {

        std::cout << "-------------------" << std::endl;
        std::cout << "--ISOMAP ANALYSIS--" << std::endl;
        std::cout << "-------------------" << std::endl;

        Eigen::MatrixXd param(settings.Ns, settings.Np);

        std::cout << "Read file of parameters" << std::endl;
        std::string file_parameters = argv[2];
        param = read_param(file_parameters, settings.Ns, settings.Np);
        std::cout << param << std::endl;
        std::cout << "done" << std::endl;

        Eigen::VectorXd norma(sn_set.cols());
        for (int i = 0; i < sn_set.cols(); i++) {
            norma(i) = sn_set.col(i).norm();
            sn_set.col(i) = sn_set.col(i) / norma(i);
        }

        int f = sn_set.rows();                                            //number of features i.e dimension of the problem.
        int s = sn_set.cols();                                            //number of snapshots.


        std::cout << "number of features:\t\t\t\t\t" << f << std::endl;
        std::cout << "number of snapshots:\t\t\t\t\t" << s << std::endl;


        std::cout << "Preliminary Step --------------->euclidean distances" << std::endl << std::endl;
        Eigen::MatrixXd euclidean = euclidean_distance(sn_set, f, s);
        //std::cout<<"distanze euclidee"<<std::endl<<euclidean<<std::endl;

        std::cout << "------------------------------------" << std::endl;
        std::cout << "Starting computation of ISOMAP steps" << std::endl;
        std::cout << "------------------------------------" << std::endl;

        int k = 0;
        double kruskal_min = 1;
        Eigen::VectorXd kruskal_stress(settings.neighbors.size());
        Eigen::MatrixXd y(euclidean.rows(), euclidean.cols());
        Eigen::MatrixXd euclidean_y(euclidean.rows(), euclidean.cols());



        for (int i = 0; i < settings.neighbors.size(); i++) {

            std::cout << "1° Step-------->grafo delle distanze per k=" << settings.neighbors[i] << std::endl;
            Eigen::MatrixXd graph = KNN(euclidean, s, settings.neighbors[i]);
            //std::cout<<"Grafo delle distanze"<<std::endl<<graph<<std::endl;


            std::cout << "2° Step-------->shortest path dijkastra's algorithm per k=" << settings.neighbors[i]
                      << std::endl;
            Eigen::MatrixXd distance_input = dijkstraalgorithm(graph, s);
            //std::cout<<"Matrix of geodesic distances "<<endl<<distance_input<<std::endl;


            std::cout << "3° Step-------->Multidimensional scaling for k=" << settings.neighbors[i] << std::endl;
            Eigen::MatrixXd reduced = multidim_scaling(distance_input, s, settings.r, settings.En, K_pc);
            // std::cout<<"autovettori in uscita dal MDS"<<endl<<reduced<<std::endl;

            std::cout << "4° Step-------->Euclidean distance for kruskal stress analysis with k="
                      << settings.neighbors[i] << std::endl;
            Eigen::MatrixXd distance_output = euclidean_distance(reduced, reduced.rows(), reduced.cols());
            //std::cout<<"La matrice delle distanze con dimensiona ridotta è: "<<endl<<distance_output<<std::endl;

            kruskal_stress(i) = kruskal_evaluation(distance_input, distance_output, s);

            if (kruskal_stress(i) < kruskal_min) {
                y = reduced;
                euclidean_y = distance_output;
                k = settings.neighbors[i];
                kruskal_min = kruskal_stress(i);
            }
            std::cout << "kruskal_stress" << std::endl << kruskal_stress(i) << std::endl << std::endl << std::endl;

        }
        std::cout << "Kruskal stress minimo \t " << kruskal_min << "\t per k=" << k << std::endl;
        std::cout << y << std::endl;

        double alpha_norm = param.col(1).norm();
        double mach_norm = param.col(0).norm();
        std::cout<<"Normalization of param" <<std::endl;
        std::cout<<alpha_norm<<std::endl<<mach_norm<<std::endl;


        Eigen::VectorXd alpha = Eigen::VectorXd::LinSpaced(100,10.0,19.0);
        //std::cout<<alpha<<std::endl;
        Eigen::VectorXd mach = Eigen::VectorXd::LinSpaced(100,0.1,0.4);
        //std::cout<<mach<<std::endl;

        Eigen::MatrixXd design_space_test(alpha.size()*alpha.size(),2);
        int interval = 0;
        for (int i = 0; i < alpha.size(); i++) {
            for (int j = 0; j < mach.size(); j++) {
                design_space_test(interval+j, 0) = mach(i);
                design_space_test(interval+j, 1) = alpha(j);
            }
            interval = interval+alpha.size();
        }

        //std::cout<<design_space_test<<std::endl;

        Eigen::MatrixXd test(alpha.size()*alpha.size(),2);
        for(int nt=0; nt < test.rows(); nt++) {
            //std::cout << "reconstruction using alpha=" << settings.alpha_rec[nt] << "\t and Mach="
            //       << settings.mach_rec[nt] << std::endl;
            std::vector<double> y_star = ISOMAP_surr_coeff(param, y, design_space_test(nt,1),
                                                           design_space_test(nt,0),
                                                           settings.flag_interp, alpha_norm, mach_norm);

            //std::cout << "values RBF" << std::endl;
            for (int i = 0; i < y_star.size(); i++) {
                test(nt, i) = y_star[i];
                // std::cout << y_star[i] << std::endl;
            }
        }
        std::cout << "min e max" << std::endl;
        std::cout << test.minCoeff() << std::endl;
        std::cout << test.maxCoeff() << std::endl;
        std::string format_solution;
        format_solution.assign(settings.out_file, settings.out_file.size() - 3, 3);

        if (format_solution == "csv") {
            std::cout << "Writing RBF interpolation data in csv format ..." << "\t";
            int nt=0;
            write_Reconstructed_fields_csv (test, design_space_test ,
                                           settings.Cols,
                                           settings.out_file,
                                           nt, test.cols());
                                // remenmber to modify the function adding a aline
            std::cout << "Done" << std::endl << std::endl;
        }
    }

    return 0;
}

