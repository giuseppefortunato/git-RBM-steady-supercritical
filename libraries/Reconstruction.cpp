#include "Reconstruction.hpp"



#pragma clang diagnostic push
#pragma ide diagnostic ignored "OCUnusedGlobalDeclarationInspection"
smartuq::surrogate::RBF_FUNCTION get_key_rbf ( const std::string &key_string ) {

    if (key_string == "LINEAR")
        return smartuq::surrogate::LINEAR;
    else if (key_string == "CUBIC")
        return smartuq::surrogate::CUBIC;
    else if (key_string == "GAUSSIAN")
        return smartuq::surrogate::GAUSSIAN;
    else if (key_string == "THIN_PLATE")
        return smartuq::surrogate::THIN_PLATE;
    else if (key_string == "MULTIQUADRATICS")
        return smartuq::surrogate::MULTIQUADRATICS;

    return smartuq::surrogate::LINEAR;

}


std::vector<double> ISOMAP_surr_coeff( const Eigen::MatrixXd &param,
                                       const Eigen::RowVectorXd &norme,
                                       const Eigen::MatrixXd &Coeffs,
                                       const double param_rec_1,
                                       const double param_rec_2,
                                       const double param_rec_3,
                                       const double param_rec_4,
                                       const double param_rec_5,
                                       const double param_rec_6,
                                       const int Np,
                                       std::string flag_interp,
                                       const Eigen::RowVectorXd supercritical) {

    std::vector <std::vector<double>> T(param.rows(), std::vector<double>(param.cols()));
    //for supercritical
    std::vector<double> all_rec(Np);
    for (int i=0; i<Np; i++)
        all_rec[i]=supercritical(i);

    /*std::vector<double> all_rec(6);
    all_rec[0] = param_rec_1;
    all_rec[1] = param_rec_2;
    all_rec[2] = param_rec_3;
    all_rec[3] = param_rec_4;
    all_rec[4] = param_rec_5;
    all_rec[5] = param_rec_6;*/
    for (int i = 0; i < param.rows(); i++) {
        for(int j=0; j< param.cols(); j++)
            T[i][j] = param(i, j)/ norme(j);

    }
    //std::cout<<"inside RBF function, T="<<std::endl;
    //std::cout<<T<<std::endl;

    std::vector<double> t(Np);
        for (int i=0; i<Np; i++)
            t[i] = all_rec[i]/ norme(i);

    //Vector of surrogate coefficients
    std::vector<double> coefs_intrp(Coeffs.rows());

    double avgDt=0.0;

    // Create surrogates for coefficients
    std::vector <rbf> surr_coefs = {};
    RBF_CONSTANTS rbf_const{avgDt, 0.0};

    for (int i = 0; i < Coeffs.rows(); i++) {

        std::vector<double> coefs;
        for (int j = 0; j < param.rows(); j++)
            coefs.push_back(Coeffs(i, j));

        surr_coefs.push_back(rbf(T, coefs, get_key_rbf(flag_interp), rbf_const));
        surr_coefs[i].build();
        surr_coefs[i].evaluate(t, coefs_intrp[i]);
    }

    return coefs_intrp;

}




Eigen::MatrixXd Reconstruction_POD ( const Eigen::MatrixXd &param,
                                     const Eigen::RowVectorXd &norme,
                                     const Eigen::VectorXd &K_pc,
                                     const Eigen::VectorXd &lam,
                                     const Eigen::MatrixXd &Coeffs,
                                     const Eigen::MatrixXd &phi,
                                     const double param_rec_1,
                                     const double param_rec_2,
                                     const double param_rec_3,
                                     const double param_rec_4,
                                     const double param_rec_5,
                                     const double param_rec_6,
                                     int Nrec,
                                     int Np,
                                     std::string flag_interp,
                                     const Eigen::RowVectorXd &supercritical) {

    std::vector<std::vector<double> > T(param.rows(), std::vector<double>(param.cols()));
    //for supercritical
    //std::vector<double> all_rec(Np);
    //for (int i=0; i<Np; i++)
    //    all_rec[i]=supercritical(i);

    std::vector<double> all_rec(6);
    all_rec[0] = param_rec_1;
    all_rec[1] = param_rec_2;
    all_rec[2] = param_rec_3;
    all_rec[3] = param_rec_4;
    all_rec[4] = param_rec_5;
    all_rec[5] = param_rec_6;

    double avgDt = 0.0;
    for (int i = 0; i < param.rows(); i++) {
        for (int j = 0; j < param.cols(); j++) {
            T[i][j] = param(i, j) / norme(j);
        }
    }


    std::vector<double> t(Np);
        for (int i = 0; i < Np; i++)
            t[i] = all_rec[i] / norme(i);

        //for ( int i = 1; i < t_vec.size(); i++ ) {
        //    avgDt += t_vec[i] - t_vec[i-1];
        //}

        //avgDt = avgDt/(double)(t_vec.size()-1);

        //Vector of surrogate coefficients
        std::vector<double> coefs_intrp(Nrec);

        // Create surrogates for coefficients
        std::vector<rbf> surr_coefs = {};
        RBF_CONSTANTS rbf_const{avgDt, 0.0};


        for (int i = 0; i < Nrec; i++) {

            std::vector<double> coefs;
            for (int j = 0; j < param.rows(); j++)
                coefs.push_back(Coeffs(i, j));

            surr_coefs.push_back(rbf(T, coefs, get_key_rbf(flag_interp), rbf_const));
            surr_coefs[i].build();
            surr_coefs[i].evaluate(t, coefs_intrp[i]);
        }

        Eigen::VectorXd coefs_t(Nrec);

        for (int i = 0; i < Nrec; i++)
            coefs_t(i) = coefs_intrp[i];

        Eigen::MatrixXd Sig = Eigen::MatrixXd::Zero(Nrec, Nrec);

        for (int i = 0; i < Nrec; i++)
            Sig(i, i) = std::sqrt(lam(i));

        Eigen::MatrixXd Rec_field = phi.leftCols(Nrec) * Sig * coefs_t.head(Nrec);

        return Rec_field;

    }


    Eigen::MatrixXcd Reconstruction_DMD(const double time, const double dt,
                                        const Eigen::VectorXcd &alfa,
                                        const Eigen::MatrixXcd &Phi,
                                        const Eigen::VectorXcd &lam,
                                        const std::string flag_prob) {

        int Nm = lam.size();
        Eigen::VectorXcd omega(Nm);

        for (int i = 0; i < Nm; i++) {
            omega(i) = std::log(lam(i)) / dt;
        }

        Eigen::MatrixXcd v_rec = Eigen::MatrixXcd::Zero(Phi.rows(), 1);

        for (int i = 0; i < Nm; i++) {
            v_rec += std::exp(omega(i) * time) * alfa(i) * Phi.col(i);
        }

        if (flag_prob == "SCALAR") {
            return v_rec;
        } else if (flag_prob == "VECTOR-2D" || flag_prob == "VELOCITY-2D") {
            Eigen::MatrixXcd Rec(Phi.rows() / 2, 2);
            Rec.col(0) = v_rec.topRows(Phi.rows() / 2);
            Rec.col(1) = v_rec.bottomRows(Phi.rows() / 2);
            return Rec;
        } else if (flag_prob == "VECTOR-3D" || flag_prob == "VELOCITY-3D") {
            Eigen::MatrixXcd Rec(Phi.rows() / 3, 3);
            Rec.col(0) = v_rec.topRows(Phi.rows() / 3);
            Rec.col(1) = v_rec.middleRows(Phi.rows() / 3, Phi.rows() / 3);
            Rec.col(2) = v_rec.bottomRows(Phi.rows() / 3);
            return Rec;
        } else {
            std::cout << "Set well problem flag! Exiting ... " << std::endl;
            exit(EXIT_FAILURE);
        }

        return v_rec;
    }

    Eigen::MatrixXd TimeEvo_SPOD(const std::vector<double> &t_vec,
                                 const Eigen::VectorXd &time,
                                 const Eigen::MatrixXd &Coefs,
                                 const Eigen::MatrixXd &Phi,
                                 const Eigen::VectorXd &lam,
                                 const std::string &flag_interp) {

        int Nm = Phi.cols();
        int Nrec = time.size();
        Eigen::MatrixXd Coefs_interp(Nrec, Nm);
        std::vector<rbf> surr_coefs = getSurrCoefs(t_vec, Coefs, flag_interp);


        for (int i = 0; i < Nrec; i++) {
            std::vector<double> t(1, time(i));

            for (int j = 0; j < Nm; j++)
                surr_coefs[j].evaluate(t, Coefs_interp(i, j));

        }

        Eigen::MatrixXd Sig = Eigen::MatrixXd::Zero(Nm, Nm);

        for (int i = 0; i < Nm; i++)
            Sig(i, i) = sqrt(lam(i));

        return Phi * Sig * Coefs_interp.transpose();

    }


    Eigen::MatrixXcd TimeEvo_DMD(Eigen::VectorXd &time,
                                 double dt,
                                 const Eigen::VectorXcd &alfa,
                                 const Eigen::MatrixXcd &Phi,
                                 const Eigen::VectorXcd &lam) {

        int Nm = lam.size();
        Eigen::VectorXcd omega(Nm);

        for (int i = 0; i < Nm; i++) {
/*Remove after check*/          //    std::cout << " Lambda_" << i << " = " << lam(i) << "\t"; 
            omega(i) = std::log(lam(i)) / dt;
/*Remove after check*/          //    std::cout << " omega_" << i << " = " << omega(i) << std::endl;
        }

        Eigen::MatrixXcd v_rec = Eigen::MatrixXcd::Zero(Phi.rows(), time.size());

        for (int k = 0; k < time.size(); k++) {
            for (int i = 0; i < Nm; i++) {
                v_rec.col(k) += std::exp(omega(i) * time(k)) * alfa(i) * Phi.col(i);
                // std::cout << " exp(omega_" << i << "*t) : " << std::exp(omega(i)*time(k)) << std::endl;
            }
        }

        return v_rec;

    }


    Eigen::MatrixXd TimeEvo_RDMD(const std::vector<double> &t_vec,
                                 const Eigen::VectorXd &time,
                                 const Eigen::MatrixXd &Coefs,
                                 const Eigen::MatrixXd &Phi,
                                 const std::string &flag_interp) {

        int Nm = Phi.cols();
        int Nrec = time.size();
        Eigen::MatrixXd Coefs_interp(Nrec, Nm);
        std::vector<rbf> surr_coefs = getSurrCoefs(t_vec, Coefs, flag_interp);


        for (int i = 0; i < Nrec; i++) {
            std::vector<double> t(1, time(i));

            for (int j = 0; j < Nm; j++)
                surr_coefs[j].evaluate(t, Coefs_interp(i, j));

        }

        return Phi * Coefs_interp.transpose();

    }


    void nodes_mrDMD_sort(std::vector<node_mrDMD> &nodes) {

        int swap_count = 1;
        int temp;
        node_mrDMD temp_node;

        while (swap_count > 0) {

            swap_count = 0;

            for (int index = 1; index < nodes.size(); index++) {

                if (nodes[index].l < nodes[index - 1].l) {

                    temp_node = nodes[index - 1];
                    nodes[index - 1] = nodes[index];
                    nodes[index] = temp_node;

                    swap_count++;

                }
            }
        }

    }


    Eigen::MatrixXcd Reconstruction_mrDMD(const double t_inst, const double dts,
                                          const std::vector<node_mrDMD> &nodes,
                                          const std::string flag_prob) {

        Eigen::MatrixXcd v_rec = Eigen::MatrixXcd::Zero(nodes[0].Modes.rows(), 1);
        double dt;
        std::complex<double> omegaj;

        for (int i = 0; i < nodes.size(); i++) {

            if ((t_inst >= nodes[i].t_begin) && (t_inst < nodes[i].t_end)) {

                dt = dts * (double) nodes[i].step;
                for (int j = 0; j < nodes[i].lam.size(); j++) {

                    omegaj = std::log(nodes[i].lam(j)) / dt;
                    v_rec += std::exp(omegaj * (t_inst - nodes[i].t_begin)) * nodes[i].Coefs(j) * nodes[i].Modes.col(j);
                    // v_rec += std::exp(omegaj*t_inst)*nodes[i].Coefs(j)*nodes[i].Modes.col(j);
                }

            }

        }

        if (flag_prob == "SCALAR") {

            return v_rec;

        } else if (flag_prob == "VECTOR-2D" || flag_prob == "VELOCITY-2D") {

            Eigen::MatrixXcd Rec(v_rec.rows() / 2, 2);
            Rec.col(0) = v_rec.topRows(v_rec.rows() / 2);
            Rec.col(1) = v_rec.bottomRows(v_rec.rows() / 2);
            return Rec;

        } else if (flag_prob == "VECTOR-3D" || flag_prob == "VELOCITY-3D") {

            Eigen::MatrixXcd Rec(v_rec.rows() / 3, 3);
            Rec.col(0) = v_rec.topRows(v_rec.rows() / 3);
            Rec.col(1) = v_rec.middleRows(v_rec.rows() / 3, v_rec.rows() / 3);
            Rec.col(2) = v_rec.bottomRows(v_rec.rows() / 3);
            return Rec;

        } else {

            std::cout << "Set well problem flag! Exiting ... " << std::endl;
            exit(EXIT_FAILURE);

        }

    }


    std::vector<rbf> getSurrCoefs(const std::vector<double> &t_vec,
                                  const Eigen::MatrixXd &Coeffs,
                                  std::string flag_interp) {
        int Ns = Coeffs.rows();
        int Nrec = Coeffs.cols();

        std::vector<std::vector<double> > T(t_vec.size(), std::vector<double>(1));
        double avgDt = 0.0;

        for (int i = 0; i < t_vec.size(); i++) {
            T[i][0] = t_vec[i];
        }

        for (int i = 1; i < t_vec.size(); i++) {
            avgDt += t_vec[i] - t_vec[i - 1];
        }

        avgDt = avgDt / (double) (t_vec.size() - 1);

        // Create surrogates for coefficients
        std::vector<rbf> surr_coefs{};
        RBF_CONSTANTS rbf_const{avgDt, 0.0};

        for (int i = 0; i < Nrec; i++) {

            std::vector<double> coefs{};
            for (int j = 0; j < t_vec.size(); j++)
                coefs.push_back(Coeffs(j, i));

            surr_coefs.push_back(rbf(T, coefs, get_key_rbf(flag_interp), rbf_const));
            surr_coefs[i].build();

        }

        return surr_coefs;

    }


    Eigen::MatrixXcd Reconstruction_Hybrid_DMD(const double time,
                                               const std::vector<double> t_vec,
                                               const Eigen::MatrixXcd &alfa,
                                               const Eigen::MatrixXcd &Phi,
                                               const Eigen::VectorXcd &omega,
                                               const std::string flag_prob,
                                               const std::string flag_interp) {

        int Nm = omega.size();


        Eigen::MatrixXcd v_rec = Eigen::MatrixXcd::Zero(Phi.rows(), 1);
        Eigen::MatrixXd Coefs_real = alfa.real();
        Eigen::MatrixXd Coefs_imag = alfa.imag();
        std::vector<rbf> surr_coefs_real = getSurrCoefs(t_vec, Coefs_real, flag_interp);
        std::vector<rbf> surr_coefs_imag = getSurrCoefs(t_vec, Coefs_imag, flag_interp);

        std::vector<double> t(1, time);
        for (int i = 0; i < Nm; i++) {
            double R, I;
            surr_coefs_real[i].evaluate(t, R);
            surr_coefs_imag[i].evaluate(t, I);
            std::complex<double> coef_t(R, I);

            v_rec += std::exp(omega(i) * time) * coef_t * Phi.col(i);

        }

        if (flag_prob == "SCALAR") {

            return v_rec;

        } else if (flag_prob == "VECTOR-2D" || flag_prob == "VELOCITY-2D") {

            Eigen::MatrixXcd Rec(Phi.rows() / 2, 2);
            Rec.col(0) = v_rec.topRows(Phi.rows() / 2);
            Rec.col(1) = v_rec.bottomRows(Phi.rows() / 2);
            return Rec;

        } else if (flag_prob == "VECTOR-3D" || flag_prob == "VELOCITY-3D") {

            Eigen::MatrixXcd Rec(Phi.rows() / 3, 3);
            Rec.col(0) = v_rec.topRows(Phi.rows() / 3);
            Rec.col(1) = v_rec.middleRows(Phi.rows() / 3, Phi.rows() / 3);
            Rec.col(2) = v_rec.bottomRows(Phi.rows() / 3);
            return Rec;

        } else {

            std::cout << "Set well problem flag! Exiting ... " << std::endl;
            exit(EXIT_FAILURE);

        }

    }


    Eigen::MatrixXd Reconstruction_RDMD(const double time,
                                        const std::vector<double> t_vec,
                                        const Eigen::MatrixXd &alfa,
                                        const Eigen::MatrixXd &Phi,
                                        const std::string flag_prob,
                                        const std::string flag_interp) {

        int Nm = alfa.rows();
        std::vector<rbf> surr_coefs = getSurrCoefs(t_vec, alfa.transpose(), flag_interp);
        Eigen::VectorXd coef_t(Nm);

        std::vector<double> t(1, time);
        for (int i = 0; i < Nm; i++)
            surr_coefs[i].evaluate(t, coef_t(i));

        Eigen::MatrixXd v_rec = Phi * coef_t;

        if (flag_prob == "SCALAR") {
            return v_rec;
        } else if (flag_prob == "VECTOR-2D" || flag_prob == "VELOCITY-2D") {

            Eigen::MatrixXd Rec(Phi.rows() / 2, 2);
            Rec.col(0) = v_rec.topRows(Phi.rows() / 2);
            Rec.col(1) = v_rec.bottomRows(Phi.rows() / 2);
            return Rec;
        } else if (flag_prob == "VECTOR-3D" || flag_prob == "VELOCITY-3D") {

            Eigen::MatrixXd Rec(Phi.rows() / 3, 3);
            Rec.col(0) = v_rec.topRows(Phi.rows() / 3);
            Rec.col(1) = v_rec.middleRows(Phi.rows() / 3, Phi.rows() / 3);
            Rec.col(2) = v_rec.bottomRows(Phi.rows() / 3);
            return Rec;
        } else {

            std::cout << "Set well problem flag! Exiting ... " << std::endl;
            return Eigen::MatrixXd::Zero(1, 1);
            exit(EXIT_FAILURE);
        }

    }


    Eigen::MatrixXd Reconstruction_DMD_Interp(const double time,
                                              const std::vector<double> t_vec,
                                              const Eigen::MatrixXcd &alfa,
                                              const Eigen::MatrixXcd &Phi,
                                              const std::string flag_prob,
                                              const std::string flag_interp) {

        int Nm = alfa.rows();
        std::vector<rbf> surr_coefs_real = getSurrCoefs(t_vec, alfa.real().transpose(), flag_interp);
        std::vector<rbf> surr_coefs_imag = getSurrCoefs(t_vec, alfa.imag().transpose(), flag_interp);
        Eigen::VectorXcd coef_t(Nm);

        std::vector<double> t(1, time);
        for (int i = 0; i < Nm; i++) {
            double tmp_i, tmp_r;
            surr_coefs_real[i].evaluate(t, tmp_r);
            surr_coefs_imag[i].evaluate(t, tmp_i);
            std::complex<double> c(tmp_r, tmp_i);
            coef_t(i) = c;
        }

        Eigen::MatrixXcd v_rec = Phi * coef_t;

        if (flag_prob == "SCALAR") {
            return v_rec.real();
        } else if (flag_prob == "VECTOR-2D" || flag_prob == "VELOCITY-2D") {

            Eigen::MatrixXcd Rec(Phi.rows() / 2, 2);
            Rec.col(0) = v_rec.topRows(Phi.rows() / 2);
            Rec.col(1) = v_rec.bottomRows(Phi.rows() / 2);
            return Rec.real();
        } else if (flag_prob == "VECTOR-3D" || flag_prob == "VELOCITY-3D") {

            Eigen::MatrixXcd Rec(Phi.rows() / 3, 3);
            Rec.col(0) = v_rec.topRows(Phi.rows() / 3);
            Rec.col(1) = v_rec.middleRows(Phi.rows() / 3, Phi.rows() / 3);
            Rec.col(2) = v_rec.bottomRows(Phi.rows() / 3);
            return Rec.real();
        } else {

            std::cout << "Set well problem flag! Exiting ... " << std::endl;
            exit(EXIT_FAILURE);
        }

        return Eigen::MatrixXd::Zero(1, 1);

    }


#pragma clang diagnostic pop
