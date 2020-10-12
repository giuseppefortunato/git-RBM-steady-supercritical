#ifndef WRITE_OUTPUTS_HPP
#define WRITE_OUTPUTS_HPP

#include "Reconstruction.hpp"
#include "read_Inputs.hpp"

const int CGNS_STRING_SIZE = 33;

void Config_stream ( prob_settings settings );

void write_modes_sPOD ( const Eigen::MatrixXd &Phi_cut, 
                        const Eigen::MatrixXd &Coords,
                        std::string flag_prob );

void write_modes_DMD ( const Eigen::MatrixXcd &Phi_cut,
                    const Eigen::MatrixXd &Coords, 
                    std::string flag_prob );

void write_alfa_lam_DMD( Eigen::VectorXcd Alfa,
                         Eigen::VectorXcd Lambda);

void write_coeffs_sPOD ( const Eigen::MatrixXd &Coeffs,
                        const std::vector<double> &t_vec,
                        const Eigen::VectorXd &lam );

void write_TimeDynamics_DMD ( Eigen::VectorXcd omega,
                              Eigen::VectorXcd alfa,
                              Eigen::VectorXd t);

void write_CoefsDynamics_mrDMD( std::vector<node_mrDMD> &nodes, 
                                int level,
                                int ns,
                                int max_levels );



string write_reconstruction_file ( const Eigen::MatrixXd &Rec,
                                   const Eigen::MatrixXd &Coords,
                                   int nt,
                                   int nC,
                                   const int Nm,
                                   prob_settings settings);



void Write_Restart_Cons_Time ( const Eigen::MatrixXd &Rec,
                                    const Eigen::MatrixXd &Coords,
                                    std::string filename,
                                    int nt,
                                    int nC,
                                    double alpha,
                                    double beta,
                                    std::string flag,
                                    std::string method);

void Write_Restart_Cons_Time_steady ( const Eigen::MatrixXd &Rec,
                                      const Eigen::MatrixXd &Coords,
                                      std::string filename,
                                      int nC,
                                      double alpha,
                                      double beta,
                                      double mach,
                                      std::string flag,
                                      std::string method);


void write_modes ( const Eigen::MatrixXd &Phi_cut ); //write modes RDMD

void write_coefs ( const Eigen::MatrixXd &Coefs ); //write coefs RDMD

void write_err_j ( Eigen::MatrixXd data, std::string filename ); //write error/jaccard surface for RBM method

void Write_Plot3d_Modes( Eigen::MatrixXd Phi,           //Modes have to be scalar functions
                        std::string filename, 
                        plot3d_info Info );

#endif // WRITE_OUTPUTS_HPP