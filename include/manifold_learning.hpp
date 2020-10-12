#ifndef ISOMAP_HPP
#define ISOMAP_HPP

#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>


#include "Eigen/Dense"
#include "Eigen/Eigenvalues"
#include "Extract_Basis.hpp"

using namespace std;
using namespace Eigen;

MatrixXd isomap(MatrixXd &snap, int f, int s );

MatrixXd euclidean_distance(MatrixXd &snapshots, int righe, int colonne);


MatrixXd KNN(MatrixXd &link, int s, int k);
void sort_vec(VectorXd &vec, VectorXd &sorted_vec, VectorXi &ind);



void init(int nodo_corrente, int V, VectorXd &parent, VectorXd &dist);
int getNearest(VectorXd nodo_visitato_ciclo, int V, VectorXd &dist);
void dijkstra(VectorXd nodo_visitato, int V, VectorXd &dist, MatrixXd &graph_distance, VectorXd &parent);
MatrixXd dijkstraalgorithm(MatrixXd& graph_distance, int V);


inline void extractNLargestEigens(unsigned n, VectorXd &S, MatrixXd &V);
MatrixXd multidim_scaling(MatrixXd& distance, int s, int d, double En, VectorXd K_pc);

double kruskal_evaluation(MatrixXd &distance_input, MatrixXd &distance_output, int s) ;




void ISOMAP(MatrixXd &sn_set, prob_settings &settings, int Nr, int q, VectorXd &K_pc, VectorXi &k, VectorXd &kruskal_min, MatrixXd &y);





VectorXi embedding_nearest(vector<double> &y_star, MatrixXd &y, VectorXd &embedding_sorted);


void penalty_term(MatrixXd &C, int k_reconstruction, VectorXd &embedding_sorted);
void gram_matrix(MatrixXd &G, int k_reconstruction, VectorXi &y_star_index_nearest, MatrixXd &y, vector<double> &y_star);


VectorXd lagrange(MatrixXd &G_tilde);


MatrixXd create_rec_fields(prob_settings &settings, std::vector<double> &y_star, Eigen::MatrixXd &y_part,
                           Eigen::MatrixXd &sn_set_partial, Eigen::VectorXd &norma);

#endif //ISOMAP_HPP
