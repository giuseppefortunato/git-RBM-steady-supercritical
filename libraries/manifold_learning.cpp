#include "manifold_learning.hpp"
#include "Extract_Basis.hpp"
#include <cmath>


using namespace std;
using namespace Eigen;



// function to calculate euclidean distance.
MatrixXd euclidean_distance(MatrixXd &snapshots, int righe, int colonne){

    MatrixXd  weight = MatrixXd::Zero(colonne,colonne);
    for(int snap=0; snap<righe ; snap++) {                                                          //ciclo per processare tutte le righe della matrice
        VectorXd features=snapshots.row(snap);
        for (int i = 0; i < colonne; i++) {
            for (int j = 0; j < colonne; j++) {
                double distance = 0;
                distance = (features(i) - features(j))*(features(i) - features(j));
                weight(i, j) = weight(i, j) + distance;
            }
        }
    }

    for (int i = 0; i < colonne; i++) {
        for (int j = 0; j < colonne; j++)
            weight(i, j) = std::sqrt(weight(i, j));
    }



    return weight;
}

// functions to perform KNN algorithm
MatrixXd KNN(MatrixXd &link, int s, int k){


    MatrixXd distance(s, s) ;
    for(int i=0; i < s; i++) {
        for (int j = 0; j < s; j++)
            distance(i,j) = numeric_limits<double>::infinity();
    }

    for(int  i=0; i < s; i++ ){
        VectorXd x(s);
        x = link.row(i);
        VectorXi ind;
        VectorXd sorted;
        sort_vec(x, sorted, ind);
        for (int neigh=0; neigh< k+1; neigh++ ){
            distance(i,ind(neigh)) = sorted(neigh);
        }
    }
     for(int i=0; i<s; i++){
         for(int j=0; j<s; j++){
             if(i!=j){
                 if(distance(i,j)!= numeric_limits<double>::infinity() && distance(j,i)!=distance(i,j))
                     distance(j,i)=distance(i,j);}
         }
     }

    //matlab prove to make the graph symmetric

   // MatrixXd trans= distance.transpose();
    //distance =(distance+trans)/2;


    return distance;
}
void sort_vec(Eigen::VectorXd &vec, Eigen::VectorXd &sorted_vec, Eigen::VectorXi &ind){
    ind=VectorXi::LinSpaced(vec.size(),0,vec.size()-1);                                                           //[0 1 2 3 ... N-1]
    auto rule=[vec](int i, int j)->bool{
        return vec(i)<vec(j);
    };                                                                                                                  // regular expression, as a predicate of sort
    std::sort(ind.data(),ind.data()+ind.size(),rule);
    //The data member function returns a pointer to the first element of VectorXd, similar to begin()
    sorted_vec.resize(vec.size());
    for(int i=0;i<vec.size();i++){
        sorted_vec(i)=vec(ind(i));
    }
}


//functions to perform dijkstra algorithm
void init(int nodo_corrente, int V, VectorXd &parent, VectorXd &dist){
    for(int i=0; i<V; i++) {
        parent(i) = i;
        dist(i) = numeric_limits<double>::infinity();
    }
    dist(nodo_corrente) = 0;
}
int getNearest(VectorXd nodo_visitato_ciclo, int V, VectorXd &dist){
    double minvalue = numeric_limits<double>::infinity();
    int minnode = 0 ;
    for(int i=0; i<V; i++){
        if(!nodo_visitato_ciclo(i) && dist(i) < minvalue){
            minvalue = dist(i);
            minnode = i;

        }
    }
    return minnode;
}
void dijkstra(VectorXd nodo_visitato, int V, VectorXd &dist, MatrixXd &graph_distance, VectorXd &parent){
    for(int i=0; i<V; i++){
        int nearest = getNearest(nodo_visitato, V, dist);
        nodo_visitato(nearest) = 1;

        for (int adj =0; adj <V; adj++){
            if(graph_distance(nearest,adj) != numeric_limits<double>::infinity() && dist(adj) > dist(nearest)+graph_distance(nearest,adj) ){
                dist(adj) = dist(nearest)+graph_distance(nearest,adj);
                parent(adj) = nearest;
            }
        }

    }

}
MatrixXd dijkstraalgorithm(MatrixXd& graph_distance, int V) {

    //const int V=9 ;                                             // number of vertices
    //int cost[V][V];                                             // matrix that rappresent our graph
    VectorXd dist(V);                                                // vettore distance a partire da un singolo nodo verso li altri
    int src;                                                    // indicatore del nodo che sto considerando
    VectorXd visited = VectorXd::Zero(V);                                       // inizializzo tutti i nodi cone non visitati...booleano
    VectorXd parent(V);
    MatrixXd final_distance= MatrixXd::Zero(V,V);


    /*cout << "numero di vertici = " << V << endl;
    cout << "la matrice dei costi associata ai vertici è: " << endl;
    for (int i = 0; i < V; i++) {
        cout << endl;
        for (int j = 0; j < V; j++) {
            cout << graph_distance(i,j) << "  ";
        }
    }
    cout << endl << endl;*/

    for(src=0; src<V; src++) {
        //cout<<"algoritmo per il nodo:"<<src<<endl;
        //VectorXd visited = VectorXd::Zero(V);
        init(src, V, parent, dist);
        dijkstra(visited, V, dist, graph_distance, parent);
        for(int i=0; i<V; i++)
            final_distance(src,i)= dist(i);
    }



    return final_distance;
}


//functions to perform multidimensional scaling
inline void extractNLargestEigens(unsigned n, VectorXd &S, MatrixXd &V)
{
    // Note: m is the original dimension

    const unsigned m = S.rows();

    // Copy the original matrix
    const MatrixXd origV = V;

    // Sort by eigenvalue
    constexpr double epsilon = 1e-16;
    vector<pair<double, unsigned>> sortS(m);
    for (unsigned i = 0; i < m; ++ i) sortS[i] = make_pair(max(S(i), epsilon), i);
    partial_sort(sortS.begin(), sortS.begin() + n, sortS.end(), greater<pair<double, unsigned>>());



    // Resize matrices
    S.resize(n);
    V.resize(m, n);

    // Set values
    for (unsigned i = 0; i < n; ++ i)
    {
        S(i)     = sortS[i].first;
        V.col(i) = origV.col(sortS[i].second);
    }


}
MatrixXd multidim_scaling(MatrixXd& distance, int s, int d, double En, VectorXd K_pc){

    const MatrixXd H = MatrixXd::Identity(s, s) - (1.0 / s) * VectorXd::Ones(s) * VectorXd::Ones(s).transpose();           //H = orthogonal projection onto span {1}
    MatrixXd D(s,s);
    for(int i=0; i< s; i++){
        for (int j=0; j<s; j++){
            double tmp=0;
            tmp=distance(i,j)*distance(i,j);
            D(i,j)= tmp;
        }
    }

    MatrixXd B = -0.5*H*D*H;                                       // double centered matrix, rank B give me the embedding dimensional space
    //cout<<"B"<<endl<<B<<endl;
    EigenSolver<MatrixXd> solver(B);
    VectorXd S = solver.eigenvalues().real();
    MatrixXd V = solver.eigenvectors().real();
    extractNLargestEigens(s, S, V);
    //cout<<"autovalori"<<endl<<S.transpose()<<endl;
    //cout<<"autovettori"<<endl<<V<<endl;

    int Nrec;
    if ( d == 0 ) {
            Nrec = Nmod( En, K_pc);
            //std::cout << "Number of dimension for the desired energy content : " << Nrec << std::endl;
        } else {
            int Nmax = V.cols();
            Nrec = std::min(d,Nmax);
            //std::cout << " Number of dimension fixed : " << Nrec << std::endl;
        }

    VectorXd lambda(Nrec) ;
    MatrixXd autovet(s,Nrec);
    for(int i=0; i<Nrec; i++){
        lambda(i)= S(i);
        autovet.col(i)=V.col(i);
    }

    //cout<<"autovalori dopo Nrec "<<endl<<lambda.transpose()<<endl;

    for (int i=0; i<Nrec; i++){
        for(int j=0; j<s; j++){
            autovet(j,i)= autovet(j,i)*(std::sqrt(lambda(i)));
        }
    }

    //verifica normalizzazione--------> verificata
    //int ver=0;
    //ver = autovet.col(0).transpose()*autovet.col(0);
    //cout<<"verific normalizzazione: velore ottenuto "<<ver<<endl;
    //cout<<"lambda associato"<<lambda(0)<<endl;



    return autovet.transpose();
}


//function to evaluate the kruskal-stress
double kruskal_evaluation(MatrixXd &distance_input, MatrixXd &distance_output, int s){
    double num = 0;
    double den = 0;
    double val = 0;
    for (int w = 0; w < s; w++) {
        for (int j = 0; j < s; j++) {
            num = num + ((distance_input(w, j) - distance_output(w, j))*(distance_input(w, j) - distance_output(w, j)));
            den = den + (distance_input(w, j)*distance_input(w,j));
        }

        val=sqrt(num/den);


    }

    return val;
}


//FUNCTION TO PERFORM -ISOMAP-
void ISOMAP(MatrixXd &sn_set, prob_settings &settings, int Nr, int q, VectorXd &K_pc, VectorXi &k, VectorXd &kruskal_min, MatrixXd &y){

    MatrixXd sn_set_partial(Nr, sn_set.cols());
    sn_set_partial = sn_set.block(Nr * q, 0, Nr, sn_set.cols());
    VectorXd norma(sn_set.cols());

    for (int i = 0; i < sn_set.cols(); i++) {
        norma(i) = sn_set_partial.col(i).norm();
        sn_set_partial.col(i) = sn_set_partial.col(i) / norma(i);                       //---> Normalization of the conserved variable
    }

    MatrixXd euclidean = euclidean_distance(sn_set_partial,                  //--> function to compute the euclidean distance between each pair of snapshot.
                                                   sn_set_partial.rows(),              // euclidean is a distance matrix(Ns,Ns) , symmetric with zero on diagonal.
                                                   sn_set_partial.cols());             // Ns number of snapshots.

    VectorXd kruskal_stress(settings.neighbors.size());
    MatrixXd euclidean_y(euclidean.rows(),euclidean.cols());

    for( int i=0; i<settings.neighbors.size(); i++ ){

        MatrixXd graph = KNN(euclidean,                                     // take as input the distance matrix, and the number of neighbors that you want to connect to create your graph.
                             sn_set_partial.cols(),                             //graph is a symmetric matrix of distance, where each number represent the two snapshot connected. If two snapshot are
                             settings.neighbors[i]);                            // not connected, there is anf INF in the matrix.

        MatrixXd distance_input = dijkstraalgorithm(graph,                   //compute all the distance between snapshots following the shortest paths highligthed in the matrix graph.
                                                    sn_set_partial.cols());

        MatrixXd reduced = multidim_scaling(distance_input,                         // Dimensionality reduction of the system.
                                            sn_set_partial.cols(),              //r_isomap is the number of eigenvectors that i want to use
                                            settings.r_isomap,                  // reduced (Ns, r_isomap )
                                            settings.En,
                                            K_pc);

        MatrixXd distance_output = euclidean_distance(reduced,
                                                      reduced.rows(),
                                                      reduced.cols());

        kruskal_stress(i)= kruskal_evaluation(distance_input,
                                              distance_output,
                                              sn_set_partial.cols());

        if (kruskal_stress(i)<kruskal_min(q)) {
            for (int j = 0; j < settings.r_isomap; j++)
                y.row((q * settings.r_isomap) + j) = reduced.row(j);

            k(q) = settings.neighbors[i];
            kruskal_min(q) = kruskal_stress(i);
        }

    }

}









//function to evaluate the nearest point for y*, that represent the position of the snapshot that i want to reconstruct.
VectorXi embedding_nearest(vector<double> &y_star, MatrixXd &y, VectorXd &embedding_sorted) {

    VectorXd embedding_distance= VectorXd::Zero(y.cols());
    for (int i = 0; i < y.rows(); i++) {
        for (int j = 0; j < y.cols(); j++) {
            double val = 0;
            val = (y_star[i] - y(i, j)) * (y_star[i] - y(i, j));
            embedding_distance(j) = embedding_distance(j) + val;
        }
    }

    for (int i = 0; i < embedding_distance.size(); i++){
        double tmp=0;
        tmp=sqrt(embedding_distance(i));
        embedding_distance(i) =  tmp;}

    //std::cout << "embedding_distance" << std::endl << embedding_distance << std::endl;
    VectorXi ind;
    VectorXd sorted;
    sort_vec(embedding_distance, sorted, ind);
    embedding_sorted=sorted;
    //std::cout << "embedding_distance_sorted_from functions" << std::endl << sorted << std::endl;
    //std::cout << "index" << std::endl << ind << std::endl;

    return ind;

}


void penalty_term(MatrixXd &C, int k_reconstruction, VectorXd &embedding_sorted){

    double eps= 0.01;
    double omega= 4;
    for (int i=0; i<k_reconstruction; i++) {
        double val=0;
        val= (embedding_sorted(i)/embedding_sorted(k_reconstruction-1));     //embedding_sorted.rows()-1
        C(i,i)= pow(val,omega)*eps;
    }
}

void gram_matrix(MatrixXd &G, int k_reconstruction, VectorXi &y_star_index_nearest, MatrixXd &y, vector<double> &y_star){

    for(int i=0; i < k_reconstruction; i++)
    {
        int snap1 = y_star_index_nearest(i);                                    //snap assume il valore dello snapshot che mi serve
        VectorXd vett1(y.rows());
        for(int index=0; index<y.rows(); index++) {
            vett1(index) = y_star[index] - y(index, snap1);
        }
        for(int j=0; j < k_reconstruction; j++)
        {
            VectorXd vett2(y.rows());
            int snap2 = y_star_index_nearest(j);
            for(int index=0; index<y.rows(); index++){
                vett2(index)=y_star[index]-y(index,snap2);
            }
            G(i,j)= vett1.transpose()*vett2;
        }
    }

}

VectorXd lagrange(MatrixXd &G_tilde) {

    MatrixXd linear_system = MatrixXd::Zero(G_tilde.rows() + 1, G_tilde.rows() + 1);
    for (int i = 0; i < G_tilde.rows(); i++) {
        for (int j = 0; j <= G_tilde.rows(); j++) {
            linear_system(i, j) = 2 * G_tilde(i, j);
            if (j == G_tilde.rows())
                linear_system(i, j) = 1;
        }
    }

    for (int j = 0; j < G_tilde.rows(); j++)
        linear_system(G_tilde.rows(), j) = 1;


    Eigen::VectorXd rhs = Eigen::VectorXd::Zero(G_tilde.rows() + 1);
    rhs(G_tilde.rows()) = 1;

    //std::cout << "linear_system" << std::endl << linear_system << std::endl;
    //std::cout << "rhs" << std::endl << rhs << std::endl;
    Eigen::VectorXd w(G_tilde.rows() + 1);
    w = linear_system.colPivHouseholderQr().solve(rhs);
    //std::cout << "Vettore pesi" << std::endl << w << std::endl;

    return w;
}

//GENERAL FUNCTION TO COMPUTE THE RECONSTRUCTED FIELDS USING THE FUNCTIONS IMPLEMENTED ABOVE.
MatrixXd create_rec_fields(prob_settings &settings, std::vector<double> &y_star, Eigen::MatrixXd &y_part,
                           Eigen::MatrixXd &sn_set_partial, Eigen::VectorXd &norma){

    //  NEAREST POINT TO RECONSTRUCTION ON y*
    Eigen::VectorXd embedding_sorted = Eigen::VectorXd::Zero(settings.r_isomap);
    Eigen::VectorXi y_star_index_nearest = embedding_nearest(y_star, y_part, embedding_sorted);
    std::cout<<"embedding distance"<<std::endl<<embedding_sorted.transpose()<<std::endl;
    std::cout<<"index of embedding neighbors"<<std::endl<<y_star_index_nearest.transpose()<<std::endl;

    // EVALUATION OF PENALTY TERM
    Eigen::MatrixXd C = Eigen::MatrixXd::Zero(settings.k_rec, settings.k_rec);
    penalty_term(C, settings.k_rec, embedding_sorted);

    //EVALUATION OF GRAM MATRIX
    Eigen::MatrixXd G(settings.k_rec, settings.k_rec);
    gram_matrix(G, settings.k_rec, y_star_index_nearest, y_part, y_star);

    //EVALUATION OF WEIGHTS
    Eigen::MatrixXd G_tilde = G + C;
    Eigen::VectorXd weigths;
    weigths = lagrange(G_tilde);
    std::cout<<"vector of weights "<<std::endl<<weigths.transpose()<<std::endl<<std::endl;

    Eigen::MatrixXd snapshot_star = Eigen::MatrixXd::Zero(sn_set_partial.rows(), 1);
    for (int index = 0; index < settings.k_rec; index++) {
        int snap = y_star_index_nearest(index);
        snapshot_star = snapshot_star + (weigths(index) * sn_set_partial.col(snap) * norma(snap));
    }

    return snapshot_star;
}


