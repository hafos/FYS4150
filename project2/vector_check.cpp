#include <iostream>
#include <armadillo>
#include <cmath>

using namespace arma;

bool check_eigenvectors(const arma::mat& U, const arma::mat& V, const double& tolerance){
  // This function checks wether the eigenvectors in U and V are similar.
  bool sim = 1; // Wether all the vectors are the similar
  int vecmatched = 0; // Count how many of U have been matched with columns in V
  arma::vec order(V.n_cols, fill::zeros); // The order of V that fits to U
  arma::vec inVismatched(V.n_cols, fill::zeros); // Counts which columns in V have been matced
  for (int i=0; i < U.n_cols; i++){ // Go over column i in U
    for (int j=0; j < V.n_cols; j++){ // Go over column j in V
      if (inVismatched(j)==0){ // Only check the columns in V that are not matched yet
        if (approx_equal(U.col(i), V.col(j), "absdiff", tolerance)){
          inVismatched(j) = 1;
          order(i) = j;
          vecmatched = vecmatched + 1;
          break;
        }
        else if (approx_equal(U.col(i), -V.col(j), "absdiff", tolerance)){
          inVismatched(j) = 1;
          order(i) = j;
          vecmatched = vecmatched + 1;
          break;
        }
      }
     }
    }
  // Check the number of matched columns in U
  if (vecmatched != U.n_cols){
    sim = 0;
  }
  // Print order of vectors:
  std::cout << "Vector order: ";
  for (int i=0; i < order.n_elem; i++){
    std::cout << order(i) << " ";
  }
  std::cout << endl;
  return sim;
}

bool check_eigenvalues(const arma::vec& U, const arma::vec& V, const double& tolerance){
  bool sim = 1;
  int vecmatched = 0;
  arma::vec order(V.n_elem, fill::zeros);
  arma::vec inVismatched(V.n_elem, fill::zeros);
  for (int i=0; i<U.n_elem; i++){
    for (int j=0; j<V.n_elem; j++){
      //std::cout<< i << " " << j << endl;
      if (inVismatched(j)==0){
        if (std::abs(U(i)-V(j)) <= tolerance){
          //std::cout << U(i) << " " << V(j) << endl;
          inVismatched(j) = 1;
          order(i) = j;
          vecmatched += 1;
          break;
        }
      }
    }
    //inVismatched.print();
  }
  if (vecmatched != U.n_elem){
    sim = 0;
  }
  // Print order of values:
  std::cout << "Value order: ";
  for (int i=0; i<order.n_elem; i++){
    std::cout << order(i) << " ";
  }
  std::cout << endl;
  return sim;
}
