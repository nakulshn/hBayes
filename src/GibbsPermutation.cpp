#include <RcppArmadillo.h>
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]

using namespace arma;
using namespace Rcpp;
using namespace arma;

double binomaloglik_(const arma::vec& Y, const arma::vec& Xbeta, const arma::vec& expXbeta)
{
    return sum(Y % Xbeta - log1p(expXbeta));
}

double normalloglik(const arma::vec& res, const double sigma)
{
    double lik = - ((double) res.size() )* log(sigma);
    lik += -0.5 *  sum(res % res)/pow(sigma,2);
    return lik;
}


// [[Rcpp::export]]
Rcpp::List permute_gibbs_beta_normal_cpp(const arma::vec& Y,
                                        arma::vec& Xbeta,
                                        const arma::mat& X,
                                        const double sigma,
                                        arma::vec& beta,
                                             arma::ivec& betaind) {

    betaind -= 1;


    const arma::uword p = X.n_cols;

    arma::vec beta_new(p, fill::zeros);
    arma::vec acc_vec(p, fill::zeros);
    arma::vec residual  = Y - Xbeta;

    double log_lik        =   normalloglik(residual, sigma);

    for (uword i = 0; i < p; ++i) {

        arma::vec  Xi = X.col(i);
        double U_  = R::runif(-0.5+1e-8, p-0.5+1e-8);

        //
        int   proposal = (int) std::round(U_);

        arma::vec  Xj = X.col(proposal);
        //log_lik_star
        arma::vec  Xbeta_star =  Xbeta + Xi * (beta(proposal) - beta(i)) - Xj * (beta(proposal) - beta(i));

        arma::vec residual_star     = Y - Xbeta_star;
        double log_lik_star     =   normalloglik(residual_star, sigma);
        double log_U = std::log(R::runif(0,1));
        double MH_ratio = log_lik_star - log_lik ;
        if(log_U  < MH_ratio){
            double beta_swap_star =beta(proposal);
            int beta_ind_prop = betaind(proposal);
            beta(proposal)  = beta(i);

            betaind(proposal) = betaind(i);
            betaind(i) = beta_ind_prop;
            beta(i)     = beta_swap_star;
            log_lik     = log_lik_star;
            residual    = residual_star;
            Xbeta = Xbeta_star;
            acc_vec(i) += 1.;
        }
    }


    betaind += 1;
    return List::create(Named("acc.vec") = wrap(acc_vec),
                        Named("beta.ind") = wrap(betaind));
}


// [[Rcpp::export]]
Rcpp::List permute_gibbs_beta_logistic_cpp(const arma::vec& Y,
                                         arma::vec& Xbeta,
                                         const arma::mat& X,
                                         arma::vec& beta) {



    const arma::uword p = X.n_cols;

    arma::vec beta_new(p, fill::zeros);
    arma::vec acc_vec(p, fill::zeros);
    arma::vec expXbeta    = exp(Xbeta);
    double log_lik        =  binomaloglik_(Y, Xbeta, expXbeta);

    for (uword i = 0; i < p; ++i) {

        arma::vec  Xi = X.col(i);
        double U_  = R::runif(-0.5+1e-8, p-0.5+1e-8);

        //
        int   proposal = (int) std::round(U_);

        arma::vec  Xj = X.col(proposal);
        //log_lik_star
        arma::vec  Xbeta_star =  Xbeta + Xi * (beta(proposal) - beta(i)) - Xj * (beta(proposal) - beta(i));


        arma::vec expXbeta_star    = exp(Xbeta_star);
        double log_lik_star     = binomaloglik_(Y, Xbeta_star, expXbeta_star);
        double log_U = std::log(R::runif(0,1));
        double MH_ratio = log_lik_star - log_lik ;
        if(log_U  < MH_ratio){
            double beta_swap_star =beta(proposal);
            beta(proposal)  = beta(i);

            beta(i)     = beta_swap_star;
            log_lik     = log_lik_star;
            Xbeta = Xbeta_star;
            acc_vec(i) += 1.;
        }
    }


    return List::create(Named("acc.vec") = wrap(acc_vec));
}



// Binary search helper function for an Armadillo row vector.
int binary_search_sample(const arma::rowvec& cum_row, double u) {
    int lo = 0;
    int hi = cum_row.n_elem - 1;

    while (lo < hi) {
        int mid = (lo + hi) / 2;
        if (cum_row(mid) < u) {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }
    return lo;
}



// [[Rcpp::export]]
int sample_j_cumulative(const arma::mat& cumP, arma::uword i) {
    double u = R::runif(0.0, 1.0);
    arma::rowvec cum_row = cumP.row(i);
    int sampled_index = binary_search_sample(cum_row, u);

    return sampled_index ;
}


// [[Rcpp::export]]
Rcpp::List permute_gibbs_beta_normal_corr_cpp(const arma::vec& Y,
                                              arma::vec& Xbeta,
                                              const arma::mat& X,
                                              const double sigma,
                                              arma::vec& beta,
                                              arma::ivec& betaind,
                                              const arma::mat& cumP) {

    betaind -= 1;


    const arma::uword p = X.n_cols;

    arma::vec beta_new(p, fill::zeros);
    arma::vec acc_vec(p, fill::zeros);
    arma::vec residual  = Y - Xbeta;

    double log_lik        =   normalloglik(residual, sigma);

    for (uword i = 0; i < p; ++i) {

        arma::vec  Xi = X.col(i);

        //
        int    proposal = sample_j_cumulative(cumP, i) ;
        arma::vec  Xj = X.col(proposal);
        //log_lik_star
        arma::vec  Xbeta_star =  Xbeta + Xi * (beta(proposal) - beta(i)) - Xj * (beta(proposal) - beta(i));

        arma::vec residual_star     = Y - Xbeta_star;
        double log_lik_star     =   normalloglik(residual_star, sigma);
        double log_U = std::log(R::runif(0,1));
        double MH_ratio = log_lik_star - log_lik ;
        if(log_U  < MH_ratio){
            double beta_swap_star =beta(proposal);
            int beta_ind_prop = betaind(proposal);
            beta(proposal)  = beta(i);

            betaind(proposal) = betaind(i);
            betaind(i) = beta_ind_prop;
            beta(i)     = beta_swap_star;
            log_lik     = log_lik_star;
            residual    = residual_star;
            Xbeta = Xbeta_star;
            acc_vec(i) += 1.;
        }
    }


    betaind += 1;
    return List::create(Named("acc.vec") = wrap(acc_vec),
                        Named("beta.ind") = wrap(betaind));
}
